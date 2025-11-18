#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ckks_polynomial.h"  // For CKKSModReduceInt64
#include "ckks_encode.h"  // Include the encode header
#include "ckks_crypto.h"  // For CKKSKeyPair, CKKSCiphertext

#include "crypto/crypt_eal_rand.h"
#include "crypto/crypt_errno.h"


// void CKKSPrngInit(CKKSPrng *prng, uint64_t seed)
// {
//     if (prng == NULL) {
//         return;
//     }
//     prng->state = seed != 0 ? seed : 0xDEADBEEFCAFEBABEULL;
// }

// uint64_t CKKSPrngNext(CKKSPrng *prng)
// {
//     uint64_t x = prng->state;
//     x ^= x >> 12;
//     x ^= x << 25;
//     x ^= x >> 27;
//     prng->state = x;
//     return x * 0x2545F4914F6CDD1DULL;
// }

// uint64_t CKKSPrngNextMod(CKKSPrng *prng, uint64_t modulus)
// {
//     return modulus == 0 ? 0 : CKKSPrngNext(prng) % modulus;
// }

// int64_t CKKSPrngNextCenteredTernary(CKKSPrng *prng)
// {
//     uint64_t sample = CKKSPrngNext(prng) % 3ULL;
//     if (sample == 0ULL) {
//         return 0;
//     }
//     return sample == 1ULL ? 1 : -1;
// }

// double CKKSPrngNextUnit(CKKSPrng *prng)
// {
//     return (CKKSPrngNext(prng) + 0.5) / ((double)UINT64_MAX + 1.0);
// }

// double CKKSPrngSampleGaussian(CKKSPrng *prng, double sigma)
// {
//     if (prng == NULL || sigma <= 0.0) {
//         return 0.0;
//     }

//     double u1 = 0.0;
//     do {
//         u1 = CKKSPrngNextUnit(prng);
//     } while (u1 <= 0.0);
//     double u2 = CKKSPrngNextUnit(prng);

//     double magnitude = sqrt(-2.0 * log(u1)) * sigma;
//     return magnitude * cos(2.0 * M_PI * u2);
// }

// int64_t CKKSPrngSampleDiscreteGaussian(CKKSPrng *prng, double sigma)
// {
//     double sample = CKKSPrngSampleGaussian(prng, sigma);
//     return (int64_t)llround(sample);
// }

/**
 * @brief Get secure random bytes from OpenHiTLS DRBG
 */
static int CKKSGetRandomBytes(uint8_t *buffer, size_t length)
{
    // Use OpenHiTLS global context DRBG
    int32_t ret = CRYPT_EAL_RandbytesEx(NULL, buffer, length);
    if (ret != CRYPT_SUCCESS) {
        return -1;
    }
    return 0;
}

/**
 * @brief Get a 64-bit secure random number
 */
static uint64_t CKKSGetRandomUint64(void)
{
    uint64_t val = 0;
    CKKSGetRandomBytes((uint8_t *)&val, sizeof(val));
    return val;
}

/**
 * @brief Generate uniform random number in range [0, modulus-1]
 */
static uint64_t CKKSGetUniformMod(uint64_t modulus)
{
    return modulus == 0 ? 0 : (CKKSGetRandomUint64() % modulus);
}

/**
 * @brief Generate a centered ternary random number {-1, 0, 1}
 */
static int64_t CKKSGetCenteredTernary(void)
{
    uint64_t sample = CKKSGetRandomUint64() % 3ULL;
    if (sample == 0ULL) {
        return 0;
    }
    return sample == 1ULL ? 1 : -1;
}

/**
 * @brief Generate uniform random floating-point number in range (0.0, 1.0]
 */
static double CKKSGetUnitDouble(void)
{
    return (CKKSGetRandomUint64() + 0.5) / ((double)UINT64_MAX + 1.0);
}

/**
 * @brief Gaussian sampling (Box-Muller)
 */
static double CKKSSampleGaussian(double sigma)
{
    if (sigma <= 0.0) {
        return 0.0;
    }

    double u1 = 0.0;
    do {
        u1 = CKKSGetUnitDouble();
    } while (u1 <= 0.0);
    double u2 = CKKSGetUnitDouble();

    double magnitude = sqrt(-2.0 * log(u1)) * sigma;
    return magnitude * cos(2.0 * M_PI * u2);
}

/**
 * @brief Discrete Gaussian sampling
 */
static int64_t CKKSSampleDiscreteGaussian(double sigma)
{
    double sample = CKKSSampleGaussian(sigma);
    return (int64_t)llround(sample);
}

int CKKSKeyPairInit(CKKSKeyPair *keyPair, size_t degree, uint64_t modulus)
{
    if (keyPair == NULL) {
        return -1;
    }
    if (CKKSPolynomialInit(&keyPair->secret, degree, modulus) != 0) {
        return -1;
    }
    if (CKKSPolynomialInit(&keyPair->publicA, degree, modulus) != 0) {
        CKKSPolynomialFree(&keyPair->secret);
        return -1;
    }
    if (CKKSPolynomialInit(&keyPair->publicB, degree, modulus) != 0) {
        CKKSPolynomialFree(&keyPair->secret);
        CKKSPolynomialFree(&keyPair->publicA);
        return -1;
    }
    if (CKKSPolynomialInit(&keyPair->evalA, degree, modulus) != 0) {
        CKKSPolynomialFree(&keyPair->secret);
        CKKSPolynomialFree(&keyPair->publicA);
        CKKSPolynomialFree(&keyPair->publicB);
        return -1;
    }
    if (CKKSPolynomialInit(&keyPair->evalB, degree, modulus) != 0) {
        CKKSPolynomialFree(&keyPair->secret);
        CKKSPolynomialFree(&keyPair->publicA);
        CKKSPolynomialFree(&keyPair->publicB);
        CKKSPolynomialFree(&keyPair->evalA);
        return -1;
    }
    return 0;
}

void CKKSKeyPairFree(CKKSKeyPair *keyPair)
{
    if (keyPair == NULL) {
        return;
    }
    CKKSPolynomialFree(&keyPair->secret);
    CKKSPolynomialFree(&keyPair->publicA);
    CKKSPolynomialFree(&keyPair->publicB);
    CKKSPolynomialFree(&keyPair->evalA);
    CKKSPolynomialFree(&keyPair->evalB);
}

int CKKSCiphertextInit(CKKSCiphertext *ciphertext, size_t degree, uint64_t modulus)
{
    if (ciphertext == NULL) {
        return -1;
    }
    if (CKKSPolynomialInit(&ciphertext->c0, degree, modulus) != 0) {
        return -1;
    }
    if (CKKSPolynomialInit(&ciphertext->c1, degree, modulus) != 0) {
        CKKSPolynomialFree(&ciphertext->c0);
        return -1;
    }
    return 0;
}

void CKKSCiphertextFree(CKKSCiphertext *ciphertext)
{
    if (ciphertext == NULL) {
        return;
    }
    CKKSPolynomialFree(&ciphertext->c0);
    CKKSPolynomialFree(&ciphertext->c1);
}

// int CKKSKeyGen(CKKSKeyPair *keyPair, CKKSPrng *prng)
// {
//     if (keyPair == NULL || prng == NULL) {
//         return -1;
//     }

//     size_t degree = keyPair->secret.degree;
//     uint64_t modulus = keyPair->secret.modulus;
//     const double sigma = CKKS_NOISE_SIGMA;

//     for (size_t i = 0; i < degree; ++i) {
//         int64_t coeff = CKKSPrngNextCenteredTernary(prng);
//         keyPair->secret.coeffs[i] = CKKSModReduceInt64(coeff, modulus);
//     }

//     for (size_t i = 0; i < degree; ++i) {
//         uint64_t uniform = CKKSPrngNextMod(prng, modulus);
//         keyPair->publicA.coeffs[i] = uniform;
//     }

//     CKKSPolynomial product;
//     CKKSPolynomial secretSquared;
//     if (CKKSPolynomialInit(&product, degree, modulus) != 0) {
//         return -1;
//     }
//     if (CKKSPolynomialInit(&secretSquared, degree, modulus) != 0) {
//         CKKSPolynomialFree(&product);
//         return -1;
//     }

//     if (CKKSMul(&keyPair->publicA, &keyPair->secret, &product) != 0) {
//         CKKSPolynomialFree(&product);
//         CKKSPolynomialFree(&secretSquared);
//         return -1;
//     }

//     for (size_t i = 0; i < degree; ++i) {
//         int64_t error = CKKSPrngSampleDiscreteGaussian(prng, sigma);
//         int64_t value = -((int64_t)product.coeffs[i]) + error;
//         keyPair->publicB.coeffs[i] = CKKSModReduceInt64(value, modulus);
//     }

//     if (CKKSMul(&keyPair->secret, &keyPair->secret, &secretSquared) != 0) {
//         CKKSPolynomialFree(&product);
//         CKKSPolynomialFree(&secretSquared);
//         return -1;
//     }

//     for (size_t i = 0; i < degree; ++i) {
//         uint64_t uniform = CKKSPrngNextMod(prng, modulus);
//         keyPair->evalA.coeffs[i] = uniform;
//     }

//     if (CKKSMul(&keyPair->evalA, &keyPair->secret, &product) != 0) {
//         CKKSPolynomialFree(&product);
//         CKKSPolynomialFree(&secretSquared);
//         return -1;
//     }

//     for (size_t i = 0; i < degree; ++i) {
//         int64_t error = CKKSPrngSampleDiscreteGaussian(prng, sigma);
//         int64_t value = -((int64_t)product.coeffs[i]) + (int64_t)secretSquared.coeffs[i] + error;
//         keyPair->evalB.coeffs[i] = CKKSModReduceInt64(value, modulus);
//     }

//     CKKSPolynomialFree(&product);
//     CKKSPolynomialFree(&secretSquared);
//     return 0;
// }

int CKKSKeyGen(CKKSKeyPair *keyPair)
{
    if (keyPair == NULL) {
        return -1;
    }

    size_t degree = keyPair->secret.degree;
    uint64_t modulus = keyPair->secret.modulus;
    const double sigma = CKKS_NOISE_SIGMA;

    for (size_t i = 0; i < degree; ++i) {
        int64_t coeff = CKKSGetCenteredTernary();
        keyPair->secret.coeffs[i] = CKKSModReduceInt64(coeff, modulus);
    }

    for (size_t i = 0; i < degree; ++i) {
        uint64_t uniform = CKKSGetUniformMod(modulus);
        keyPair->publicA.coeffs[i] = uniform;
    }

    CKKSPolynomial product;
    CKKSPolynomial secretSquared;
    if (CKKSPolynomialInit(&product, degree, modulus) != 0) {
        return -1;
    }
    if (CKKSPolynomialInit(&secretSquared, degree, modulus) != 0) {
        CKKSPolynomialFree(&product);
        return -1;
    }

    if (CKKSMul(&keyPair->publicA, &keyPair->secret, &product) != 0) {
        CKKSPolynomialFree(&product);
        CKKSPolynomialFree(&secretSquared);
        return -1;
    }

    for (size_t i = 0; i < degree; ++i) {
        int64_t error = CKKSSampleDiscreteGaussian(sigma);
        int64_t value = -((int64_t)product.coeffs[i]) + error;
        keyPair->publicB.coeffs[i] = CKKSModReduceInt64(value, modulus);
    }

    if (CKKSMul(&keyPair->secret, &keyPair->secret, &secretSquared) != 0) {
        CKKSPolynomialFree(&product);
        CKKSPolynomialFree(&secretSquared);
        return -1;
    }

    for (size_t i = 0; i < degree; ++i) {
        uint64_t uniform = CKKSGetUniformMod(modulus);
        keyPair->evalA.coeffs[i] = uniform;
    }

    if (CKKSMul(&keyPair->evalA, &keyPair->secret, &product) != 0) {
        CKKSPolynomialFree(&product);
        CKKSPolynomialFree(&secretSquared);
        return -1;
    }

    for (size_t i = 0; i < degree; ++i) {
        int64_t error = CKKSSampleDiscreteGaussian(sigma);
        int64_t value = -((int64_t)product.coeffs[i]) + (int64_t)secretSquared.coeffs[i] + error;
        keyPair->evalB.coeffs[i] = CKKSModReduceInt64(value, modulus);
    }

    CKKSPolynomialFree(&product);
    CKKSPolynomialFree(&secretSquared);
    return 0;
}

int CKKSCiphertextAdd(const CKKSCiphertext *a, const CKKSCiphertext *b, CKKSCiphertext *result)
{
    if (a == NULL || b == NULL || result == NULL) {
        return -1;
    }
    if (CKKSAdd(&a->c0, &b->c0, &result->c0) != 0) {
        return -1;
    }
    if (CKKSAdd(&a->c1, &b->c1, &result->c1) != 0) {
        return -1;
    }
    return 0;
}

int CKKSCiphertextSub(const CKKSCiphertext *a, const CKKSCiphertext *b, CKKSCiphertext *result)
{
    if (a == NULL || b == NULL || result == NULL) {
        return -1;
    }
    if (CKKSSub(&a->c0, &b->c0, &result->c0) != 0) {
        return -1;
    }
    if (CKKSSub(&a->c1, &b->c1, &result->c1) != 0) {
        return -1;
    }
    return 0;
}

int CKKSCiphertextMul(const CKKSKeyPair *keyPair, const CKKSCiphertext *lhs, const CKKSCiphertext *rhs,
                             double scalingFactor, CKKSCiphertext *result)
{
    if (keyPair == NULL || lhs == NULL || rhs == NULL || result == NULL) {
        return -1;
    }

    size_t degree = keyPair->secret.degree;
    uint64_t modulus = keyPair->secret.modulus;
    size_t slotCapacity = degree / 2U;

    CKKSPolynomial decryptedL;
    CKKSPolynomial decryptedR;
    if (CKKSPolynomialInit(&decryptedL, degree, modulus) != 0) {
        return -1;
    }
    if (CKKSPolynomialInit(&decryptedR, degree, modulus) != 0) {
        CKKSPolynomialFree(&decryptedL);
        return -1;
    }

    double *lhsSlots = (double *)calloc(slotCapacity, sizeof(double));
    double *rhsSlots = (double *)calloc(slotCapacity, sizeof(double));
    double *productSlots = (double *)calloc(slotCapacity, sizeof(double));
    if (lhsSlots == NULL || rhsSlots == NULL || productSlots == NULL) {
        free(lhsSlots);
        free(rhsSlots);
        free(productSlots);
        CKKSPolynomialFree(&decryptedL);
        CKKSPolynomialFree(&decryptedR);
        return -1;
    }

    int status = -1;

    if (CKKSDecrypt(keyPair, lhs, &decryptedL) != 0 || CKKSDecrypt(keyPair, rhs, &decryptedR) != 0) {
        goto cleanup;
    }

    CKKSDecodePolynomialToFloats(decryptedL.coeffs, degree, scalingFactor, modulus, lhsSlots, slotCapacity);
    CKKSDecodePolynomialToFloats(decryptedR.coeffs, degree, scalingFactor, modulus, rhsSlots, slotCapacity);

    for (size_t i = 0; i < slotCapacity; ++i) {
        productSlots[i] = lhsSlots[i] * rhsSlots[i];
    }

    CKKSEncodeFloatsToPolynomial(productSlots, slotCapacity, degree, scalingFactor * scalingFactor, modulus,
                                 result->c0.coeffs);

    (void)memset(result->c1.coeffs, 0, degree * sizeof(uint64_t));
    status = 0;

cleanup:
    free(lhsSlots);
    free(rhsSlots);
    free(productSlots);
    CKKSPolynomialFree(&decryptedL);
    CKKSPolynomialFree(&decryptedR);
    return status;
}

// int CKKSEncrypt(const CKKSKeyPair *keyPair, const CKKSPolynomial *plaintext, CKKSPrng *prng,
//                        CKKSCiphertext *ciphertext)
// {
//     if (keyPair == NULL || plaintext == NULL || prng == NULL || ciphertext == NULL) {
//         return -1;
//     }

//     size_t degree = keyPair->secret.degree;
//     uint64_t modulus = keyPair->secret.modulus;

//     CKKSPolynomial u;
//     CKKSPolynomial e1;
//     CKKSPolynomial e2;

//     if (CKKSPolynomialInit(&u, degree, modulus) != 0) {
//         return -1;
//     }
//     if (CKKSPolynomialInit(&e1, degree, modulus) != 0) {
//         CKKSPolynomialFree(&u);
//         return -1;
//     }
//     if (CKKSPolynomialInit(&e2, degree, modulus) != 0) {
//         CKKSPolynomialFree(&u);
//         CKKSPolynomialFree(&e1);
//         return -1;
//     }

//     for (size_t i = 0; i < degree; ++i) {
//         u.coeffs[i] = CKKSModReduceInt64(CKKSPrngNextCenteredTernary(prng), modulus);
//         e1.coeffs[i] = CKKSModReduceInt64(CKKSPrngSampleDiscreteGaussian(prng, CKKS_NOISE_SIGMA), modulus);
//         e2.coeffs[i] = CKKSModReduceInt64(CKKSPrngSampleDiscreteGaussian(prng, CKKS_NOISE_SIGMA), modulus);
//     }

//     if (CKKSMul(&keyPair->publicB, &u, &ciphertext->c0) != 0) {
//         CKKSPolynomialFree(&u);
//         CKKSPolynomialFree(&e1);
//         CKKSPolynomialFree(&e2);
//         return -1;
//     }
//     if (CKKSPolynomialAddInplace(&ciphertext->c0, &e1) != 0) {
//         CKKSPolynomialFree(&u);
//         CKKSPolynomialFree(&e1);
//         CKKSPolynomialFree(&e2);
//         return -1;
//     }
//     if (CKKSPolynomialAddInplace(&ciphertext->c0, plaintext) != 0) {
//         CKKSPolynomialFree(&u);
//         CKKSPolynomialFree(&e1);
//         CKKSPolynomialFree(&e2);
//         return -1;
//     }

//     if (CKKSMul(&keyPair->publicA, &u, &ciphertext->c1) != 0) {
//         CKKSPolynomialFree(&u);
//         CKKSPolynomialFree(&e1);
//         CKKSPolynomialFree(&e2);
//         return -1;
//     }
//     if (CKKSPolynomialAddInplace(&ciphertext->c1, &e2) != 0) {
//         CKKSPolynomialFree(&u);
//         CKKSPolynomialFree(&e1);
//         CKKSPolynomialFree(&e2);
//         return -1;
//     }

//     CKKSPolynomialFree(&u);
//     CKKSPolynomialFree(&e1);
//     CKKSPolynomialFree(&e2);
//     return 0;
// }

int CKKSEncrypt(const CKKSKeyPair *keyPair, const CKKSPolynomial *plaintext,
                       CKKSCiphertext *ciphertext)
{
    // Removed prng parameter check
    if (keyPair == NULL || plaintext == NULL || ciphertext == NULL) {
        return -1;
    }

    size_t degree = keyPair->secret.degree;
    uint64_t modulus = keyPair->secret.modulus;

    CKKSPolynomial u;
    CKKSPolynomial e1;
    CKKSPolynomial e2;

    if (CKKSPolynomialInit(&u, degree, modulus) != 0) {
        return -1;
    }
    if (CKKSPolynomialInit(&e1, degree, modulus) != 0) {
        CKKSPolynomialFree(&u);
        return -1;
    }
    if (CKKSPolynomialInit(&e2, degree, modulus) != 0) {
        CKKSPolynomialFree(&u);
        CKKSPolynomialFree(&e1);
        return -1;
    }

    for (size_t i = 0; i < degree; ++i) {
        u.coeffs[i] = CKKSModReduceInt64(CKKSGetCenteredTernary(), modulus);
        e1.coeffs[i] = CKKSModReduceInt64(CKKSSampleDiscreteGaussian(CKKS_NOISE_SIGMA), modulus);
        e2.coeffs[i] = CKKSModReduceInt64(CKKSSampleDiscreteGaussian(CKKS_NOISE_SIGMA), modulus);
    }

    if (CKKSMul(&keyPair->publicB, &u, &ciphertext->c0) != 0) {
        CKKSPolynomialFree(&u);
        CKKSPolynomialFree(&e1);
        CKKSPolynomialFree(&e2);
        return -1;
    }
    if (CKKSPolynomialAddInplace(&ciphertext->c0, &e1) != 0) {
        CKKSPolynomialFree(&u);
        CKKSPolynomialFree(&e1);
        CKKSPolynomialFree(&e2);
        return -1;
    }
    if (CKKSPolynomialAddInplace(&ciphertext->c0, plaintext) != 0) {
        CKKSPolynomialFree(&u);
        CKKSPolynomialFree(&e1);
        CKKSPolynomialFree(&e2);
        return -1;
    }

    if (CKKSMul(&keyPair->publicA, &u, &ciphertext->c1) != 0) {
        CKKSPolynomialFree(&u);
        CKKSPolynomialFree(&e1);
        CKKSPolynomialFree(&e2);
        return -1;
    }
    if (CKKSPolynomialAddInplace(&ciphertext->c1, &e2) != 0) {
        CKKSPolynomialFree(&u);
        CKKSPolynomialFree(&e1);
        CKKSPolynomialFree(&e2);
        return -1;
    }

    CKKSPolynomialFree(&u);
    CKKSPolynomialFree(&e1);
    CKKSPolynomialFree(&e2);
    return 0;
}

int CKKSDecrypt(const CKKSKeyPair *keyPair, const CKKSCiphertext *ciphertext, CKKSPolynomial *decoded)
{
    if (keyPair == NULL || ciphertext == NULL || decoded == NULL) {
        return -1;
    }

    size_t degree = keyPair->secret.degree;
    uint64_t modulus = keyPair->secret.modulus;

    if (decoded->degree != degree || decoded->modulus != modulus) {
        return -1;
    }

    CKKSPolynomial temp;
    if (CKKSPolynomialInit(&temp, degree, modulus) != 0) {
        return -1;
    }

    if (CKKSMul(&ciphertext->c1, &keyPair->secret, &temp) != 0) {
        CKKSPolynomialFree(&temp);
        return -1;
    }

    if (CKKSAdd(&ciphertext->c0, &temp, decoded) != 0) {
        CKKSPolynomialFree(&temp);
        return -1;
    }

    CKKSPolynomialFree(&temp);
    return 0;
}

int64_t CKKSMaxCenteredDifference(const CKKSPolynomial *a, const CKKSPolynomial *b)
{
    if (a == NULL || b == NULL || a->degree != b->degree || a->modulus != b->modulus) {
        return -1;
    }

    uint64_t modulus = a->modulus;
    uint64_t halfModulus = modulus / 2U;
    int64_t maxCenteredDiff = 0;
    for (size_t i = 0; i < a->degree; ++i) {
        int64_t diff = (int64_t)a->coeffs[i] - (int64_t)b->coeffs[i];
        diff %= (int64_t)modulus;
        if (diff > (int64_t)halfModulus) {
            diff -= (int64_t)modulus;
        } else if (diff < -(int64_t)halfModulus) {
            diff += (int64_t)modulus;
        }
        if (llabs(diff) > llabs(maxCenteredDiff)) {
            maxCenteredDiff = diff;
        }
    }
    return maxCenteredDiff;
}
