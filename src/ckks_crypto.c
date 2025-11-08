#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ckks_polynomial.h"  // For CKKSModReduceInt64
#include "ckks_encode.h"  // Include the encode header
#include "ckks_crypto.h"  // For CKKSKeyPair, CKKSCiphertext


void CKKSPrngInit(CKKSPrng *prng, uint64_t seed)
{
    if (prng == NULL) {
        return;
    }
    prng->state = seed != 0 ? seed : 0xDEADBEEFCAFEBABEULL;
}

 uint64_t CKKSPrngNext(CKKSPrng *prng)
{
    uint64_t x = prng->state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    prng->state = x;
    return x * 0x2545F4914F6CDD1DULL;
}

 uint64_t CKKSPrngNextMod(CKKSPrng *prng, uint64_t modulus)
{
    return modulus == 0 ? 0 : CKKSPrngNext(prng) % modulus;
}

static int64_t CKKSPrngNextCenteredTernary(CKKSPrng *prng)
{
    uint64_t sample = CKKSPrngNext(prng) % 3ULL;
    if (sample == 0ULL) {
        return 0;
    }
    return sample == 1ULL ? 1 : -1;
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

 int CKKSKeyGen(CKKSKeyPair *keyPair, CKKSPrng *prng)
{
    if (keyPair == NULL || prng == NULL) {
        return -1;
    }

    size_t degree = keyPair->secret.degree;
    uint64_t modulus = keyPair->secret.modulus;

    for (size_t i = 0; i < degree; ++i) {
        int64_t coeff = CKKSPrngNextCenteredTernary(prng);
        keyPair->secret.coeffs[i] = CKKSModReduceInt64(coeff, modulus);
    }

    for (size_t i = 0; i < degree; ++i) {
        keyPair->publicA.coeffs[i] = CKKSPrngNextMod(prng, modulus);
    }

    CKKSPolynomial product;
    if (CKKSPolynomialInit(&product, degree, modulus) != 0) {
        return -1;
    }

    if (CKKSMul(&keyPair->publicA, &keyPair->secret, &product) != 0) {
        CKKSPolynomialFree(&product);
        return -1;
    }

    for (size_t i = 0; i < degree; ++i) {
        int64_t error = CKKSPrngNextCenteredTernary(prng);
        int64_t value = -((int64_t)product.coeffs[i]) + error;
        keyPair->publicB.coeffs[i] = CKKSModReduceInt64(value, modulus);
    }

    CKKSPolynomialFree(&product);
    return 0;
}

 int CKKSEncrypt(const CKKSKeyPair *keyPair, const CKKSPolynomial *plaintext, CKKSPrng *prng,
                       CKKSCiphertext *ciphertext)
{
    if (keyPair == NULL || plaintext == NULL || prng == NULL || ciphertext == NULL) {
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
        u.coeffs[i] = CKKSModReduceInt64(CKKSPrngNextCenteredTernary(prng), modulus);
        e1.coeffs[i] = CKKSModReduceInt64(CKKSPrngNextCenteredTernary(prng), modulus);
        e2.coeffs[i] = CKKSModReduceInt64(CKKSPrngNextCenteredTernary(prng), modulus);
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

    if (CKKSPolynomialAdd(&ciphertext->c0, &temp, decoded) != 0) {
        CKKSPolynomialFree(&temp);
        return -1;
    }

    CKKSPolynomialFree(&temp);
    return 0;
}
