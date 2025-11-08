#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ckks_polynomial.h"
#include "ckks_encode.h"
#include "ckks_crypto.h"

int test_polynomial(){

 const size_t degree =4098; /* Dimension of the cyclotomic polynomial x^1024 + 1 */
    const uint64_t modulus = 65537; /* 16-bit Fermat prime used for the demo */

    CKKSPolynomial a;
    CKKSPolynomial b;
    CKKSPolynomial result;

    assert(CKKSPolynomialInit(&a, degree, modulus) == 0);
    assert(CKKSPolynomialInit(&b, degree, modulus) == 0);
    assert(CKKSPolynomialInit(&result, degree, modulus) == 0);

    uint64_t *coeffsA = (uint64_t *)malloc(degree * sizeof(uint64_t));
    uint64_t *coeffsB = (uint64_t *)malloc(degree * sizeof(uint64_t));
    uint64_t *expectedAdd = (uint64_t *)malloc(degree * sizeof(uint64_t));
    uint64_t *expectedSub = (uint64_t *)malloc(degree * sizeof(uint64_t));
    uint64_t *expectedMul = (uint64_t *)malloc(degree * sizeof(uint64_t));
    uint64_t *expectedDiv = (uint64_t *)malloc(degree * sizeof(uint64_t));
    assert(coeffsA != NULL && coeffsB != NULL && expectedAdd != NULL && expectedSub != NULL && expectedMul != NULL &&
           expectedDiv != NULL);

    for (size_t i = 0; i < degree; ++i) {
        uint64_t ai = (17ULL * (uint64_t)i + 19ULL) % modulus;
        uint64_t bi = (23ULL * (uint64_t)(i + 1ULL) * (uint64_t)(i + 1ULL) + 29ULL) % modulus;
        coeffsA[i] = ai;
        coeffsB[i] = bi;
        uint64_t add = ai + bi;
        if (add >= modulus) {
            add -= modulus;
        }
        expectedAdd[i] = add;
        expectedSub[i] = CKKSModReduceInt64((int64_t)ai - (int64_t)bi, modulus);
    }

    assert(CKKSPolynomialSet(&a, coeffsA, degree) == 0);
    assert(CKKSPolynomialSet(&b, coeffsB, degree) == 0);

    CKKSNaiveMulReduce(coeffsA, coeffsB, degree, modulus, expectedMul);

    uint64_t inverseThree = CKKSPolynomialModInverse(3, modulus);
    assert(inverseThree != 0);
    for (size_t i = 0; i < degree; ++i) {
        expectedDiv[i] = CKKSPolynomialMulMod(coeffsA[i], inverseThree, modulus);
    }

    /* Addition */
    assert(CKKSPolynomialAdd(&a, &b, &result) == 0);
    assert(CKKSMatchCoefficients(&result, expectedAdd) == 1);
    printf("CKKS polynomial addition passed for degree %zu.\n", degree);

    /* Subtraction */
    assert(CKKSPolynomialSub(&a, &b, &result) == 0);
    assert(CKKSMatchCoefficients(&result, expectedSub) == 1);
    printf("CKKS polynomial subtraction passed for degree %zu.\n", degree);

    /* Multiplication modulo x^degree + 1 */
    assert(CKKSMul(&a, &b, &result) == 0);
    assert(CKKSMatchCoefficients(&result, expectedMul) == 1);
    printf("CKKS polynomial multiplication passed for degree %zu.\n", degree);

    /* Division by the scalar 3 */
    assert(CKKSScalarDivide(&a, 3, &result) == 0);
    assert(CKKSMatchCoefficients(&result, expectedDiv) == 1);
    printf("CKKS polynomial scalar division passed for degree %zu.\n", degree);

    const double plaintextSamples[] = {1.125, -0.75, 3.5, -2.125, 0.625, 4.0, -1.5, 2.25};
    const size_t plaintextCount = sizeof(plaintextSamples) / sizeof(plaintextSamples[0]);
    size_t slotCapacity = degree / 2U;
    CKKSComplex *slots = (CKKSComplex *)calloc(slotCapacity, sizeof(CKKSComplex));
    uint64_t *encodedCoeffs = (uint64_t *)calloc(degree, sizeof(uint64_t));
    assert(slots != NULL && encodedCoeffs != NULL);

    CKKSCopyFloatsToComplexSlots(plaintextSamples, plaintextCount, slots, slotCapacity);
    printf("Converted floating-point inputs into CKKS complex slots (first 6 entries):\n");
    for (size_t i = 0; i < 6 && i < plaintextCount; ++i) {
        printf("  slot[%zu] = %.6f + %.6fi\n", i, slots[i].real, slots[i].imag);
    }

    const double scalingFactor = 4096.0;
    CKKSEncodeComplexSlotsToPolynomial(slots, plaintextCount, degree, scalingFactor, modulus, encodedCoeffs);

    CKKSPolynomial encodedPoly;
    assert(CKKSPolynomialInit(&encodedPoly, degree, modulus) == 0);
    assert(CKKSPolynomialSet(&encodedPoly, encodedCoeffs, degree) == 0);

    printf("Encoded CKKS polynomial coefficients derived from the complex slots (first 8 entries):\n");
    for (size_t i = 0; i < 8 && i < degree; ++i) {
        printf("  coeff[%zu] = %llu\n", i, (unsigned long long)encodedPoly.coeffs[i]);
    }

    CKKSComplex *decodedSlots = (CKKSComplex *)calloc(slotCapacity, sizeof(CKKSComplex));
    double *recoveredFloatsFromComplex = (double *)calloc(plaintextCount, sizeof(double));
    double *recoveredFloatsDirect = (double *)calloc(plaintextCount, sizeof(double));
    assert(decodedSlots != NULL && recoveredFloatsFromComplex != NULL && recoveredFloatsDirect != NULL);

    CKKSDecodePolynomialToComplex(encodedPoly.coeffs, degree, scalingFactor, modulus, decodedSlots, plaintextCount);
    printf("Decoded CKKS complex slots recovered from the polynomial (first 6 entries):\n");
    for (size_t i = 0; i < 6 && i < plaintextCount; ++i) {
        printf("  slot[%zu] = %.6f + %.6fi\n", i, decodedSlots[i].real, decodedSlots[i].imag);
    }

    CKKSComplexSlotsToFloats(decodedSlots, plaintextCount, recoveredFloatsFromComplex, plaintextCount);
    CKKSDecodePolynomialToFloats(encodedPoly.coeffs, degree, scalingFactor, modulus, recoveredFloatsDirect,
                                 plaintextCount);

    const double comparisonTolerance = 5e-3;
    double maxDiffComplex = 0.0;
    double maxDiffDirect = 0.0;
    for (size_t i = 0; i < plaintextCount; ++i) {
        double diffComplex = fabs(plaintextSamples[i] - recoveredFloatsFromComplex[i]);
        double diffDirect = fabs(plaintextSamples[i] - recoveredFloatsDirect[i]);
        if (diffComplex > maxDiffComplex) {
            maxDiffComplex = diffComplex;
        }
        if (diffDirect > maxDiffDirect) {
            maxDiffDirect = diffDirect;
        }
    }

    printf("Maximum reconstruction error via complex slots: %.10f\n", maxDiffComplex);
    printf("Maximum reconstruction error via direct decoding: %.10f\n", maxDiffDirect);

    assert(maxDiffComplex <= comparisonTolerance);
    assert(maxDiffDirect <= comparisonTolerance);

    printf("Recovered floating-point slots via complex projection (first 6 entries):\n");
    for (size_t i = 0; i < 6 && i < plaintextCount; ++i) {
        printf("  value[%zu] = %.6f\n", i, recoveredFloatsFromComplex[i]);
    }

    printf("Recovered floating-point slots via direct polynomial decoding (first 6 entries):\n");
    for (size_t i = 0; i < 6 && i < plaintextCount; ++i) {
        printf("  value[%zu] = %.6f\n", i, recoveredFloatsDirect[i]);
    }

    printf("✅ CKKS float-to-polynomial round-trip passed for %zu inputs.\n", plaintextCount);

    printf("\n=== CKKS conversion sanity sample ===\n");
    const double conversionDemo[] = {0.5, -1.25, 2.75, -3.5, 4.0};
    const size_t conversionDemoCount = sizeof(conversionDemo) / sizeof(conversionDemo[0]);
    const double conversionTolerance = 5e-3;

    memset(slots, 0, slotCapacity * sizeof(CKKSComplex));
    memset(encodedCoeffs, 0, degree * sizeof(uint64_t));
    memset(decodedSlots, 0, slotCapacity * sizeof(CKKSComplex));
    memset(recoveredFloatsFromComplex, 0, plaintextCount * sizeof(double));
    memset(recoveredFloatsDirect, 0, plaintextCount * sizeof(double));

    CKKSCopyFloatsToComplexSlots(conversionDemo, conversionDemoCount, slots, slotCapacity);
    CKKSEncodeComplexSlotsToPolynomial(slots, conversionDemoCount, degree, scalingFactor, modulus, encodedCoeffs);
    assert(CKKSPolynomialSet(&encodedPoly, encodedCoeffs, degree) == 0);

    CKKSDecodePolynomialToComplex(encodedPoly.coeffs, degree, scalingFactor, modulus, decodedSlots,
                                  conversionDemoCount);
    CKKSComplexSlotsToFloats(decodedSlots, conversionDemoCount, recoveredFloatsFromComplex, conversionDemoCount);
    CKKSDecodePolynomialToFloats(encodedPoly.coeffs, degree, scalingFactor, modulus, recoveredFloatsDirect,
                                 conversionDemoCount);

    double maxConversionDiffComplex = 0.0;
    double maxConversionDiffDirect = 0.0;
    for (size_t i = 0; i < conversionDemoCount; ++i) {
        double diffComplex = fabs(conversionDemo[i] - recoveredFloatsFromComplex[i]);
        double diffDirect = fabs(conversionDemo[i] - recoveredFloatsDirect[i]);
        if (diffComplex > maxConversionDiffComplex) {
            maxConversionDiffComplex = diffComplex;
        }
        if (diffDirect > maxConversionDiffDirect) {
            maxConversionDiffDirect = diffDirect;
        }
        printf("  sample[%zu] original = % .6f | via complex = % .6f | via direct = % .6f\n", i, conversionDemo[i],
               recoveredFloatsFromComplex[i], recoveredFloatsDirect[i]);
    }

    printf("Maximum conversion error via complex slots: %.10f\n", maxConversionDiffComplex);
    printf("Maximum conversion error via direct decode: %.10f\n", maxConversionDiffDirect);
    assert(maxConversionDiffComplex <= conversionTolerance);
    assert(maxConversionDiffDirect <= conversionTolerance);
    printf("✅ CKKS conversion sanity sample passed for %zu inputs.\n", conversionDemoCount);

    CKKSNaiveMulReduce(coeffsA, encodedCoeffs, degree, modulus, expectedMul);
    assert(CKKSMul(&a, &encodedPoly, &result) == 0);
    assert(CKKSMatchCoefficients(&result, expectedMul) == 1);
    printf("CKKS multiplication with encoded polynomial payload passed for degree %zu.\n", degree);

    const CKKSComplex complexSamples[] = {
        {1.125, 0.000},
        {-0.875, 1.250},
        {2.250, -0.625},
        {-1.500, -1.125},
        {0.750, 0.375},
        {3.000, -1.750},
        {-2.250, 0.875},
        {1.125, -0.250}
    };
    const size_t complexSampleCount = sizeof(complexSamples) / sizeof(complexSamples[0]);
    const double complexComparisonTolerance = 1e-2;

    memset(slots, 0, slotCapacity * sizeof(CKKSComplex));
    for (size_t i = 0; i < complexSampleCount && i < slotCapacity; ++i) {
        slots[i] = complexSamples[i];
    }

    memset(encodedCoeffs, 0, degree * sizeof(uint64_t));
    CKKSEncodeComplexSlotsToPolynomial(slots, complexSampleCount, degree, scalingFactor, modulus, encodedCoeffs);
    assert(CKKSPolynomialSet(&encodedPoly, encodedCoeffs, degree) == 0);

    memset(decodedSlots, 0, slotCapacity * sizeof(CKKSComplex));
    CKKSDecodePolynomialToComplex(encodedPoly.coeffs, degree, scalingFactor, modulus, decodedSlots,
                                  complexSampleCount);

    double maxComplexRealDiff = 0.0;
    double maxComplexImagDiff = 0.0;
    for (size_t i = 0; i < complexSampleCount; ++i) {
        double diffReal = fabs(slots[i].real - decodedSlots[i].real);
        double diffImag = fabs(slots[i].imag - decodedSlots[i].imag);
        if (diffReal > maxComplexRealDiff) {
            maxComplexRealDiff = diffReal;
        }
        if (diffImag > maxComplexImagDiff) {
            maxComplexImagDiff = diffImag;
        }
    }

    double maxTailMagnitude = 0.0;
    for (size_t i = complexSampleCount; i < slotCapacity; ++i) {
        double magnitude = hypot(decodedSlots[i].real, decodedSlots[i].imag);
        if (magnitude > maxTailMagnitude) {
            maxTailMagnitude = magnitude;
        }
    }

    printf("\nCKKS complex slot encode/decode validation (first 6 entries):\n");
    for (size_t i = 0; i < 6 && i < complexSampleCount; ++i) {
        printf("  original[%zu] = %.6f + %.6fi, decoded[%zu] = %.6f + %.6fi\n", i, slots[i].real, slots[i].imag,
               i, decodedSlots[i].real, decodedSlots[i].imag);
    }
    printf("Maximum complex real component error: %.10f\n", maxComplexRealDiff);
    printf("Maximum complex imaginary component error: %.10f\n", maxComplexImagDiff);
    printf("Maximum tail slot magnitude after decoding: %.10f\n", maxTailMagnitude);

    assert(maxComplexRealDiff <= complexComparisonTolerance);
    assert(maxComplexImagDiff <= complexComparisonTolerance);
    assert(maxTailMagnitude <= complexComparisonTolerance);

    memset(recoveredFloatsFromComplex, 0, plaintextCount * sizeof(double));
    memset(recoveredFloatsDirect, 0, plaintextCount * sizeof(double));

    CKKSComplexSlotsToFloats(decodedSlots, complexSampleCount, recoveredFloatsFromComplex, complexSampleCount);
    CKKSDecodePolynomialToFloats(encodedPoly.coeffs, degree, scalingFactor, modulus, recoveredFloatsDirect,
                                 complexSampleCount);

    double maxFloatDiffComplex = 0.0;
    double maxFloatDiffDirect = 0.0;
    for (size_t i = 0; i < complexSampleCount; ++i) {
        double diffComplex = fabs(slots[i].real - recoveredFloatsFromComplex[i]);
        double diffDirect = fabs(slots[i].real - recoveredFloatsDirect[i]);
        if (diffComplex > maxFloatDiffComplex) {
            maxFloatDiffComplex = diffComplex;
        }
        if (diffDirect > maxFloatDiffDirect) {
            maxFloatDiffDirect = diffDirect;
        }
    }

    printf("Maximum float reconstruction error via complex projection: %.10f\n", maxFloatDiffComplex);
    printf("Maximum float reconstruction error via direct decode: %.10f\n", maxFloatDiffDirect);

    assert(maxFloatDiffComplex <= complexComparisonTolerance);
    assert(maxFloatDiffDirect <= complexComparisonTolerance);

    free(decodedSlots);
    free(recoveredFloatsFromComplex);
    free(recoveredFloatsDirect);
    CKKSPolynomialFree(&encodedPoly);
    free(slots);
    free(encodedCoeffs);

    free(coeffsA);
    free(coeffsB);
    free(expectedAdd);
    free(expectedSub);
    free(expectedMul);
    free(expectedDiv);

    CKKSPolynomialFree(&a);
    CKKSPolynomialFree(&b);
    CKKSPolynomialFree(&result);
    return 0;


}





int test_ckks_crypto(void)
{

 const size_t degree = 1024;
    const uint64_t modulus = 65537;
    const double scalingFactor = 4096.0;

    CKKSPrng prng;
    CKKSPrngInit(&prng, 42ULL);

    CKKSKeyPair keyPair;
    if (CKKSKeyPairInit(&keyPair, degree, modulus) != 0) {
        fprintf(stderr, "Failed to initialize key pair.\n");
        return 1;
    }

    if (CKKSKeyGen(&keyPair, &prng) != 0) {
        fprintf(stderr, "Key generation failed.\n");
        CKKSKeyPairFree(&keyPair);
        return 1;
    }
    

    printf("CKKS key generation completed for degree %zu.\n", degree);

    const double plaintextValues[] = {1.125, -0.75, 3.5, 0.25, -1.0, 2.75, -3.125, 0.875};
    const size_t plaintextCount = sizeof(plaintextValues) / sizeof(plaintextValues[0]);

    CKKSPolynomial plaintext;
    if (CKKSPolynomialInit(&plaintext, degree, modulus) != 0) {
        fprintf(stderr, "Failed to initialize plaintext polynomial.\n");
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    CKKSEncodeFloatsToPolynomial(plaintextValues, plaintextCount, degree, scalingFactor, modulus, plaintext.coeffs);

    printf("Plaintext polynomial (first 8 coefficients):\n");
    CKKSPolynomialPrint(&plaintext, "m", 8);

    CKKSCiphertext ciphertext;
    if (CKKSCiphertextInit(&ciphertext, degree, modulus) != 0) {
        fprintf(stderr, "Failed to initialize ciphertext.\n");
        CKKSPolynomialFree(&plaintext);
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    if (CKKSEncrypt(&keyPair, &plaintext, &prng, &ciphertext) != 0) {
        fprintf(stderr, "Encryption failed.\n");
        CKKSCiphertextFree(&ciphertext);
        CKKSPolynomialFree(&plaintext);
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    printf("Encryption produced ciphertext polynomials (first 8 coefficients each):\n");
    CKKSPolynomialPrint(&ciphertext.c0, "c0", 8);
    CKKSPolynomialPrint(&ciphertext.c1, "c1", 8);

    CKKSPolynomial decrypted;
    if (CKKSPolynomialInit(&decrypted, degree, modulus) != 0) {
        fprintf(stderr, "Failed to initialize decrypted polynomial.\n");
        CKKSCiphertextFree(&ciphertext);
        CKKSPolynomialFree(&plaintext);
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    if (CKKSDecrypt(&keyPair, &ciphertext, &decrypted) != 0) {
        fprintf(stderr, "Decryption failed.\n");
        CKKSPolynomialFree(&decrypted);
        CKKSCiphertextFree(&ciphertext);
        CKKSPolynomialFree(&plaintext);
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    printf("Decrypted polynomial (first 8 coefficients):\n");
    CKKSPolynomialPrint(&decrypted, "m'", 8);

    uint64_t halfModulus = modulus / 2U;
    int64_t maxCenteredDiff = 0;
    for (size_t i = 0; i < degree; ++i) {
        int64_t diff = (int64_t)plaintext.coeffs[i] - (int64_t)decrypted.coeffs[i];
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

    printf("Maximum centered coefficient delta after decrypt: %lld\n", (long long)maxCenteredDiff);

    double recovered[plaintextCount];
    double original[plaintextCount];
    memcpy(original, plaintextValues, sizeof(original));

    CKKSDecodePolynomialToFloats(decrypted.coeffs, decrypted.degree, scalingFactor, decrypted.modulus, recovered,
                                 plaintextCount);

    double maxError = 0.0;
    for (size_t i = 0; i < plaintextCount; ++i) {
        double diff = fabs(recovered[i] - original[i]);
        if (diff > maxError) {
            maxError = diff;
        }
        printf("Slot %zu -> original %.6f, recovered %.6f\n", i, original[i], recovered[i]);
    }

    printf("Maximum decoding error after encrypt/decrypt: %.10f\n", maxError);
    int success = (fabsl((long double)maxCenteredDiff) <= 64 && maxError <= 5e-3) ? 1 : 0;
    if (success) {
        printf("✅ Decryption maintained plaintext within acceptable CKKS noise bounds.\n");
    } else {
        fprintf(stderr, "⚠️  Decryption noise exceeded expected bounds.\n");
    }

    CKKSPolynomialFree(&decrypted);
    CKKSCiphertextFree(&ciphertext);
    CKKSPolynomialFree(&plaintext);
    CKKSKeyPairFree(&keyPair);
    return success ? 0 : 1;

}

int main(void)
{
  test_polynomial();
  test_ckks_crypto();
 return 0;
}
