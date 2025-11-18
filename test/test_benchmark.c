/*
 * CKKS Homomorphic Encryption System Performance Test Suite
 * 
 * This file contains all performance benchmark code for evaluating the performance
 * of various operations in the CKKS implementation
 * 
 * Test items include:
 * 1. Polynomial operation performance (addition, subtraction, multiplication, division)
 * 2. Encode/decode performance
 * 3. Encryption/decryption performance
 * 4. Homomorphic operation performance (ciphertext addition, subtraction, multiplication)
 * 
 * Supports multiple data scale tests for evaluating performance at different dimensions
 */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ckks_polynomial.h"
#include "ckks_encode.h"
#include "ckks_crypto.h"

/**
 * Get current time in seconds
 * Uses CLOCK_MONOTONIC to ensure monotonic increase, suitable for performance testing
 */
static double get_time_in_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}

/**
 * Print performance statistics
 * @param operation Operation name
 * @param time_seconds Execution time (seconds)
 * @param iterations Number of iterations
 */
static void print_performance(const char *operation, double time_seconds, size_t iterations) {
    printf("%s: %.6f seconds", operation, time_seconds);
    if (iterations > 1) {
        printf(" (%.2f ops/sec, avg %.6f sec/op)\n", 
               iterations / time_seconds, time_seconds / iterations);
    } else {
        printf("\n");
    }
}

/* ============================================
 * Polynomial Operation Performance Tests
 * ============================================ */

int test_polynomial_benchmark(size_t degree) {
    double start_time, end_time;

    const uint64_t modulus = 65537; /* 16-bit Fermat prime */

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
    assert(coeffsA != NULL && coeffsB != NULL && expectedAdd != NULL && expectedSub != NULL && 
           expectedMul != NULL && expectedDiv != NULL);

    // Initialize test data
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

    /* Polynomial addition performance test */
    const size_t add_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < add_iterations; i++) {
        assert(CKKSAdd(&a, &b, &result) == 0);
    }
    end_time = get_time_in_seconds();
    assert(CKKSMatchCoefficients(&result, expectedAdd) == 1);
    printf("CKKS polynomial addition test passed (degree %zu)\n", degree);
    print_performance("Polynomial addition", end_time - start_time, add_iterations);

    /* Polynomial subtraction performance test */
    const size_t sub_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < sub_iterations; i++) {
        assert(CKKSSub(&a, &b, &result) == 0);
    }
    end_time = get_time_in_seconds();
    assert(CKKSMatchCoefficients(&result, expectedSub) == 1);
    printf("CKKS polynomial subtraction test passed (degree %zu)\n", degree);
    print_performance("Polynomial subtraction", end_time - start_time, sub_iterations);

    /* Polynomial multiplication performance test */
    const size_t mul_iterations = 1000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < mul_iterations; i++) {
        assert(CKKSMul(&a, &b, &result) == 0);
    }
    end_time = get_time_in_seconds();
    assert(CKKSMatchCoefficients(&result, expectedMul) == 1);
    printf("CKKS polynomial multiplication test passed (degree %zu)\n", degree);
    print_performance("Polynomial multiplication", end_time - start_time, mul_iterations);

    /* Polynomial scalar division performance test */
    const size_t div_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < div_iterations; i++) {
        assert(CKKSScalarDivide(&a, 3, &result) == 0);
    }
    end_time = get_time_in_seconds();
    assert(CKKSMatchCoefficients(&result, expectedDiv) == 1);
    printf("CKKS polynomial scalar division test passed (degree %zu)\n", degree);
    print_performance("Polynomial scalar division", end_time - start_time, div_iterations);

    // Clean up resources
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

/* ============================================
 * Encode/Decode Performance Tests
 * ============================================ */

int test_encode_decode_benchmark(size_t degree) {
    double start_time, end_time;

    const uint64_t modulus = 65537;
    const double scalingFactor = 4096.0;

    const double plaintextSamples[] = {1.125, -0.75, 3.5, -2.125, 0.625, 4.0, -1.5, 2.25};
    const size_t plaintextCount = sizeof(plaintextSamples) / sizeof(plaintextSamples[0]);
    size_t slotCapacity = degree / 2U;

    CKKSComplex *slots = (CKKSComplex *)calloc(slotCapacity, sizeof(CKKSComplex));
    uint64_t *encodedCoeffs = (uint64_t *)calloc(degree, sizeof(uint64_t));
    assert(slots != NULL && encodedCoeffs != NULL);

    CKKSCopyFloatsToComplexSlots(plaintextSamples, plaintextCount, slots, slotCapacity);

    /* Encoding performance test: complex slots to polynomial */
    const size_t encode_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < encode_iterations; i++) {
        CKKSEncodeComplexSlotsToPolynomial(slots, plaintextCount, degree, scalingFactor, modulus, encodedCoeffs);
    }
    end_time = get_time_in_seconds();
    print_performance("Encoding (complex slots→polynomial)", end_time - start_time, encode_iterations);

    CKKSPolynomial encodedPoly;
    assert(CKKSPolynomialInit(&encodedPoly, degree, modulus) == 0);
    assert(CKKSPolynomialSet(&encodedPoly, encodedCoeffs, degree) == 0);

    CKKSComplex *decodedSlots = (CKKSComplex *)calloc(slotCapacity, sizeof(CKKSComplex));
    double *recoveredFloats = (double *)calloc(plaintextCount, sizeof(double));
    assert(decodedSlots != NULL && recoveredFloats != NULL);

    /* Decoding performance test: polynomial to complex slots */
    const size_t decode_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < decode_iterations; i++) {
        CKKSDecodePolynomialToComplex(encodedPoly.coeffs, degree, scalingFactor, modulus, 
                                      decodedSlots, plaintextCount);
    }
    end_time = get_time_in_seconds();
    print_performance("Decoding (polynomial→complex slots)", end_time - start_time, decode_iterations);

    /* Direct decoding performance test: polynomial to floating-point */
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < decode_iterations; i++) {
        CKKSDecodePolynomialToFloats(encodedPoly.coeffs, degree, scalingFactor, modulus, 
                                     recoveredFloats, plaintextCount);
    }
    end_time = get_time_in_seconds();
    print_performance("Decoding (polynomial→floating-point)", end_time - start_time, decode_iterations);

    printf("Encode/decode performance test completed\n");

    // Clean up resources
    free(decodedSlots);
    free(recoveredFloats);
    CKKSPolynomialFree(&encodedPoly);
    free(slots);
    free(encodedCoeffs);

    return 0;
}

/* ============================================
 * Encryption/Decryption Performance Tests
 * ============================================ */

int test_encryption_benchmark(size_t degree) {
    double start_time, end_time;

    const uint64_t modulus = 65537;
    const double scalingFactor = 128.0;

    CKKSKeyPair keyPair;
    if (CKKSKeyPairInit(&keyPair, degree, modulus) != 0) {
        fprintf(stderr, "Key pair initialization failed.\n");
        return 1;
    }

    /* Key generation performance test */
    start_time = get_time_in_seconds();
    if (CKKSKeyGen(&keyPair) != 0) {
        fprintf(stderr, "Key generation failed.\n");
        CKKSKeyPairFree(&keyPair);
        return 1;
    }
    end_time = get_time_in_seconds();
    printf("CKKS key generation completed (degree %zu)\n", degree);
    print_performance("Key generation", end_time - start_time, 1);

    const double plaintextValues[] = {1.125, -0.75, 3.5, 0.25, -1.0, 2.75, -3.125, 0.875};
    const size_t plaintextCount = sizeof(plaintextValues) / sizeof(plaintextValues[0]);

    CKKSPolynomial plaintext;
    if (CKKSPolynomialInit(&plaintext, degree, modulus) != 0) {
        fprintf(stderr, "Plaintext polynomial initialization failed.\n");
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    CKKSEncodeFloatsToPolynomial(plaintextValues, plaintextCount, degree, scalingFactor, 
                                 modulus, plaintext.coeffs);

    CKKSCiphertext ciphertext;
    if (CKKSCiphertextInit(&ciphertext, degree, modulus) != 0) {
        fprintf(stderr, "Ciphertext initialization failed.\n");
        CKKSPolynomialFree(&plaintext);
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    /* Encryption performance test */
    const size_t encrypt_iterations = 1000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < encrypt_iterations; i++) {
        if (CKKSEncrypt(&keyPair, &plaintext, &ciphertext) != 0) {
            fprintf(stderr, "Encryption failed.\n");
            CKKSCiphertextFree(&ciphertext);
            CKKSPolynomialFree(&plaintext);
            CKKSKeyPairFree(&keyPair);
            return 1;
        }
    }
    end_time = get_time_in_seconds();
    printf("Encryption operation completed\n");
    print_performance("Encryption operation", end_time - start_time, encrypt_iterations);

    CKKSPolynomial decrypted;
    if (CKKSPolynomialInit(&decrypted, degree, modulus) != 0) {
        fprintf(stderr, "Decrypted polynomial initialization failed.\n");
        CKKSCiphertextFree(&ciphertext);
        CKKSPolynomialFree(&plaintext);
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    /* Decryption performance test */
    const size_t decrypt_iterations = 1000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < decrypt_iterations; i++) {
        if (CKKSDecrypt(&keyPair, &ciphertext, &decrypted) != 0) {
            fprintf(stderr, "Decryption failed.\n");
            CKKSPolynomialFree(&decrypted);
            CKKSCiphertextFree(&ciphertext);
            CKKSPolynomialFree(&plaintext);
            CKKSKeyPairFree(&keyPair);
            return 1;
        }
    }
    end_time = get_time_in_seconds();
    printf("Decryption operation completed\n");
    print_performance("Decryption operation", end_time - start_time, decrypt_iterations);

    // Clean up resources
    CKKSPolynomialFree(&decrypted);
    CKKSCiphertextFree(&ciphertext);
    CKKSPolynomialFree(&plaintext);
    CKKSKeyPairFree(&keyPair);

    return 0;
}

/* ============================================
 * Homomorphic Operation Performance Tests
 * ============================================ */

int test_homomorphic_ops_benchmark(size_t degree) {
    double start_time, end_time;

    const uint64_t modulus = 65537;
    const double scalingFactor = 64.0;

    CKKSKeyPair keyPair;
    if (CKKSKeyPairInit(&keyPair, degree, modulus) != 0) {
        fprintf(stderr, "Key pair initialization failed.\n");
        return 1;
    }

    start_time = get_time_in_seconds();
    if (CKKSKeyGen(&keyPair) != 0) {
        fprintf(stderr, "Key generation failed.\n");
        CKKSKeyPairFree(&keyPair);
        return 1;
    }
    end_time = get_time_in_seconds();
    printf("CKKS key generation completed (degree %zu)\n", degree);
    print_performance("Key generation", end_time - start_time, 1);

    const double plaintextValuesA[] = {1.125, -0.75, 3.5, 0.25, -1.0, 2.75, -3.125, 0.875};
    const double plaintextValuesB[] = {0.5, 1.25, -2.0, 0.75, -1.5, 1.875, 2.25, -0.625};
    const size_t plaintextCount = sizeof(plaintextValuesA) / sizeof(plaintextValuesA[0]);

    CKKSPolynomial plaintextA = {0};
    CKKSPolynomial plaintextB = {0};
    if (CKKSPolynomialInit(&plaintextA, degree, modulus) != 0 ||
        CKKSPolynomialInit(&plaintextB, degree, modulus) != 0) {
        fprintf(stderr, "Plaintext polynomial initialization failed.\n");
        CKKSPolynomialFree(&plaintextA);
        CKKSPolynomialFree(&plaintextB);
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    CKKSEncodeFloatsToPolynomial(plaintextValuesA, plaintextCount, degree, scalingFactor, 
                                 modulus, plaintextA.coeffs);
    CKKSEncodeFloatsToPolynomial(plaintextValuesB, plaintextCount, degree, scalingFactor, 
                                 modulus, plaintextB.coeffs);

    CKKSCiphertext ciphertextA = {0};
    CKKSCiphertext ciphertextB = {0};
    if (CKKSCiphertextInit(&ciphertextA, degree, modulus) != 0 ||
        CKKSCiphertextInit(&ciphertextB, degree, modulus) != 0) {
        fprintf(stderr, "Ciphertext initialization failed.\n");
        CKKSCiphertextFree(&ciphertextA);
        CKKSCiphertextFree(&ciphertextB);
        CKKSPolynomialFree(&plaintextA);
        CKKSPolynomialFree(&plaintextB);
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    if (CKKSEncrypt(&keyPair, &plaintextA, &ciphertextA) != 0 ||
        CKKSEncrypt(&keyPair, &plaintextB, &ciphertextB) != 0) {
        fprintf(stderr, "Encryption failed.\n");
        CKKSCiphertextFree(&ciphertextA);
        CKKSCiphertextFree(&ciphertextB);
        CKKSPolynomialFree(&plaintextA);
        CKKSPolynomialFree(&plaintextB);
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    CKKSCiphertext cipherAdd = {0};
    CKKSCiphertext cipherSub = {0};
    CKKSCiphertext cipherMul = {0};
    if (CKKSCiphertextInit(&cipherAdd, degree, modulus) != 0 ||
        CKKSCiphertextInit(&cipherSub, degree, modulus) != 0 ||
        CKKSCiphertextInit(&cipherMul, degree, modulus) != 0) {
        fprintf(stderr, "Ciphertext operation result initialization failed.\n");
        return 1;
    }

    /* Homomorphic addition performance test */
    const size_t homo_add_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < homo_add_iterations; i++) {
        if (CKKSCiphertextAdd(&ciphertextA, &ciphertextB, &cipherAdd) != 0) {
            fprintf(stderr, "Ciphertext addition failed.\n");
            return 1;
        }
    }
    end_time = get_time_in_seconds();
    print_performance("Ciphertext homomorphic addition", end_time - start_time, homo_add_iterations);

    /* Homomorphic subtraction performance test */
    const size_t homo_sub_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < homo_sub_iterations; i++) {
        if (CKKSCiphertextSub(&ciphertextA, &ciphertextB, &cipherSub) != 0) {
            fprintf(stderr, "Ciphertext subtraction failed.\n");
            return 1;
        }
    }
    end_time = get_time_in_seconds();
    print_performance("Ciphertext homomorphic subtraction", end_time - start_time, homo_sub_iterations);

    /* Homomorphic multiplication performance test */
    const size_t homo_mul_iterations = 1000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < homo_mul_iterations; i++) {
        if (CKKSCiphertextMul(&keyPair, &ciphertextA, &ciphertextB, scalingFactor, &cipherMul) != 0) {
            fprintf(stderr, "Ciphertext multiplication failed.\n");
            return 1;
        }
    }
    end_time = get_time_in_seconds();
    print_performance("Ciphertext homomorphic multiplication", end_time - start_time, homo_mul_iterations);

    printf("Homomorphic operation performance test completed\n");

    // Clean up resources
    CKKSCiphertextFree(&cipherAdd);
    CKKSCiphertextFree(&cipherSub);
    CKKSCiphertextFree(&cipherMul);
    CKKSCiphertextFree(&ciphertextA);
    CKKSCiphertextFree(&ciphertextB);
    CKKSPolynomialFree(&plaintextA);
    CKKSPolynomialFree(&plaintextB);
    CKKSKeyPairFree(&keyPair);

    return 0;
}

/* ============================================
 * Main Program Entry
 * ============================================ */

int main(void) {
    printf("\n");
    printf("╔════════════════════════════════════════════════════════════╗\n");
    printf("║          CKKS Homomorphic Encryption Performance Suite      ║\n");
    printf("╠════════════════════════════════════════════════════════════╣\n");
    printf("║  Test Items: Polynomial ops, Encode/Decode, Encrypt/Decrypt, Homomorphic ops ║\n");
    printf("║  Multi-scale tests: 512, 1024, 2048, 4096 dimensions      ║\n");
    printf("╚════════════════════════════════════════════════════════════╝\n");

    /* Define polynomial degrees for testing */
    const size_t test_degrees[] = {512, 1024, 2048, 4096};
    const size_t num_degrees = sizeof(test_degrees) / sizeof(test_degrees[0]);

    double total_start = get_time_in_seconds();
    int all_passed = 1;

    /* Run complete tests for each degree */
    for (size_t d = 0; d < num_degrees; d++) {
        size_t degree = test_degrees[d];
        
        printf("\n");
        printf("╔════════════════════════════════════════════════════════════╗\n");
        printf("║                Performance Test for Degree N = %-6zu        ║\n", degree);
        printf("╚════════════════════════════════════════════════════════════╝\n");

        printf("\n[Test 1/4] Polynomial operation performance test (N=%zu)...\n", degree);
        int result1 = test_polynomial_benchmark(degree);

        printf("\n[Test 2/4] Encode/decode performance test (N=%zu)...\n", degree);
        int result2 = test_encode_decode_benchmark(degree);

        printf("\n[Test 3/4] Encryption/decryption performance test (N=%zu)...\n", degree);
        int result3 = test_encryption_benchmark(degree);

        printf("\n[Test 4/4] Homomorphic operation performance test (N=%zu)...\n", degree);
        int result4 = test_homomorphic_ops_benchmark(degree);

        printf("\n");
        printf("╔════════════════════════════════════════════════════════════╗\n");
        printf("║              Test Result Summary for N=%-6zu            ║\n", degree);
        printf("╠════════════════════════════════════════════════════════════╣\n");
        printf("║  Polynomial operation test:     %s                              ║\n",
               result1 == 0 ? "Passed" : "Failed");
        printf("║  Encode/decode test:            %s                              ║\n",
               result2 == 0 ? "Passed" : "Failed");
        printf("║  Encryption/decryption test:    %s                              ║\n",
               result3 == 0 ? "Passed" : "Failed");
        printf("║  Homomorphic operation test:    %s                              ║\n",
               result4 == 0 ? "Passed" : "Failed");
        printf("╚════════════════════════════════════════════════════════════╝\n");

        if (result1 != 0 || result2 != 0 || result3 != 0 || result4 != 0) {
            all_passed = 0;
        }
    }

    double total_end = get_time_in_seconds();

    printf("\n");
    printf("╔════════════════════════════════════════════════════════════╗\n");
    printf("║                  Overall Test Summary                      ║\n");
    printf("╠════════════════════════════════════════════════════════════╣\n");
    printf("║  Number of test degrees:  %zu                             ║\n", num_degrees);
    printf("║  Overall result:           %s                              ║\n",
           all_passed ? "All passed" : "Some failed");
    printf("║  Total execution time:     %.6f seconds                   ║\n",
           total_end - total_start);
    printf("╚════════════════════════════════════════════════════════════╝\n");
    printf("\n");

    return all_passed ? 0 : 1;
}

