/*
 * CKKS 同态加密系统性能测试套件
 * 
 * 本文件包含所有性能基准测试代码，用于评估 CKKS 实现的各项操作性能
 * 
 * 测试项目包括：
 * 1. 多项式运算性能（加法、减法、乘法、除法）
 * 2. 编码/解码性能
 * 3. 加密/解密性能
 * 4. 同态运算性能（密文加法、减法、乘法）
 * 
 * 支持多种数据规模测试，便于评估不同维度下的性能表现
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
 * 获取当前时间（秒）
 * 使用 CLOCK_MONOTONIC 保证单调递增，适合性能测试
 */
static double get_time_in_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}

/**
 * 打印性能统计信息
 * @param operation 操作名称
 * @param time_seconds 执行时间（秒）
 * @param iterations 迭代次数
 */
static void print_performance(const char *operation, double time_seconds, size_t iterations) {
    printf("⏱️  %s: %.6f 秒", operation, time_seconds);
    if (iterations > 1) {
        printf(" (%.2f 次/秒, 平均 %.6f 秒/次)\n", 
               iterations / time_seconds, time_seconds / iterations);
    } else {
        printf("\n");
    }
}

/* ============================================
 * 多项式运算性能测试
 * ============================================ */

int test_polynomial_benchmark(size_t degree) {
    double start_time, end_time;

    const uint64_t modulus = 65537; /* 16位费马素数 */

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

    // 初始化测试数据
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

    /* 多项式加法性能测试 */
    const size_t add_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < add_iterations; i++) {
        assert(CKKSAdd(&a, &b, &result) == 0);
    }
    end_time = get_time_in_seconds();
    assert(CKKSMatchCoefficients(&result, expectedAdd) == 1);
    printf("✅ CKKS 多项式加法测试通过 (维度 %zu)\n", degree);
    print_performance("多项式加法", end_time - start_time, add_iterations);

    /* 多项式减法性能测试 */
    const size_t sub_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < sub_iterations; i++) {
        assert(CKKSSub(&a, &b, &result) == 0);
    }
    end_time = get_time_in_seconds();
    assert(CKKSMatchCoefficients(&result, expectedSub) == 1);
    printf("✅ CKKS 多项式减法测试通过 (维度 %zu)\n", degree);
    print_performance("多项式减法", end_time - start_time, sub_iterations);

    /* 多项式乘法性能测试 */
    const size_t mul_iterations = 1000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < mul_iterations; i++) {
        assert(CKKSMul(&a, &b, &result) == 0);
    }
    end_time = get_time_in_seconds();
    assert(CKKSMatchCoefficients(&result, expectedMul) == 1);
    printf("✅ CKKS 多项式乘法测试通过 (维度 %zu)\n", degree);
    print_performance("多项式乘法", end_time - start_time, mul_iterations);

    /* 多项式标量除法性能测试 */
    const size_t div_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < div_iterations; i++) {
        assert(CKKSScalarDivide(&a, 3, &result) == 0);
    }
    end_time = get_time_in_seconds();
    assert(CKKSMatchCoefficients(&result, expectedDiv) == 1);
    printf("✅ CKKS 多项式标量除法测试通过 (维度 %zu)\n", degree);
    print_performance("多项式标量除法", end_time - start_time, div_iterations);

    // 清理资源
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
 * 编码/解码性能测试
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

    /* 编码性能测试：复数槽到多项式 */
    const size_t encode_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < encode_iterations; i++) {
        CKKSEncodeComplexSlotsToPolynomial(slots, plaintextCount, degree, scalingFactor, modulus, encodedCoeffs);
    }
    end_time = get_time_in_seconds();
    print_performance("编码 (复数槽→多项式)", end_time - start_time, encode_iterations);

    CKKSPolynomial encodedPoly;
    assert(CKKSPolynomialInit(&encodedPoly, degree, modulus) == 0);
    assert(CKKSPolynomialSet(&encodedPoly, encodedCoeffs, degree) == 0);

    CKKSComplex *decodedSlots = (CKKSComplex *)calloc(slotCapacity, sizeof(CKKSComplex));
    double *recoveredFloats = (double *)calloc(plaintextCount, sizeof(double));
    assert(decodedSlots != NULL && recoveredFloats != NULL);

    /* 解码性能测试：多项式到复数槽 */
    const size_t decode_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < decode_iterations; i++) {
        CKKSDecodePolynomialToComplex(encodedPoly.coeffs, degree, scalingFactor, modulus, 
                                      decodedSlots, plaintextCount);
    }
    end_time = get_time_in_seconds();
    print_performance("解码 (多项式→复数槽)", end_time - start_time, decode_iterations);

    /* 直接解码性能测试：多项式到浮点数 */
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < decode_iterations; i++) {
        CKKSDecodePolynomialToFloats(encodedPoly.coeffs, degree, scalingFactor, modulus, 
                                     recoveredFloats, plaintextCount);
    }
    end_time = get_time_in_seconds();
    print_performance("解码 (多项式→浮点数)", end_time - start_time, decode_iterations);

    printf("✅ 编码/解码性能测试完成\n");

    // 清理资源
    free(decodedSlots);
    free(recoveredFloats);
    CKKSPolynomialFree(&encodedPoly);
    free(slots);
    free(encodedCoeffs);

    return 0;
}

/* ============================================
 * 加密/解密性能测试
 * ============================================ */

int test_encryption_benchmark(size_t degree) {
    double start_time, end_time;

    const uint64_t modulus = 65537;
    const double scalingFactor = 128.0;

    CKKSKeyPair keyPair;
    if (CKKSKeyPairInit(&keyPair, degree, modulus) != 0) {
        fprintf(stderr, "密钥对初始化失败。\n");
        return 1;
    }

    /* 密钥生成性能测试 */
    start_time = get_time_in_seconds();
    if (CKKSKeyGen(&keyPair) != 0) {
        fprintf(stderr, "密钥生成失败。\n");
        CKKSKeyPairFree(&keyPair);
        return 1;
    }
    end_time = get_time_in_seconds();
    printf("✅ CKKS 密钥生成完成 (维度 %zu)\n", degree);
    print_performance("密钥生成", end_time - start_time, 1);

    const double plaintextValues[] = {1.125, -0.75, 3.5, 0.25, -1.0, 2.75, -3.125, 0.875};
    const size_t plaintextCount = sizeof(plaintextValues) / sizeof(plaintextValues[0]);

    CKKSPolynomial plaintext;
    if (CKKSPolynomialInit(&plaintext, degree, modulus) != 0) {
        fprintf(stderr, "明文多项式初始化失败。\n");
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    CKKSEncodeFloatsToPolynomial(plaintextValues, plaintextCount, degree, scalingFactor, 
                                 modulus, plaintext.coeffs);

    CKKSCiphertext ciphertext;
    if (CKKSCiphertextInit(&ciphertext, degree, modulus) != 0) {
        fprintf(stderr, "密文初始化失败。\n");
        CKKSPolynomialFree(&plaintext);
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    /* 加密性能测试 */
    const size_t encrypt_iterations = 1000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < encrypt_iterations; i++) {
        if (CKKSEncrypt(&keyPair, &plaintext, &ciphertext) != 0) {
            fprintf(stderr, "加密失败。\n");
            CKKSCiphertextFree(&ciphertext);
            CKKSPolynomialFree(&plaintext);
            CKKSKeyPairFree(&keyPair);
            return 1;
        }
    }
    end_time = get_time_in_seconds();
    printf("✅ 加密操作完成\n");
    print_performance("加密操作", end_time - start_time, encrypt_iterations);

    CKKSPolynomial decrypted;
    if (CKKSPolynomialInit(&decrypted, degree, modulus) != 0) {
        fprintf(stderr, "解密多项式初始化失败。\n");
        CKKSCiphertextFree(&ciphertext);
        CKKSPolynomialFree(&plaintext);
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    /* 解密性能测试 */
    const size_t decrypt_iterations = 1000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < decrypt_iterations; i++) {
        if (CKKSDecrypt(&keyPair, &ciphertext, &decrypted) != 0) {
            fprintf(stderr, "解密失败。\n");
            CKKSPolynomialFree(&decrypted);
            CKKSCiphertextFree(&ciphertext);
            CKKSPolynomialFree(&plaintext);
            CKKSKeyPairFree(&keyPair);
            return 1;
        }
    }
    end_time = get_time_in_seconds();
    printf("✅ 解密操作完成\n");
    print_performance("解密操作", end_time - start_time, decrypt_iterations);

    // 清理资源
    CKKSPolynomialFree(&decrypted);
    CKKSCiphertextFree(&ciphertext);
    CKKSPolynomialFree(&plaintext);
    CKKSKeyPairFree(&keyPair);

    return 0;
}

/* ============================================
 * 同态运算性能测试
 * ============================================ */

int test_homomorphic_ops_benchmark(size_t degree) {
    double start_time, end_time;

    const uint64_t modulus = 65537;
    const double scalingFactor = 64.0;

    CKKSKeyPair keyPair;
    if (CKKSKeyPairInit(&keyPair, degree, modulus) != 0) {
        fprintf(stderr, "密钥对初始化失败。\n");
        return 1;
    }

    start_time = get_time_in_seconds();
    if (CKKSKeyGen(&keyPair) != 0) {
        fprintf(stderr, "密钥生成失败。\n");
        CKKSKeyPairFree(&keyPair);
        return 1;
    }
    end_time = get_time_in_seconds();
    printf("✅ CKKS 密钥生成完成 (维度 %zu)\n", degree);
    print_performance("密钥生成", end_time - start_time, 1);

    const double plaintextValuesA[] = {1.125, -0.75, 3.5, 0.25, -1.0, 2.75, -3.125, 0.875};
    const double plaintextValuesB[] = {0.5, 1.25, -2.0, 0.75, -1.5, 1.875, 2.25, -0.625};
    const size_t plaintextCount = sizeof(plaintextValuesA) / sizeof(plaintextValuesA[0]);

    CKKSPolynomial plaintextA = {0};
    CKKSPolynomial plaintextB = {0};
    if (CKKSPolynomialInit(&plaintextA, degree, modulus) != 0 ||
        CKKSPolynomialInit(&plaintextB, degree, modulus) != 0) {
        fprintf(stderr, "明文多项式初始化失败。\n");
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
        fprintf(stderr, "密文初始化失败。\n");
        CKKSCiphertextFree(&ciphertextA);
        CKKSCiphertextFree(&ciphertextB);
        CKKSPolynomialFree(&plaintextA);
        CKKSPolynomialFree(&plaintextB);
        CKKSKeyPairFree(&keyPair);
        return 1;
    }

    if (CKKSEncrypt(&keyPair, &plaintextA, &ciphertextA) != 0 ||
        CKKSEncrypt(&keyPair, &plaintextB, &ciphertextB) != 0) {
        fprintf(stderr, "加密失败。\n");
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
        fprintf(stderr, "密文运算结果初始化失败。\n");
        return 1;
    }

    /* 同态加法性能测试 */
    const size_t homo_add_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < homo_add_iterations; i++) {
        if (CKKSCiphertextAdd(&ciphertextA, &ciphertextB, &cipherAdd) != 0) {
            fprintf(stderr, "密文加法失败。\n");
            return 1;
        }
    }
    end_time = get_time_in_seconds();
    print_performance("密文同态加法", end_time - start_time, homo_add_iterations);

    /* 同态减法性能测试 */
    const size_t homo_sub_iterations = 10000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < homo_sub_iterations; i++) {
        if (CKKSCiphertextSub(&ciphertextA, &ciphertextB, &cipherSub) != 0) {
            fprintf(stderr, "密文减法失败。\n");
            return 1;
        }
    }
    end_time = get_time_in_seconds();
    print_performance("密文同态减法", end_time - start_time, homo_sub_iterations);

    /* 同态乘法性能测试 */
    const size_t homo_mul_iterations = 1000;
    start_time = get_time_in_seconds();
    for (size_t i = 0; i < homo_mul_iterations; i++) {
        if (CKKSCiphertextMul(&keyPair, &ciphertextA, &ciphertextB, scalingFactor, &cipherMul) != 0) {
            fprintf(stderr, "密文乘法失败。\n");
            return 1;
        }
    }
    end_time = get_time_in_seconds();
    print_performance("密文同态乘法", end_time - start_time, homo_mul_iterations);

    printf("✅ 同态运算性能测试完成\n");

    // 清理资源
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
 * 主程序入口
 * ============================================ */

int main(void) {
    printf("\n");
    printf("╔════════════════════════════════════════════════════════════╗\n");
    printf("║          CKKS 同态加密系统性能测试套件                    ║\n");
    printf("╠════════════════════════════════════════════════════════════╣\n");
    printf("║  测试项目：多项式运算、编码/解码、加密/解密、同态运算    ║\n");
    printf("║  多规模测试：512, 1024, 2048, 4096 维度                  ║\n");
    printf("╚════════════════════════════════════════════════════════════╝\n");

    /* 定义测试的多项式维度 */
    const size_t test_degrees[] = {512, 1024, 2048, 4096};
    const size_t num_degrees = sizeof(test_degrees) / sizeof(test_degrees[0]);

    double total_start = get_time_in_seconds();
    int all_passed = 1;

    /* 对每个维度运行完整测试 */
    for (size_t d = 0; d < num_degrees; d++) {
        size_t degree = test_degrees[d];
        
        printf("\n");
        printf("╔════════════════════════════════════════════════════════════╗\n");
        printf("║                维度 N = %-6zu 的性能测试                ║\n", degree);
        printf("╚════════════════════════════════════════════════════════════╝\n");

        printf("\n[测试 1/4] 多项式运算性能测试 (N=%zu)...\n", degree);
        int result1 = test_polynomial_benchmark(degree);

        printf("\n[测试 2/4] 编码/解码性能测试 (N=%zu)...\n", degree);
        int result2 = test_encode_decode_benchmark(degree);

        printf("\n[测试 3/4] 加密/解密性能测试 (N=%zu)...\n", degree);
        int result3 = test_encryption_benchmark(degree);

        printf("\n[测试 4/4] 同态运算性能测试 (N=%zu)...\n", degree);
        int result4 = test_homomorphic_ops_benchmark(degree);

        printf("\n");
        printf("╔════════════════════════════════════════════════════════════╗\n");
        printf("║              N=%-6zu 测试结果总结                       ║\n", degree);
        printf("╠════════════════════════════════════════════════════════════╣\n");
        printf("║  多项式运算测试:     %s                              ║\n",
               result1 == 0 ? "✅ 通过" : "❌ 失败");
        printf("║  编码/解码测试:      %s                              ║\n",
               result2 == 0 ? "✅ 通过" : "❌ 失败");
        printf("║  加密/解密测试:      %s                              ║\n",
               result3 == 0 ? "✅ 通过" : "❌ 失败");
        printf("║  同态运算测试:       %s                              ║\n",
               result4 == 0 ? "✅ 通过" : "❌ 失败");
        printf("╚════════════════════════════════════════════════════════════╝\n");

        if (result1 != 0 || result2 != 0 || result3 != 0 || result4 != 0) {
            all_passed = 0;
        }
    }

    double total_end = get_time_in_seconds();

    printf("\n");
    printf("╔════════════════════════════════════════════════════════════╗\n");
    printf("║                  全部测试总结                              ║\n");
    printf("╠════════════════════════════════════════════════════════════╣\n");
    printf("║  测试维度数量:       %zu 个                                ║\n", num_degrees);
    printf("║  总体结果:           %s                              ║\n",
           all_passed ? "✅ 全部通过" : "❌ 部分失败");
    printf("║  总执行时间:         %.6f 秒                       ║\n",
           total_end - total_start);
    printf("╚════════════════════════════════════════════════════════════╝\n");
    printf("\n");

    return all_passed ? 0 : 1;
}

