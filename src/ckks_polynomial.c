#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ckks_polynomial.h"

// /* Modular inverse function declaration */
// // static uint64_t CKKSPolynomialModInverse(uint64_t value, uint64_t modulus);

// /* 
//  * 64-bit integer modular reduction function
//  * Reduce 64-bit signed integer to range [0, modulus-1]
//  */
uint64_t CKKSModReduceInt64(int64_t value, uint64_t modulus)
{
    int64_t reduced = value % (int64_t)modulus;
    if (reduced < 0) {
        reduced += (int64_t)modulus;
    }
    return (uint64_t)reduced;
}

/* 
 * 128-bit integer modular reduction function
 * Handle modular arithmetic for large number multiplication results
 */
uint64_t CKKSModReduceInt128(__int128 value, uint64_t modulus)
{
    __int128 reduced = value % (__int128)modulus;
    if (reduced < 0) {
        reduced += (__int128)modulus;
    }
    return (uint64_t)reduced;
}

/* 
 * Initialize polynomial
 * Allocate coefficient array memory and set degree and modulus
 */
 int CKKSPolynomialInit(CKKSPolynomial *poly, size_t degree, uint64_t modulus)
{
    if (poly == NULL || degree == 0 || modulus == 0) {
        return -1;
    }
    poly->degree = degree;
    poly->modulus = modulus;
    poly->coeffs = (uint64_t *)calloc(degree, sizeof(uint64_t));
    return poly->coeffs == NULL ? -1 : 0;
}

/* Free polynomial resources */
void CKKSPolynomialFree(CKKSPolynomial *poly)
{
    if (poly == NULL) {
        return;
    }
    free(poly->coeffs);
    poly->coeffs = NULL;
    poly->degree = 0;
    poly->modulus = 0;
}

/* 
 * Set polynomial coefficients
 * Copy external coefficient array to polynomial structure and perform modular reduction
 */
 int CKKSPolynomialSet(CKKSPolynomial *poly, const uint64_t *coeffs, size_t count)
{
    if (poly == NULL || coeffs == NULL || count != poly->degree) {
        return -1;
    }
    (void)memcpy(poly->coeffs, coeffs, poly->degree * sizeof(uint64_t));
    for (size_t i = 0; i < poly->degree; ++i) {
        poly->coeffs[i] %= poly->modulus;
    }
    return 0;
}

/* 
 * Polynomial addition
 * Add corresponding coefficients and perform modular reduction
 */
 int CKKSAdd(const CKKSPolynomial *a, const CKKSPolynomial *b, CKKSPolynomial *result)
{
    if (a == NULL || b == NULL || result == NULL) {
        return -1;
    }
    if (a->degree != b->degree || a->degree != result->degree || a->modulus != b->modulus ||
        a->modulus != result->modulus) {
        return -1;
    }
    for (size_t i = 0; i < a->degree; ++i) {
        int64_t value = (int64_t)a->coeffs[i] + (int64_t)b->coeffs[i];
        result->coeffs[i] = CKKSModReduceInt64(value, a->modulus);
    }
    return 0;
}

/* 
 * Polynomial subtraction
 * Subtract corresponding coefficients and perform modular reduction
 */
 int  CKKSSub(const CKKSPolynomial *a, const CKKSPolynomial *b, CKKSPolynomial *result)
{
    if (a == NULL || b == NULL || result == NULL) {
        return -1;
    }
    if (a->degree != b->degree || a->degree != result->degree || a->modulus != b->modulus ||
        a->modulus != result->modulus) {
        return -1;
    }
    for (size_t i = 0; i < a->degree; ++i) {
        int64_t value = (int64_t)a->coeffs[i] - (int64_t)b->coeffs[i];
        result->coeffs[i] = CKKSModReduceInt64(value, a->modulus);
    }
    return 0;
}

int CKKSPolynomialAddInplace(CKKSPolynomial *accumulator, const CKKSPolynomial *summand)
{
    if (accumulator == NULL || summand == NULL || accumulator->degree != summand->degree ||
        accumulator->modulus != summand->modulus) {
        return -1;
    }
    for (size_t i = 0; i < accumulator->degree; ++i) {
        uint64_t sum = accumulator->coeffs[i] + summand->coeffs[i];
        if (sum >= accumulator->modulus) {
            sum -= accumulator->modulus;
        }
        accumulator->coeffs[i] = sum;
    }
    return 0;
}

/* 
 * Modular exponentiation (fast power algorithm)
 * Compute base^exponent mod modulus
 */
 uint64_t CKKSPolynomialPowMod(uint64_t base, uint64_t exponent, uint64_t modulus)
{
    uint64_t result = 1;
    base %= modulus;
    while (exponent > 0) {
        if (exponent & 1U) {
            result = CKKSModReduceInt128((__int128)result * base, modulus);
        }
        base = CKKSModReduceInt128((__int128)base * base, modulus);
        exponent >>= 1U;
    }
    return result;
}

/* Modular multiplication: compute a * b mod modulus */
 uint64_t CKKSPolynomialMulMod(uint64_t a, uint64_t b, uint64_t modulus)
{
    return CKKSModReduceInt128((__int128)a * (__int128)b, modulus);
}

/* 
 * Bit reversal function
 * Used for bit-reverse permutation in NTT
 */
size_t CKKSReverseBits(size_t value, unsigned int bitCount)
{
    size_t reversed = 0;
    for (unsigned int i = 0; i < bitCount; ++i) {
        reversed <<= 1U;
        reversed |= (value & 1U);
        value >>= 1U;
    }
    return reversed;
}

/* 
 * Bit-reverse permutation
 * Rearrange array elements in bit-reversed order
 */
 void CKKSBitReverse(uint64_t *data, size_t length)
{
    unsigned int bitCount = 0;
    size_t temp = length;
    while (temp > 1) {
        ++bitCount;
        temp >>= 1U;
    }
    for (size_t i = 0; i < length; ++i) {
        size_t j = CKKSReverseBits(i, bitCount);
        if (j > i) {
            uint64_t swap = data[i];
            data[i] = data[j];
            data[j] = swap;
        }
    }
}

/* 
 * Number Theoretic Transform (NTT)
 * Transform polynomial from time domain to frequency domain to accelerate polynomial multiplication
 */
void CKKSNTT(uint64_t *data, size_t length, uint64_t modulus, uint64_t primitiveRoot)
{
    CKKSBitReverse(data, length);
    
    // Butterfly operation
    for (size_t len = 2; len <= length; len <<= 1U) {
        size_t half = len >> 1U;
        // Compute twiddle factor
        uint64_t rootStep = CKKSPolynomialPowMod(primitiveRoot, length / len, modulus);
        
        for (size_t i = 0; i < length; i += len) {
            uint64_t w = 1;  // Initial twiddle factor
            for (size_t j = 0; j < half; ++j) {
                uint64_t u = data[i + j];
                uint64_t v = CKKSPolynomialMulMod(data[i + j + half], w, modulus);
                
                // Core butterfly operation computation
                uint64_t add = u + v;
                if (add >= modulus) {
                    add -= modulus;  // Modular reduction
                }
                uint64_t sub = u >= v ? u - v : u + modulus - v;
                
                data[i + j] = add;
                data[i + j + half] = sub;
                
                w = CKKSPolynomialMulMod(w, rootStep, modulus);  // Update twiddle factor
            }
        }
    }
}

/* 
 * Inverse Number Theoretic Transform (Inverse NTT)
 * Transform polynomial from frequency domain back to time domain
 */
 void CKKSInverseNTT(uint64_t *data, size_t length, uint64_t modulus, uint64_t primitiveRoot)
{
    CKKSBitReverse(data, length);  // Bit-reverse permutation
    
    // Butterfly operation (similar to NTT but in opposite direction)
    for (size_t len = 2; len <= length; len <<= 1U) {
        size_t half = len >> 1U;
        // Compute twiddle factor
        uint64_t rootStep = CKKSPolynomialPowMod(primitiveRoot, length / len, modulus);
        
        for (size_t i = 0; i < length; i += len) {
            uint64_t w = 1;  // Initial twiddle factor
            for (size_t j = 0; j < half; ++j) {
                uint64_t u = data[i + j];
                uint64_t v = CKKSPolynomialMulMod(data[i + j + half], w, modulus);
                
                // Core butterfly operation computation
                data[i + j] = u + v;
                if (data[i + j] >= modulus) {
                    data[i + j] -= modulus;  // Modular reduction
                }
                uint64_t sub = u >= v ? u - v : u + modulus - v;
                data[i + j + half] = sub;
                
                w = CKKSPolynomialMulMod(w, rootStep, modulus);  // Update twiddle factor
            }
        }
    }
    
    // Multiply by inverse of length to complete inverse transform
    uint64_t invLength = CKKSPolynomialModInverse((uint64_t)length, modulus);
    for (size_t i = 0; i < length; ++i) {
        data[i] = CKKSPolynomialMulMod(data[i], invLength, modulus);
    }
}

/* 
 * Polynomial multiplication (using NTT acceleration)
 * Compute a * b mod (x^degree + 1)
 */
 int CKKSMul(const CKKSPolynomial *a, const CKKSPolynomial *b, CKKSPolynomial *result)
{
    if (a == NULL || b == NULL || result == NULL) {
        return -1;  // Null pointer check
    }
    // Check if degrees and moduli match
    if (a->degree != b->degree || a->degree != result->degree || a->modulus != b->modulus ||
        a->modulus != result->modulus) {
        return -1;
    }

    const size_t degree = a->degree;
    // Calculate transform length needed for NTT (power of 2)
    size_t transformSize = 1;
    while (transformSize < degree * 2U) {
        transformSize <<= 1U;
    }

    // Allocate buffers for NTT transform
    uint64_t *fa = (uint64_t *)calloc(transformSize, sizeof(uint64_t));
    uint64_t *fb = (uint64_t *)calloc(transformSize, sizeof(uint64_t));
    if (fa == NULL || fb == NULL) {
        free(fa);
        free(fb);
        return -1;  // Memory allocation failed
    }

    // Copy polynomial coefficients to buffers
    for (size_t i = 0; i < degree; ++i) {
        fa[i] = a->coeffs[i];
        fb[i] = b->coeffs[i];
    }

    // Set NTT parameters
    const uint64_t primitiveRoot = 3; /* Primitive root modulo 65537 */
    uint64_t root = CKKSPolynomialPowMod(primitiveRoot, (a->modulus - 1) / transformSize, a->modulus);
    
    // Perform NTT transform
    CKKSNTT(fa, transformSize, a->modulus, root);
    CKKSNTT(fb, transformSize, a->modulus, root);

    // Point-wise multiplication in frequency domain
    for (size_t i = 0; i < transformSize; ++i) {
        fa[i] = CKKSPolynomialMulMod(fa[i], fb[i], a->modulus);
    }

    // Inverse NTT transform
    uint64_t invRoot = CKKSPolynomialModInverse(root, a->modulus);
    CKKSInverseNTT(fa, transformSize, a->modulus, invRoot);

    // Copy results
    for (size_t i = 0; i < degree; ++i) {
        result->coeffs[i] = fa[i];
    }

    // Handle reduction modulo x^degree + 1
    for (size_t i = degree; i < transformSize && i < degree * 2U; ++i) {
        size_t target = i - degree;
        uint64_t value = fa[i];
        if (value != 0) {
            // Perform reduction x^degree = -1
            if (result->coeffs[target] >= value) {
                result->coeffs[target] -= value;
            } else {
                result->coeffs[target] = result->coeffs[target] + a->modulus - value;
            }
        }
    }

    // Free resources
    free(fa);
    free(fb);
    return 0;
}

/* 
 * Modular inverse computation (Extended Euclidean Algorithm)
 * Compute the inverse of value modulo modulus
 */
uint64_t CKKSPolynomialModInverse(uint64_t value, uint64_t modulus)
{
    int64_t t = 0;
    int64_t newT = 1;
    int64_t r = (int64_t)modulus;
    int64_t newR = (int64_t)(value % modulus);

    // Extended Euclidean Algorithm
    while (newR != 0) {
        int64_t quotient = r / newR;
        int64_t tempT = newT;
        newT = t - quotient * newT;
        t = tempT;

        int64_t tempR = newR;
        newR = r - quotient * newR;
        r = tempR;
    }

    if (r != 1) {
        return 0; /* Inverse does not exist */
    }

    if (t < 0) {
        t += (int64_t)modulus;  // Adjust to positive range
    }
    return (uint64_t)t;
}

/* 
 * Scalar division
 * Divide each coefficient of polynomial by scalar
 */
 int CKKSScalarDivide(const CKKSPolynomial *poly, uint64_t scalar, CKKSPolynomial *result)
{
    if (poly == NULL || result == NULL || scalar == 0) {
        return -1;  // Parameter check
    }
    if (poly->degree != result->degree || poly->modulus != result->modulus) {
        return -1;  // Degree and modulus check
    }

    // Compute modular inverse of scalar
    uint64_t inverse = CKKSPolynomialModInverse(scalar, poly->modulus);
    if (inverse == 0) {
        return -1;  // Inverse does not exist
    }

    // Multiply each coefficient by inverse to implement division
    for (size_t i = 0; i < poly->degree; ++i) {
        __int128 product = (__int128)poly->coeffs[i] * (__int128)inverse;
        result->coeffs[i] = CKKSModReduceInt128(product, poly->modulus);
    }
    return 0;
}

/* 
 * Coefficient matching check
 * Verify if polynomial coefficients match expected values
 */
 int CKKSMatchCoefficients(const CKKSPolynomial *poly, const uint64_t *expected)
{
    if (poly == NULL || expected == NULL) {
        return 0;  // Null pointer check
    }
    for (size_t i = 0; i < poly->degree; ++i) {
        if (poly->coeffs[i] != expected[i]) {
            return 0;  // Coefficients do not match
        }
    }
    return 1;  // All coefficients match
}

/* 
 * Naive polynomial multiplication (for verifying correctness of NTT results)
 * Directly compute polynomial multiplication then perform modular reduction
 */
 void CKKSNaiveMulReduce(const uint64_t *a, const uint64_t *b, size_t degree, uint64_t modulus,
                               uint64_t *result)
{
    if (a == NULL || b == NULL || result == NULL) {
        return;  // Null pointer check
    }

    size_t convLength = degree * 2U;
    uint64_t *temp = (uint64_t *)calloc(convLength, sizeof(uint64_t));
    assert(temp != NULL);

    // Direct convolution computation
    for (size_t i = 0; i < degree; ++i) {
        for (size_t j = 0; j < degree; ++j) {
            size_t idx = i + j;
            uint64_t prod = CKKSPolynomialMulMod(a[i], b[j], modulus);
            uint64_t sum = temp[idx] + prod;
            if (sum >= modulus) {
                sum -= modulus;  // Modular reduction
            }
            temp[idx] = sum;
        }
    }

    // Copy lower degree terms
    for (size_t i = 0; i < degree; ++i) {
        result[i] = temp[i];
    }

    // Reduction modulo x^degree + 1
    for (size_t idx = degree; idx < convLength; ++idx) {
        uint64_t value = temp[idx];
        if (value == 0) {
            continue;
        }
        size_t target = idx - degree;
        // Perform reduction x^degree = -1
        if (result[target] >= value) {
            result[target] -= value;
        } else {
            result[target] = result[target] + modulus - value;
        }
    }

    free(temp);
}



void CKKSPolynomialCopy(const CKKSPolynomial *src, CKKSPolynomial *dst)
{
    if (src == NULL || dst == NULL || src->degree != dst->degree || src->modulus != dst->modulus) {
        return;
    }
    memcpy(dst->coeffs, src->coeffs, src->degree * sizeof(uint64_t));
}


void CKKSPolynomialSetZero(CKKSPolynomial *poly)
{
    if (poly == NULL) {
        return;
    }
    memset(poly->coeffs, 0, poly->degree * sizeof(uint64_t));
}


void CKKSPolynomialPrint(const CKKSPolynomial *poly, const char *label, size_t limit)
{
    if (poly == NULL) {
        return;
    }
    printf("%s coefficients (mod %llu):\n", label != NULL ? label : "polynomial", (unsigned long long)poly->modulus);
    size_t show = limit < poly->degree ? limit : poly->degree;
    for (size_t i = 0; i < show; ++i) {
        printf("  coeff[%zu] = %lld\n", i, (long long)poly->coeffs[i]);
    }
}


