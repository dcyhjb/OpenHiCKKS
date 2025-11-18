#ifndef CKKS_POLYNOMIAL_H
#define CKKS_POLYNOMIAL_H

#include <stdint.h>
#include <stddef.h>

/**
 * @file ckks_polynomial.h
 * @brief CKKS homomorphic encryption polynomial ring operations library
 * 
 * This library implements polynomial ring operations required for the CKKS homomorphic encryption scheme, including:
 * - Polynomial addition, subtraction, and multiplication
 * - Scalar operations
 * - Polynomial multiplication accelerated using Number Theoretic Transform (NTT)
 * - Modular arithmetic support
 * 
 * All operations are performed on the ring R = Z[X]/(X^N + 1), where N is a power of 2.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief CKKS polynomial structure
 * 
 * Represents a polynomial over a finite field defined by modulus,
 * in the form: coeffs[0] + coeffs[1]*x + ... + coeffs[degree-1]*x^(degree-1)
 */
typedef struct {
    size_t degree;      /**< Degree of the polynomial */
    uint64_t modulus;   /**< Modulus defining the finite field for coefficients */
    uint64_t *coeffs;   /**< Polynomial coefficient array of length degree */
} CKKSPolynomial;

/* ====== Memory Management Functions ====== */

/**
 * @brief Initialize polynomial
 * 
 * Allocate polynomial coefficient array and set initial parameters.
 * 
 * @param poly Pointer to the polynomial to initialize
 * @param degree Degree of the polynomial
 * @param modulus Modulus
 * @return Returns 0 on success, -1 on failure
 */
 int CKKSPolynomialInit(CKKSPolynomial *poly, size_t degree, uint64_t modulus);

/**
 * @brief Free polynomial resources
 * 
 * Free polynomial coefficient array and reset structure to initial state.
 * 
 * @param poly Pointer to the polynomial to free
 */
 void CKKSPolynomialFree(CKKSPolynomial *poly);

/**
 * @brief Set polynomial coefficients
 * 
 * Copy coefficients from external array to polynomial structure, automatically performing modular reduction.
 * 
 * @param poly Target polynomial
 * @param coeffs Coefficient array
 * @param count Number of coefficients (must equal polynomial degree)
 * @return Returns 0 on success, -1 on failure
 */
int CKKSPolynomialSet(CKKSPolynomial *poly, const uint64_t *coeffs, size_t count);

/* ====== Basic Arithmetic Operations ====== */

/**
 * @brief Polynomial addition
 * 
 * Compute result = a + b on the ring R = Z[X]/(X^N + 1).
 * 
 * @param a First polynomial
 * @param b Second polynomial
 * @param result Polynomial to store the result
 * @return Returns 0 on success, -1 on failure
 */
 int CKKSAdd(const CKKSPolynomial *a, const CKKSPolynomial *b, CKKSPolynomial *result);

/**
 * @brief Polynomial subtraction
 * 
 * Compute result = a - b on the ring R = Z[X]/(X^N + 1).
 * 
 * @param a Minuend polynomial
 * @param b Subtrahend polynomial
 * @param result Polynomial to store the result
 * @return Returns 0 on success, -1 on failure
 */

int CKKSPolynomialAddInplace(CKKSPolynomial *accumulator, const CKKSPolynomial *summand);


 int  CKKSSub(const CKKSPolynomial *a, const CKKSPolynomial *b, CKKSPolynomial *result);

/**
 * @brief Polynomial multiplication (accelerated with NTT)
 * 
 * Compute result = a * b mod (X^degree + 1) using Number Theoretic Transform for acceleration.
 * This is the most complex operation in the CKKS scheme, using NTT to reduce time complexity from O(N^2) to O(N log N).
 * 
 * @param a First polynomial
 * @param b Second polynomial
 * @param result Polynomial to store the result
 * @return Returns 0 on success, -1 on failure
 */
 int CKKSMul(const CKKSPolynomial *a, const CKKSPolynomial *b, CKKSPolynomial *result);

/**
 * @brief Scalar division
 * 
 * Compute result = poly / scalar by multiplying each coefficient by the modular inverse of the scalar.
 * 
 * @param poly Dividend polynomial
 * @param scalar Scalar divisor (must not be 0)
 * @param result Polynomial to store the result
 * @return Returns 0 on success, -1 on failure (scalar has no inverse or parameter error)
 */
int CKKSScalarDivide(const CKKSPolynomial *poly, uint64_t scalar, CKKSPolynomial *result);

/* ====== Helper Functions ====== */

/**
 * @brief Check if polynomial coefficients match expected values
 * 
 * Used for testing and verifying correctness of operation results.
 * 
 * @param poly Polynomial to check
 * @param expected Expected coefficient array
 * @return Returns 1 if matched, 0 if not matched
 */
int CKKSMatchCoefficients(const CKKSPolynomial *poly, const uint64_t *expected);

/**
 * @brief Naive polynomial multiplication (for verification)
 * 
 * Compute polynomial multiplication using direct convolution method, used to verify correctness of NTT implementation.
 * Note: This function is for testing only and has low performance.
 * 
 * @param a First coefficient array
 * @param b Second coefficient array
 * @param degree Polynomial degree
 * @param modulus Modulus
 * @param result Array to store result (length at least degree)
 */
void CKKSNaiveMulReduce(const uint64_t *a, const uint64_t *b, size_t degree, 
                        uint64_t modulus, uint64_t *result);

/* ====== Low-level Math Functions ====== */

/**
 * @brief Modular reduction function (64-bit)
 * 
 * Reduce signed 64-bit integer to the range [0, modulus-1].
 * 
 * @param value Value to reduce
 * @param modulus Modulus
 * @return Reduced value
 */
 uint64_t CKKSModReduceInt64(int64_t value, uint64_t modulus);

/**
 * @brief Modular reduction function (128-bit)
 * 
 * Handle modular arithmetic for large number multiplication results, using 128-bit integer intermediate results.
 * 
 * @param value Value to reduce (128-bit)
 * @param modulus Modulus
 * @return Reduced value
 */
uint64_t CKKSModReduceInt128(__int128 value, uint64_t modulus);

/**
 * @brief Modular exponentiation
 * 
 * Compute base^exponent mod modulus using fast exponentiation algorithm.
 * 
 * @param base Base
 * @param exponent Exponent
 * @param modulus Modulus
 * @return Computation result
 */
uint64_t CKKSPolynomialPowMod(uint64_t base, uint64_t exponent, uint64_t modulus);

/**
 * @brief Modular multiplication
 * 
 * Compute a * b mod modulus, using 128-bit intermediate result to avoid overflow.
 * 
 * @param a First operand
 * @param b Second operand
 * @param modulus Modulus
 * @return Computation result
 */
 uint64_t CKKSPolynomialMulMod(uint64_t a, uint64_t b, uint64_t modulus);

/**
 * @brief Modular inverse computation
 * 
 * Compute the modular inverse of value modulo modulus using extended Euclidean algorithm.
 * 
 * @param value Value for which to compute the inverse
 * @param modulus Modulus
 * @return Modular inverse, returns 0 if it does not exist
 */
uint64_t CKKSPolynomialModInverse(uint64_t value, uint64_t modulus);

/* ====== Number Theoretic Transform Functions ====== */

/**
 * @brief Number Theoretic Transform (NTT)
 * 
 * Transform polynomial from time domain to frequency domain, used to accelerate polynomial multiplication.
 * 
 * @param data Input/output data array
 * @param length Data length (must be a power of 2)
 * @param modulus Modulus
 * @param primitiveRoot Primitive root
 */
 void CKKSNTT(uint64_t *data, size_t length, uint64_t modulus, uint64_t primitiveRoot);

/**
 * @brief Inverse Number Theoretic Transform (Inverse NTT)
 * 
 * Transform polynomial from frequency domain back to time domain.
 * 
 * @param data Input/output data array
 * @param length Data length (must be a power of 2)
 * @param modulus Modulus
 * @param primitiveRoot Primitive root
 */
void CKKSInverseNTT(uint64_t *data, size_t length, uint64_t modulus, uint64_t primitiveRoot);

void CKKSPolynomialPrint(const CKKSPolynomial *poly, const char *label, size_t limit);


#ifdef __cplusplus
}
#endif

#endif /* CKKS_POLYNOMIAL_H */