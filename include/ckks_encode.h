#ifndef CKKS_ENCODE_H
#define CKKS_ENCODE_H
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ckks_polynomial.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



#ifdef __cplusplus
extern "C" {
#endif



// Complex number structure for CKKS encoding
typedef struct {
    double real;
    double imag;
} CKKSComplex;


// Function declarations

/**
 * @brief Copies float array to complex slots for encoding
 * 
 * @param input Input float array
 * @param inputCount Number of elements in input array
 * @param slots Output complex slots
 * @param slotCount Number of complex slots
 */
 void CKKSCopyFloatsToComplexSlots(const double *input, size_t inputCount, 
                                  CKKSComplex *slots, size_t slotCount);

/**
 * @brief Encodes complex slots to polynomial coefficients using inverse FFT
 * 
 * @param slots Input complex slots
 * @param slotCount Number of complex slots
 * @param degree Degree of the polynomial (must be >= 2 * slotCount)
 * @param scalingFactor Scaling factor for encoding
 * @param modulus Modulus for coefficient reduction
 * @param coeffsOut Output polynomial coefficients
 */
void CKKSEncodeComplexSlotsToPolynomial(const CKKSComplex *slots, size_t slotCount,
                                        size_t degree, double scalingFactor,
                                        uint64_t modulus, uint64_t *coeffsOut);


void CKKSComplexSlotsToFloats(const CKKSComplex *slots, size_t slotCount, double *output,
                                     size_t outputCount);


void CKKSEncodeFloatsToPolynomial(const double *input, size_t inputCount, 
                                  size_t degree, double scalingFactor, 
                                  uint64_t modulus, uint64_t *coeffsOut);


void CKKSDecodePolynomialToFloats(const uint64_t *coeffs, size_t degree, double scalingFactor,
                                         uint64_t modulus, double *output, size_t outputCount);
                                    

void CKKSDecodePolynomialToComplex(const uint64_t *coeffs, size_t degree, double scalingFactor,
                                          uint64_t modulus, CKKSComplex *slotsOut, size_t slotCount);
#ifdef __cplusplus
}
#endif

#endif // CKKS_ENCODE_H