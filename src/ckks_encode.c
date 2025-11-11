#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "ckks_polynomial.h"  // For CKKSModReduceInt64
#include "ckks_encode.h"  // Include the encode header


void CKKSCopyFloatsToComplexSlots(const double *input, size_t inputCount, CKKSComplex *slots,
                                         size_t slotCount)
{
    if (slots == NULL) {
        return;
    }
    for (size_t i = 0; i < slotCount; ++i) {
        if (input != NULL && i < inputCount) {
            slots[i].real = input[i];
            slots[i].imag = 0.0;
        } else {
            slots[i].real = 0.0;
            slots[i].imag = 0.0;
        }
    }
}



void CKKSEncodeComplexSlotsToPolynomial(const CKKSComplex *slots, size_t slotCount, size_t degree,
                                               double scalingFactor, uint64_t modulus, uint64_t *coeffsOut)
{
    if (slots == NULL || coeffsOut == NULL || degree == 0) {
        return;
    }

    size_t slotCapacity = degree / 2U;
    if (slotCount > slotCapacity) {
        slotCount = slotCapacity;
    }

    for (size_t i = 0; i < degree; ++i) {
        coeffsOut[i] = 0;
    }

    for (size_t i = 0; i < slotCount; ++i) {
        long double scaledReal = (long double)slots[i].real * (long double)scalingFactor;
        long double scaledImag = (long double)slots[i].imag * (long double)scalingFactor;
        int64_t roundedReal = (int64_t)llround((double)scaledReal);
        int64_t roundedImag = (int64_t)llround((double)scaledImag);

        coeffsOut[i] = CKKSModReduceInt64(roundedReal, modulus);

        size_t imagIndex = i + slotCapacity;
        if (imagIndex < degree) {
            coeffsOut[imagIndex] = CKKSModReduceInt64(roundedImag, modulus);
        }
    }
}

void CKKSEncodeFloatsToPolynomial(const double *input, size_t inputCount, 
                                  size_t degree, double scalingFactor, 
                                  uint64_t modulus, uint64_t *coeffsOut)
{
    if (coeffsOut == NULL || degree == 0) {
        return;
    }

    // 计算可用的复数槽数量
    size_t slotCapacity = degree / 2;
    size_t slotCount = (inputCount > slotCapacity) ? slotCapacity : inputCount;
    
    // 分配复数槽数组
    CKKSComplex *slots = (CKKSComplex *)malloc(slotCapacity * sizeof(CKKSComplex));
    if (slots == NULL) {
        return; // 内存分配失败
    }
    
    // 第一步：将浮点数复制到复数槽
    CKKSCopyFloatsToComplexSlots(input, inputCount, slots, slotCapacity);
    
    // 第二步：将复数槽编码为多项式
    CKKSEncodeComplexSlotsToPolynomial(slots, slotCapacity, degree, 
                                      scalingFactor, modulus, coeffsOut);
    
    // 清理内存
    free(slots);
}

void CKKSDecodePolynomialToComplex(const uint64_t *coeffs, size_t degree, double scalingFactor,
                                          uint64_t modulus, CKKSComplex *slotsOut, size_t slotCount)
{
    if (coeffs == NULL || slotsOut == NULL || degree == 0 || scalingFactor == 0.0) {
        return;
    }

    size_t slotCapacity = degree / 2U;
    if (slotCount > slotCapacity) {
        slotCount = slotCapacity;
    }

    uint64_t halfModulus = modulus / 2U;
    for (size_t i = 0; i < slotCount; ++i) {
        uint64_t coeffReal = coeffs[i] % modulus;
        int64_t signedReal = (coeffReal > halfModulus) ? (int64_t)coeffReal - (int64_t)modulus : (int64_t)coeffReal;
        slotsOut[i].real = (double)((long double)signedReal / (long double)scalingFactor);

        size_t imagIndex = i + slotCapacity;
        if (imagIndex < degree) {
            uint64_t coeffImag = coeffs[imagIndex] % modulus;
            int64_t signedImag =
                (coeffImag > halfModulus) ? (int64_t)coeffImag - (int64_t)modulus : (int64_t)coeffImag;
            slotsOut[i].imag = (double)((long double)signedImag / (long double)scalingFactor);
        } else {
            slotsOut[i].imag = 0.0;
        }
    }
}

void CKKSComplexSlotsToFloats(const CKKSComplex *slots, size_t slotCount, double *output,
                                     size_t outputCount)
{
    if (slots == NULL || output == NULL) {
        return;
    }

    size_t limit = slotCount < outputCount ? slotCount : outputCount;
    for (size_t i = 0; i < limit; ++i) {
        output[i] = slots[i].real;
    }
    for (size_t i = limit; i < outputCount; ++i) {
        output[i] = 0.0;
    }
}

 void CKKSDecodePolynomialToFloats(const uint64_t *coeffs, size_t degree, double scalingFactor,
                                         uint64_t modulus, double *output, size_t outputCount)
{
    if (coeffs == NULL || output == NULL || degree == 0) {
        return;
    }

    size_t slotCount = degree / 2U;
    CKKSComplex *decodedSlots = (CKKSComplex *)calloc(slotCount, sizeof(CKKSComplex));
    if (decodedSlots == NULL) {
        return;
    }

    CKKSDecodePolynomialToComplex(coeffs, degree, scalingFactor, modulus, decodedSlots, slotCount);
    CKKSComplexSlotsToFloats(decodedSlots, slotCount, output, outputCount);

    free(decodedSlots);
}

