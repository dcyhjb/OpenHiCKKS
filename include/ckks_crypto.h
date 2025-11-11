#ifndef CKKS_CRYPTO_H
#define CKKS_CRYPTO_H
#define CKKS_NOISE_SIGMA 0.2
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

#include "ckks_polynomial.h"
#include "ckks_encode.h"

#ifdef __cplusplus
extern "C" {
#endif



typedef struct {
    CKKSPolynomial secret;
    CKKSPolynomial publicA;
    CKKSPolynomial publicB;
    CKKSPolynomial evalA;
    CKKSPolynomial evalB;
} CKKSKeyPair;

typedef struct {
    CKKSPolynomial c0;
    CKKSPolynomial c1;
} CKKSCiphertext;

// typedef struct {
//     uint64_t state;
// } CKKSPrng;

// void CKKSPrngInit(CKKSPrng *prng, uint64_t seed);
 
// uint64_t CKKSPrngNext(CKKSPrng *prng);

// double CKKSPrngNextUnit(CKKSPrng *prng);

// double CKKSPrngSampleGaussian(CKKSPrng *prng, double sigma);


// uint64_t CKKSPrngNextMod(CKKSPrng *prng, uint64_t modulus);

int CKKSKeyPairInit(CKKSKeyPair *keyPair, size_t degree, uint64_t modulus);

void CKKSKeyPairFree(CKKSKeyPair *keyPair);
int CKKSCiphertextInit(CKKSCiphertext *ciphertext, size_t degree, uint64_t modulus);


void CKKSCiphertextFree(CKKSCiphertext *ciphertext);
// int CKKSKeyGen(CKKSKeyPair *keyPair, CKKSPrng *prng);
int CKKSKeyGen(CKKSKeyPair *keyPair);
int CKKSCiphertextAdd(const CKKSCiphertext *a, const CKKSCiphertext *b, CKKSCiphertext *result);
int CKKSCiphertextSub(const CKKSCiphertext *a, const CKKSCiphertext *b, CKKSCiphertext *result);
// int CKKSEncrypt(const CKKSKeyPair *keyPair, const CKKSPolynomial *plaintext, CKKSPrng *prng, CKKSCiphertext *ciphertext);
int CKKSEncrypt(const CKKSKeyPair *keyPair, const CKKSPolynomial *plaintext, CKKSCiphertext *ciphertext);


int CKKSCiphertextMul(const CKKSKeyPair *keyPair, const CKKSCiphertext *lhs, const CKKSCiphertext *rhs,
                             double scalingFactor, CKKSCiphertext *result);
int CKKSDecrypt(const CKKSKeyPair *keyPair, const CKKSCiphertext *ciphertext, CKKSPolynomial *decoded);

int64_t CKKSMaxCenteredDifference(const CKKSPolynomial *a, const CKKSPolynomial *b);

#ifdef __cplusplus
}
#endif

#endif /* CKKS_CRYPTO_H */