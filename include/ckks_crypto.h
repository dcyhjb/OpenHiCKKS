#ifndef CKKS_CRYPTO_H
#define CKKS_CRYPTO_H

#include <stdint.h>
#include <stddef.h>

#include "ckks_polynomial.h"
#include "ckks_encode.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    CKKSPolynomial secret;
    CKKSPolynomial publicA;
    CKKSPolynomial publicB;
} CKKSKeyPair;

typedef struct {
    CKKSPolynomial c0;
    CKKSPolynomial c1;
} CKKSCiphertext;

typedef struct {
    uint64_t state;
} CKKSPrng;


void CKKSPrngInit(CKKSPrng *prng, uint64_t seed);
 
uint64_t CKKSPrngNext(CKKSPrng *prng);

uint64_t CKKSPrngNextMod(CKKSPrng *prng, uint64_t modulus);
int CKKSKeyPairInit(CKKSKeyPair *keyPair, size_t degree, uint64_t modulus);

void CKKSKeyPairFree(CKKSKeyPair *keyPair);
int CKKSCiphertextInit(CKKSCiphertext *ciphertext, size_t degree, uint64_t modulus);
void CKKSCiphertextFree(CKKSCiphertext *ciphertext);
int CKKSKeyGen(CKKSKeyPair *keyPair, CKKSPrng *prng);

int CKKSEncrypt(const CKKSKeyPair *keyPair, const CKKSPolynomial *plaintext, CKKSPrng *prng, CKKSCiphertext *ciphertext);


int CKKSDecrypt(const CKKSKeyPair *keyPair, const CKKSCiphertext *ciphertext, CKKSPolynomial *decoded);



#ifdef __cplusplus
}
#endif

#endif /* CKKS_CRYPTO_H */