# OpenHiCKKS

A C library implementation of the CKKS (Cheon-Kim-Kim-Song) homomorphic encryption scheme, designed for approximate arithmetic on encrypted real/complex numbers.

## Features

- **Full CKKS Scheme Implementation**: Key generation, encryption, decryption, and homomorphic operations
- **Efficient Polynomial Operations**: Number Theoretic Transform (NTT) accelerated polynomial multiplication
- **Encoding Support**: Encoding and decoding for real and complex numbers
- **Homomorphic Operations**: Support for addition, subtraction, and multiplication on encrypted data
- **Comprehensive Testing**: Unit tests and performance benchmarks included

## Architecture

The library is organized into three core modules:

### 1. CKKS Encode (`ckks_encode`)
- Converts floating-point numbers to complex slots
- Encodes complex slots into polynomial coefficients using inverse FFT
- Decodes polynomial coefficients back to floating-point numbers

### 2. CKKS Polynomial (`ckks_polynomial`)
- Polynomial ring arithmetic operations (addition, subtraction, multiplication)
- NTT-based fast polynomial multiplication
- Modular arithmetic utilities
- Operations on the ring R = Z[X]/(X^N + 1)

### 3. CKKS Crypto (`ckks_crypto`)
- Key pair generation
- Encryption and decryption
- Homomorphic ciphertext operations (add, subtract, multiply)
- Noise management

## Requirements

- CMake 3.12 or higher
- C compiler with C11 support
- OpenHiTLS libraries:
  - `hitls_crypto`
  - `hitls_bsl`
  - `boundscheck`

The OpenHiTLS libraries should be installed in `/usr/local/lib` and headers in `/usr/local/include`.

## Building

```bash
mkdir -p build
cd build
cmake ..
make
```

This will build:
- Static libraries: `libckks_encode.a`, `libckks_polynomial.a`, `libckks_crypto.a`
- Test executables: `test_main`, `test_benchmark`

## Usage

### Basic Example

```c
#include "ckks_crypto.h"
#include "ckks_encode.h"

// Initialize parameters
const size_t degree = 1024;
const uint64_t modulus = 65537;
const double scalingFactor = 4096.0;

// Generate key pair
CKKSKeyPair keyPair;
CKKSKeyPairInit(&keyPair, degree, modulus);
CKKSKeyGen(&keyPair);

// Encode plaintext
double plaintext[] = {1.5, 2.3, -0.7, 4.1};
CKKSPolynomial plaintextPoly;
CKKSPolynomialInit(&plaintextPoly, degree, modulus);
uint64_t encodedCoeffs[degree];
CKKSEncodeFloatsToPolynomial(plaintext, 4, degree, scalingFactor, modulus, encodedCoeffs);
CKKSPolynomialSet(&plaintextPoly, encodedCoeffs, degree);

// Encrypt
CKKSCiphertext ciphertext;
CKKSCiphertextInit(&ciphertext, degree, modulus);
CKKSEncrypt(&keyPair, &plaintextPoly, &ciphertext);

// Decrypt and decode
CKKSPolynomial decrypted;
CKKSPolynomialInit(&decrypted, degree, modulus);
CKKSDecrypt(&keyPair, &ciphertext, &decrypted);
double decoded[4];
CKKSDecodePolynomialToFloats(decrypted.coeffs, degree, scalingFactor, modulus, decoded, 4);

// Cleanup
CKKSCiphertextFree(&ciphertext);
CKKSKeyPairFree(&keyPair);
```

## Testing

Run the main test suite:
```bash
./build/test_main
```

Run performance benchmarks:
```bash
./build/test_benchmark
```

Or use the benchmark script for multi-scale testing:
```bash
./run_benchmark_multiscale.sh
```

## API Documentation

### Key Management
- `CKKSKeyPairInit()`: Initialize a key pair structure
- `CKKSKeyGen()`: Generate public and secret keys
- `CKKSKeyPairFree()`: Free key pair resources

### Encryption/Decryption
- `CKKSEncrypt()`: Encrypt a plaintext polynomial
- `CKKSDecrypt()`: Decrypt a ciphertext

### Homomorphic Operations
- `CKKSCiphertextAdd()`: Add two ciphertexts
- `CKKSCiphertextSub()`: Subtract two ciphertexts
- `CKKSCiphertextMul()`: Multiply two ciphertexts

### Encoding/Decoding
- `CKKSEncodeFloatsToPolynomial()`: Encode floating-point numbers to polynomial
- `CKKSDecodePolynomialToFloats()`: Decode polynomial to floating-point numbers
- `CKKSEncodeComplexSlotsToPolynomial()`: Encode complex slots to polynomial
- `CKKSDecodePolynomialToComplex()`: Decode polynomial to complex slots

### Polynomial Operations
- `CKKSAdd()`: Polynomial addition
- `CKKSSub()`: Polynomial subtraction
- `CKKSMul()`: Polynomial multiplication (NTT-accelerated)
- `CKKSNTT()`: Number Theoretic Transform
- `CKKSInverseNTT()`: Inverse NTT

See header files in `include/` for detailed API documentation.
