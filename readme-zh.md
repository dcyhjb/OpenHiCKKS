# OpenHiCKKS

一个用 C 语言实现的 CKKS（Cheon-Kim-Kim-Song）同态加密方案库，用于对加密的实数/复数进行近似算术运算。

## 特性

- **完整的 CKKS 方案实现**：密钥生成、加密、解密和同态运算
- **高效的多项式运算**：使用数论变换（NTT）加速多项式乘法
- **编码支持**：支持实数和复数的编码与解码
- **同态运算**：支持加密数据的加法、减法和乘法运算
- **完整的测试**：包含单元测试和性能基准测试

## 架构

库由三个核心模块组成：

### 1. CKKS 编码模块 (`ckks_encode`)
- 将浮点数转换为复数槽位
- 使用逆 FFT 将复数槽位编码为多项式系数
- 将多项式系数解码回浮点数

### 2. CKKS 多项式模块 (`ckks_polynomial`)
- 多项式环算术运算（加法、减法、乘法）
- 基于 NTT 的快速多项式乘法
- 模运算工具函数
- 在环 R = Z[X]/(X^N + 1) 上的运算

### 3. CKKS 密码模块 (`ckks_crypto`)
- 密钥对生成
- 加密和解密
- 同态密文运算（加、减、乘）
- 噪声管理

## 依赖要求

- CMake 3.12 或更高版本
- 支持 C11 的 C 编译器
- OpenHiTLS 库：
  - `hitls_crypto`
  - `hitls_bsl`
  - `boundscheck`

OpenHiTLS 库应安装在 `/usr/local/lib`，头文件在 `/usr/local/include`。

## 编译

```bash
mkdir -p build
cd build
cmake ..
make
```

这将构建：
- 静态库：`libckks_encode.a`、`libckks_polynomial.a`、`libckks_crypto.a`
- 测试可执行文件：`test_main`、`test_benchmark`

## 使用示例

### 基本示例

```c
#include "ckks_crypto.h"
#include "ckks_encode.h"

// 初始化参数
const size_t degree = 1024;
const uint64_t modulus = 65537;
const double scalingFactor = 4096.0;

// 生成密钥对
CKKSKeyPair keyPair;
CKKSKeyPairInit(&keyPair, degree, modulus);
CKKSKeyGen(&keyPair);

// 编码明文
double plaintext[] = {1.5, 2.3, -0.7, 4.1};
CKKSPolynomial plaintextPoly;
CKKSPolynomialInit(&plaintextPoly, degree, modulus);
uint64_t encodedCoeffs[degree];
CKKSEncodeFloatsToPolynomial(plaintext, 4, degree, scalingFactor, modulus, encodedCoeffs);
CKKSPolynomialSet(&plaintextPoly, encodedCoeffs, degree);

// 加密
CKKSCiphertext ciphertext;
CKKSCiphertextInit(&ciphertext, degree, modulus);
CKKSEncrypt(&keyPair, &plaintextPoly, &ciphertext);

// 解密和解码
CKKSPolynomial decrypted;
CKKSPolynomialInit(&decrypted, degree, modulus);
CKKSDecrypt(&keyPair, &ciphertext, &decrypted);
double decoded[4];
CKKSDecodePolynomialToFloats(decrypted.coeffs, degree, scalingFactor, modulus, decoded, 4);

// 清理资源
CKKSCiphertextFree(&ciphertext);
CKKSKeyPairFree(&keyPair);
```

## 测试

运行主测试套件：
```bash
./build/test_main
```

运行性能基准测试：
```bash
./build/test_benchmark
```

或使用基准测试脚本进行多规模测试：
```bash
./run_benchmark_multiscale.sh
```

## API 文档

### 密钥管理
- `CKKSKeyPairInit()`: 初始化密钥对结构
- `CKKSKeyGen()`: 生成公钥和私钥
- `CKKSKeyPairFree()`: 释放密钥对资源

### 加密/解密
- `CKKSEncrypt()`: 加密明文多项式
- `CKKSDecrypt()`: 解密密文

### 同态运算
- `CKKSCiphertextAdd()`: 密文加法
- `CKKSCiphertextSub()`: 密文减法
- `CKKSCiphertextMul()`: 密文乘法

### 编码/解码
- `CKKSEncodeFloatsToPolynomial()`: 将浮点数编码为多项式
- `CKKSDecodePolynomialToFloats()`: 将多项式解码为浮点数
- `CKKSEncodeComplexSlotsToPolynomial()`: 将复数槽位编码为多项式
- `CKKSDecodePolynomialToComplex()`: 将多项式解码为复数槽位

### 多项式运算
- `CKKSAdd()`: 多项式加法
- `CKKSSub()`: 多项式减法
- `CKKSMul()`: 多项式乘法（NTT 加速）
- `CKKSNTT()`: 数论变换
- `CKKSInverseNTT()`: 逆数论变换

详细 API 文档请参见 `include/` 目录中的头文件。
