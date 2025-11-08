#ifndef CKKS_POLYNOMIAL_H
#define CKKS_POLYNOMIAL_H

#include <stdint.h>
#include <stddef.h>

/**
 * @file ckks_polynomial.h
 * @brief CKKS同态加密多项式环运算库
 * 
 * 这个库实现了CKKS同态加密方案中所需的多项式环运算，包括：
 * - 多项式加法、减法、乘法
 * - 标量运算
 * - 使用数论变换(NTT)加速的多项式乘法
 * - 模运算支持
 * 
 * 所有运算在环 R = Z[X]/(X^N + 1) 上进行，其中N是2的幂次方。
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief CKKS多项式结构体
 * 
 * 表示在模数modulus定义的有界域上的多项式，
 * 形式为: coeffs[0] + coeffs[1]*x + ... + coeffs[degree-1]*x^(degree-1)
 */
typedef struct {
    size_t degree;      /**< 多项式的度数 */
    uint64_t modulus;   /**< 模数，定义系数所在的有限域 */
    uint64_t *coeffs;   /**< 多项式系数数组，长度为degree */
} CKKSPolynomial;

/* ====== 内存管理函数 ====== */

/**
 * @brief 初始化多项式
 * 
 * 分配多项式系数数组并设置初始参数。
 * 
 * @param poly 要初始化的多项式指针
 * @param degree 多项式的度数
 * @param modulus 模数
 * @return 成功返回0，失败返回-1
 */
 int CKKSPolynomialInit(CKKSPolynomial *poly, size_t degree, uint64_t modulus);

/**
 * @brief 释放多项式资源
 * 
 * 释放多项式系数数组并将结构体重置为初始状态。
 * 
 * @param poly 要释放的多项式指针
 */
 void CKKSPolynomialFree(CKKSPolynomial *poly);

/**
 * @brief 设置多项式系数
 * 
 * 从外部数组复制系数到多项式结构中，并自动进行模约简。
 * 
 * @param poly 目标多项式
 * @param coeffs 系数数组
 * @param count 系数数量（必须等于多项式的度数）
 * @return 成功返回0，失败返回-1
 */
int CKKSPolynomialSet(CKKSPolynomial *poly, const uint64_t *coeffs, size_t count);

/* ====== 基本算术运算 ====== */

/**
 * @brief 多项式加法
 * 
 * 计算 result = a + b，在环 R = Z[X]/(X^N + 1) 上进行。
 * 
 * @param a 第一个多项式
 * @param b 第二个多项式
 * @param result 存储结果的多项式
 * @return 成功返回0，失败返回-1
 */
 int CKKSPolynomialAdd(const CKKSPolynomial *a, const CKKSPolynomial *b, CKKSPolynomial *result);

/**
 * @brief 多项式减法
 * 
 * 计算 result = a - b，在环 R = Z[X]/(X^N + 1) 上进行。
 * 
 * @param a 被减多项式
 * @param b 减数多项式
 * @param result 存储结果的多项式
 * @return 成功返回0，失败返回-1
 */

int CKKSPolynomialAddInplace(CKKSPolynomial *accumulator, const CKKSPolynomial *summand);


 int CKKSPolynomialSub(const CKKSPolynomial *a, const CKKSPolynomial *b, CKKSPolynomial *result);

/**
 * @brief 多项式乘法（使用NTT加速）
 * 
 * 计算 result = a * b mod (X^degree + 1)，使用数论变换进行加速。
 * 这是CKKS方案中最复杂的运算，使用NTT将时间复杂度从O(N^2)降低到O(N log N)。
 * 
 * @param a 第一个多项式
 * @param b 第二个多项式
 * @param result 存储结果的多项式
 * @return 成功返回0，失败返回-1
 */
 int CKKSMul(const CKKSPolynomial *a, const CKKSPolynomial *b, CKKSPolynomial *result);

/**
 * @brief 标量除法
 * 
 * 计算 result = poly / scalar，通过对每个系数乘以标量的模逆元实现。
 * 
 * @param poly 被除多项式
 * @param scalar 标量除数（不能为0）
 * @param result 存储结果的多项式
 * @return 成功返回0，失败返回-1（标量没有逆元或参数错误）
 */
int CKKSScalarDivide(const CKKSPolynomial *poly, uint64_t scalar, CKKSPolynomial *result);

/* ====== 辅助函数 ====== */

/**
 * @brief 检查多项式系数是否匹配期望值
 * 
 * 用于测试和验证运算结果的正确性。
 * 
 * @param poly 要检查的多项式
 * @param expected 期望的系数数组
 * @return 匹配返回1，不匹配返回0
 */
int CKKSMatchCoefficients(const CKKSPolynomial *poly, const uint64_t *expected);

/**
 * @brief 朴素多项式乘法（用于验证）
 * 
 * 使用直接卷积方法计算多项式乘法，用于验证NTT实现的正确性。
 * 注意：此函数仅用于测试，性能较低。
 * 
 * @param a 第一个系数数组
 * @param b 第二个系数数组
 * @param degree 多项式度数
 * @param modulus 模数
 * @param result 存储结果的数组（长度至少为degree）
 */
void CKKSNaiveMulReduce(const uint64_t *a, const uint64_t *b, size_t degree, 
                        uint64_t modulus, uint64_t *result);

/* ====== 底层数学函数 ====== */

/**
 * @brief 模约简函数（64位）
 * 
 * 将有符号64位整数约简到[0, modulus-1]范围内。
 * 
 * @param value 要约简的值
 * @param modulus 模数
 * @return 约简后的值
 */
 uint64_t CKKSModReduceInt64(int64_t value, uint64_t modulus);

/**
 * @brief 模约简函数（128位）
 * 
 * 处理大数乘法结果的模运算，使用128位整数中间结果。
 * 
 * @param value 要约简的值（128位）
 * @param modulus 模数
 * @return 约简后的值
 */
uint64_t CKKSModReduceInt128(__int128 value, uint64_t modulus);

/**
 * @brief 模幂运算
 * 
 * 使用快速幂算法计算 base^exponent mod modulus。
 * 
 * @param base 底数
 * @param exponent 指数
 * @param modulus 模数
 * @return 计算结果
 */
uint64_t CKKSPolynomialPowMod(uint64_t base, uint64_t exponent, uint64_t modulus);

/**
 * @brief 模乘法
 * 
 * 计算 a * b mod modulus，使用128位中间结果避免溢出。
 * 
 * @param a 第一个操作数
 * @param b 第二个操作数
 * @param modulus 模数
 * @return 计算结果
 */
 uint64_t CKKSPolynomialMulMod(uint64_t a, uint64_t b, uint64_t modulus);

/**
 * @brief 模逆元计算
 * 
 * 使用扩展欧几里得算法计算value在模modulus下的逆元。
 * 
 * @param value 要求逆元的值
 * @param modulus 模数
 * @return 逆元，如果不存在则返回0
 */
uint64_t CKKSPolynomialModInverse(uint64_t value, uint64_t modulus);

/* ====== 数论变换函数 ====== */

/**
 * @brief 数论变换 (NTT)
 * 
 * 将多项式从时域转换到频域，用于加速多项式乘法。
 * 
 * @param data 输入/输出数据数组
 * @param length 数据长度（必须是2的幂次）
 * @param modulus 模数
 * @param primitiveRoot 原根
 */
 void CKKSNTT(uint64_t *data, size_t length, uint64_t modulus, uint64_t primitiveRoot);

/**
 * @brief 逆数论变换 (Inverse NTT)
 * 
 * 将多项式从频域转换回时域。
 * 
 * @param data 输入/输出数据数组
 * @param length 数据长度（必须是2的幂次）
 * @param modulus 模数
 * @param primitiveRoot 原根
 */
void CKKSInverseNTT(uint64_t *data, size_t length, uint64_t modulus, uint64_t primitiveRoot);

void CKKSPolynomialPrint(const CKKSPolynomial *poly, const char *label, size_t limit);


#ifdef __cplusplus
}
#endif

#endif /* CKKS_POLYNOMIAL_H */