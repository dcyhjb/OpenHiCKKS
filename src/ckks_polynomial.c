#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ckks_polynomial.h"

// /* 模逆元函数声明 */
// // static uint64_t CKKSPolynomialModInverse(uint64_t value, uint64_t modulus);

// /* 
//  * 64位整数模约简函数
//  * 将64位有符号整数约简到[0, modulus-1]范围内
//  */
uint64_t CKKSModReduceInt64(int64_t value, uint64_t modulus)
{
    int64_t reduced = value % (int64_t)modulus;
    if (reduced < 0) {
        reduced += (int64_t)modulus;  // 处理负数情况
    }
    return (uint64_t)reduced;
}

/* 
 * 128位整数模约简函数
 * 处理大数乘法结果的模运算
 */
uint64_t CKKSModReduceInt128(__int128 value, uint64_t modulus)
{
    __int128 reduced = value % (__int128)modulus;
    if (reduced < 0) {
        reduced += (__int128)modulus;  // 处理负数情况
    }
    return (uint64_t)reduced;
}

/* 
 * 初始化多项式
 * 分配系数数组内存并设置度数和模数
 */
 int CKKSPolynomialInit(CKKSPolynomial *poly, size_t degree, uint64_t modulus)
{
    if (poly == NULL || degree == 0 || modulus == 0) {
        return -1;  // 参数检查
    }
    poly->degree = degree;
    poly->modulus = modulus;
    poly->coeffs = (uint64_t *)calloc(degree, sizeof(uint64_t));  // 分配并初始化为0
    return poly->coeffs == NULL ? -1 : 0;  // 内存分配检查
}

/* 释放多项式资源 */
void CKKSPolynomialFree(CKKSPolynomial *poly)
{
    if (poly == NULL) {
        return;
    }
    free(poly->coeffs);      // 释放系数数组
    poly->coeffs = NULL;     // 避免野指针
    poly->degree = 0;
    poly->modulus = 0;
}

/* 
 * 设置多项式系数
 * 将外部系数数组复制到多项式结构中，并进行模约简
 */
 int CKKSPolynomialSet(CKKSPolynomial *poly, const uint64_t *coeffs, size_t count)
{
    if (poly == NULL || coeffs == NULL || count != poly->degree) {
        return -1;  // 参数检查
    }
    (void)memcpy(poly->coeffs, coeffs, poly->degree * sizeof(uint64_t));
    // 确保所有系数都在模数范围内
    for (size_t i = 0; i < poly->degree; ++i) {
        poly->coeffs[i] %= poly->modulus;
    }
    return 0;
}

/* 
 * 多项式加法
 * 对应系数相加并模约简
 */
 int CKKSPolynomialAdd(const CKKSPolynomial *a, const CKKSPolynomial *b, CKKSPolynomial *result)
{
    if (a == NULL || b == NULL || result == NULL) {
        return -1;  // 空指针检查
    }
    // 检查度数和模数是否匹配
    if (a->degree != b->degree || a->degree != result->degree || a->modulus != b->modulus ||
        a->modulus != result->modulus) {
        return -1;
    }
    // 逐系数相加
    for (size_t i = 0; i < a->degree; ++i) {
        int64_t value = (int64_t)a->coeffs[i] + (int64_t)b->coeffs[i];
        result->coeffs[i] = CKKSModReduceInt64(value, a->modulus);
    }
    return 0;
}

/* 
 * 多项式减法
 * 对应系数相减并模约简
 */
 int CKKSPolynomialSub(const CKKSPolynomial *a, const CKKSPolynomial *b, CKKSPolynomial *result)
{
    if (a == NULL || b == NULL || result == NULL) {
        return -1;  // 空指针检查
    }
    // 检查度数和模数是否匹配
    if (a->degree != b->degree || a->degree != result->degree || a->modulus != b->modulus ||
        a->modulus != result->modulus) {
        return -1;
    }
    // 逐系数相减
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
 * 模幂运算 (快速幂算法)
 * 计算 base^exponent mod modulus
 */
 uint64_t CKKSPolynomialPowMod(uint64_t base, uint64_t exponent, uint64_t modulus)
{
    uint64_t result = 1;
    base %= modulus;  // 初始模约简
    while (exponent > 0) {
        if (exponent & 1U) {  // 如果当前位为1
            result = CKKSModReduceInt128((__int128)result * base, modulus);
        }
        base = CKKSModReduceInt128((__int128)base * base, modulus);  // 平方
        exponent >>= 1U;  // 移到下一位
    }
    return result;
}

/* 模乘法：计算 a * b mod modulus */
 uint64_t CKKSPolynomialMulMod(uint64_t a, uint64_t b, uint64_t modulus)
{
    return CKKSModReduceInt128((__int128)a * (__int128)b, modulus);
}

/* 
 * 位反转函数
 * 用于NTT中的位反转置换
 */
size_t CKKSReverseBits(size_t value, unsigned int bitCount)
{
    size_t reversed = 0;
    for (unsigned int i = 0; i < bitCount; ++i) {
        reversed <<= 1U;
        reversed |= (value & 1U);  // 取最低位
        value >>= 1U;  // 右移
    }
    return reversed;
}

/* 
 * 位反转置换
 * 将数组元素按位反转顺序重新排列
 */
 void CKKSBitReverse(uint64_t *data, size_t length)
{
    unsigned int bitCount = 0;
    size_t temp = length;
    // 计算需要的位数
    while (temp > 1) {
        ++bitCount;
        temp >>= 1U;
    }
    // 执行位反转置换
    for (size_t i = 0; i < length; ++i) {
        size_t j = CKKSReverseBits(i, bitCount);
        if (j > i) {  // 避免重复交换
            uint64_t swap = data[i];
            data[i] = data[j];
            data[j] = swap;
        }
    }
}

/* 
 * 数论变换 (NTT)
 * 将多项式从时域转换到频域，加速多项式乘法
 */
void CKKSNTT(uint64_t *data, size_t length, uint64_t modulus, uint64_t primitiveRoot)
{
    CKKSBitReverse(data, length);  // 位反转置换
    
    // 蝴蝶操作
    for (size_t len = 2; len <= length; len <<= 1U) {
        size_t half = len >> 1U;
        // 计算旋转因子
        uint64_t rootStep = CKKSPolynomialPowMod(primitiveRoot, length / len, modulus);
        
        for (size_t i = 0; i < length; i += len) {
            uint64_t w = 1;  // 初始旋转因子
            for (size_t j = 0; j < half; ++j) {
                uint64_t u = data[i + j];
                uint64_t v = CKKSPolynomialMulMod(data[i + j + half], w, modulus);
                
                // 蝴蝶操作核心计算
                uint64_t add = u + v;
                if (add >= modulus) {
                    add -= modulus;  // 模约简
                }
                uint64_t sub = u >= v ? u - v : u + modulus - v;
                
                data[i + j] = add;
                data[i + j + half] = sub;
                
                w = CKKSPolynomialMulMod(w, rootStep, modulus);  // 更新旋转因子
            }
        }
    }
}

/* 
 * 逆数论变换 (Inverse NTT)
 * 将多项式从频域转换回时域
 */
 void CKKSInverseNTT(uint64_t *data, size_t length, uint64_t modulus, uint64_t primitiveRoot)
{
    CKKSBitReverse(data, length);  // 位反转置换
    
    // 蝴蝶操作（与NTT类似但方向相反）
    for (size_t len = 2; len <= length; len <<= 1U) {
        size_t half = len >> 1U;
        // 计算旋转因子
        uint64_t rootStep = CKKSPolynomialPowMod(primitiveRoot, length / len, modulus);
        
        for (size_t i = 0; i < length; i += len) {
            uint64_t w = 1;  // 初始旋转因子
            for (size_t j = 0; j < half; ++j) {
                uint64_t u = data[i + j];
                uint64_t v = CKKSPolynomialMulMod(data[i + j + half], w, modulus);
                
                // 蝴蝶操作核心计算
                data[i + j] = u + v;
                if (data[i + j] >= modulus) {
                    data[i + j] -= modulus;  // 模约简
                }
                uint64_t sub = u >= v ? u - v : u + modulus - v;
                data[i + j + half] = sub;
                
                w = CKKSPolynomialMulMod(w, rootStep, modulus);  // 更新旋转因子
            }
        }
    }
    
    // 乘以长度的逆元完成逆变换
    uint64_t invLength = CKKSPolynomialModInverse((uint64_t)length, modulus);
    for (size_t i = 0; i < length; ++i) {
        data[i] = CKKSPolynomialMulMod(data[i], invLength, modulus);
    }
}

/* 
 * 多项式乘法 (使用NTT加速)
 * 计算 a * b mod (x^degree + 1)
 */
 int CKKSMul(const CKKSPolynomial *a, const CKKSPolynomial *b, CKKSPolynomial *result)
{
    if (a == NULL || b == NULL || result == NULL) {
        return -1;  // 空指针检查
    }
    // 检查度数和模数是否匹配
    if (a->degree != b->degree || a->degree != result->degree || a->modulus != b->modulus ||
        a->modulus != result->modulus) {
        return -1;
    }

    const size_t degree = a->degree;
    // 计算NTT需要的变换长度（2的幂次）
    size_t transformSize = 1;
    while (transformSize < degree * 2U) {
        transformSize <<= 1U;
    }

    // 分配NTT变换的缓冲区
    uint64_t *fa = (uint64_t *)calloc(transformSize, sizeof(uint64_t));
    uint64_t *fb = (uint64_t *)calloc(transformSize, sizeof(uint64_t));
    if (fa == NULL || fb == NULL) {
        free(fa);
        free(fb);
        return -1;  // 内存分配失败
    }

    // 复制多项式系数到缓冲区
    for (size_t i = 0; i < degree; ++i) {
        fa[i] = a->coeffs[i];
        fb[i] = b->coeffs[i];
    }

    // 设置NTT参数
    const uint64_t primitiveRoot = 3; /* 模65537的原根 */
    uint64_t root = CKKSPolynomialPowMod(primitiveRoot, (a->modulus - 1) / transformSize, a->modulus);
    
    // 执行NTT变换
    CKKSNTT(fa, transformSize, a->modulus, root);
    CKKSNTT(fb, transformSize, a->modulus, root);

    // 频域点乘
    for (size_t i = 0; i < transformSize; ++i) {
        fa[i] = CKKSPolynomialMulMod(fa[i], fb[i], a->modulus);
    }

    // 逆NTT变换
    uint64_t invRoot = CKKSPolynomialModInverse(root, a->modulus);
    CKKSInverseNTT(fa, transformSize, a->modulus, invRoot);

    // 复制结果
    for (size_t i = 0; i < degree; ++i) {
        result->coeffs[i] = fa[i];
    }

    // 处理模 x^degree + 1 的约简
    for (size_t i = degree; i < transformSize && i < degree * 2U; ++i) {
        size_t target = i - degree;
        uint64_t value = fa[i];
        if (value != 0) {
            // 执行 x^degree = -1 的约简
            if (result->coeffs[target] >= value) {
                result->coeffs[target] -= value;
            } else {
                result->coeffs[target] = result->coeffs[target] + a->modulus - value;
            }
        }
    }

    // 释放资源
    free(fa);
    free(fb);
    return 0;
}

/* 
 * 模逆元计算 (扩展欧几里得算法)
 * 计算 value 在模 modulus 下的逆元
 */
uint64_t CKKSPolynomialModInverse(uint64_t value, uint64_t modulus)
{
    int64_t t = 0;
    int64_t newT = 1;
    int64_t r = (int64_t)modulus;
    int64_t newR = (int64_t)(value % modulus);

    // 扩展欧几里得算法
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
        return 0; /* 逆元不存在 */
    }

    if (t < 0) {
        t += (int64_t)modulus;  // 调整到正数范围
    }
    return (uint64_t)t;
}

/* 
 * 标量除法
 * 多项式每个系数除以标量
 */
 int CKKSScalarDivide(const CKKSPolynomial *poly, uint64_t scalar, CKKSPolynomial *result)
{
    if (poly == NULL || result == NULL || scalar == 0) {
        return -1;  // 参数检查
    }
    if (poly->degree != result->degree || poly->modulus != result->modulus) {
        return -1;  // 度数和模数检查
    }

    // 计算标量的模逆元
    uint64_t inverse = CKKSPolynomialModInverse(scalar, poly->modulus);
    if (inverse == 0) {
        return -1;  // 逆元不存在
    }

    // 每个系数乘以逆元实现除法
    for (size_t i = 0; i < poly->degree; ++i) {
        __int128 product = (__int128)poly->coeffs[i] * (__int128)inverse;
        result->coeffs[i] = CKKSModReduceInt128(product, poly->modulus);
    }
    return 0;
}

/* 
 * 系数匹配检查
 * 验证多项式系数是否与期望值匹配
 */
 int CKKSMatchCoefficients(const CKKSPolynomial *poly, const uint64_t *expected)
{
    if (poly == NULL || expected == NULL) {
        return 0;  // 空指针检查
    }
    for (size_t i = 0; i < poly->degree; ++i) {
        if (poly->coeffs[i] != expected[i]) {
            return 0;  // 系数不匹配
        }
    }
    return 1;  // 所有系数匹配
}

/* 
 * 朴素多项式乘法（用于验证NTT结果的正确性）
 * 直接计算多项式乘法然后进行模约简
 */
 void CKKSNaiveMulReduce(const uint64_t *a, const uint64_t *b, size_t degree, uint64_t modulus,
                               uint64_t *result)
{
    if (a == NULL || b == NULL || result == NULL) {
        return;  // 空指针检查
    }

    size_t convLength = degree * 2U;
    uint64_t *temp = (uint64_t *)calloc(convLength, sizeof(uint64_t));
    assert(temp != NULL);

    // 直接卷积计算
    for (size_t i = 0; i < degree; ++i) {
        for (size_t j = 0; j < degree; ++j) {
            size_t idx = i + j;
            uint64_t prod = CKKSPolynomialMulMod(a[i], b[j], modulus);
            uint64_t sum = temp[idx] + prod;
            if (sum >= modulus) {
                sum -= modulus;  // 模约简
            }
            temp[idx] = sum;
        }
    }

    // 复制低次项
    for (size_t i = 0; i < degree; ++i) {
        result[i] = temp[i];
    }

    // 模 x^degree + 1 约简
    for (size_t idx = degree; idx < convLength; ++idx) {
        uint64_t value = temp[idx];
        if (value == 0) {
            continue;
        }
        size_t target = idx - degree;
        // 执行 x^degree = -1 的约简
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



