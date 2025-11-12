#!/bin/bash

# OpenHiCKKS 多规模性能测试脚本
# 
# 本脚本自动运行性能测试并保存结果，便于性能对比和分析

echo "╔════════════════════════════════════════════════════════════╗"
echo "║      OpenHiCKKS 多规模性能测试脚本                        ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# 检查是否在正确的目录
if [ ! -d "build" ]; then
    echo "❌ 错误：找不到 build 目录"
    echo "请先运行："
    echo "  mkdir -p build"
    echo "  cd build"
    echo "  cmake .."
    echo "  make"
    exit 1
fi

# 检查可执行文件是否存在
if [ ! -f "build/test_benchmark" ]; then
    echo "❌ 错误：找不到 test_benchmark 可执行文件"
    echo "请先编译项目："
    echo "  cd build"
    echo "  cmake .."
    echo "  make"
    exit 1
fi

# 创建结果目录
RESULTS_DIR="benchmark_results"
mkdir -p "$RESULTS_DIR"

# 生成时间戳
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
RESULT_FILE="$RESULTS_DIR/benchmark_${TIMESTAMP}.txt"

echo "📊 开始运行多规模性能测试..."
echo "测试规模：512, 1024, 2048, 4096"
echo ""

# 运行测试并保存结果
./build/test_benchmark | tee "$RESULT_FILE"

# 检查测试是否成功
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo ""
    echo "✅ 测试完成！结果已保存到："
    echo "   $RESULT_FILE"
    echo ""
    echo "💡 提示："
    echo "   - 查看结果：cat $RESULT_FILE"
    echo "   - 对比结果：diff $RESULTS_DIR/benchmark_<时间戳1>.txt $RESULTS_DIR/benchmark_<时间戳2>.txt"
    echo ""
else
    echo ""
    echo "❌ 测试失败！请检查错误信息。"
    exit 1
fi

