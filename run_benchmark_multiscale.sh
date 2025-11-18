#!/bin/bash

# OpenHiCKKS Multi-scale Performance Benchmark Script
# 
# This script automatically runs performance tests and saves results for easy comparison and analysis

echo "╔════════════════════════════════════════════════════════════╗"
echo "║      OpenHiCKKS Multi-scale Performance Benchmark         ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Check if in the correct directory
if [ ! -d "build" ]; then
    echo "Error: build directory not found"
    echo "Please run first:"
    echo "  mkdir -p build"
    echo "  cd build"
    echo "  cmake .."
    echo "  make"
    exit 1
fi

# Check if executable exists
if [ ! -f "build/test_benchmark" ]; then
    echo "Error: test_benchmark executable not found"
    echo "Please build the project first:"
    echo "  cd build"
    echo "  cmake .."
    echo "  make"
    exit 1
fi

# Create results directory
RESULTS_DIR="benchmark_results"
mkdir -p "$RESULTS_DIR"

# Generate timestamp
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
RESULT_FILE="$RESULTS_DIR/benchmark_${TIMESTAMP}.txt"

echo "Starting multi-scale performance benchmark..."
echo "Test scales: 512, 1024, 2048, 4096"
echo ""

# Run tests and save results
./build/test_benchmark | tee "$RESULT_FILE"

# Check if tests succeeded
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo ""
    echo "Tests completed! Results saved to:"
    echo "   $RESULT_FILE"
    echo ""
    echo "Tips:"
    echo "   - View results: cat $RESULT_FILE"
    echo "   - Compare results: diff $RESULTS_DIR/benchmark_<timestamp1>.txt $RESULTS_DIR/benchmark_<timestamp2>.txt"
    echo ""
else
    echo ""
    echo "Tests failed! Please check error messages."
    exit 1
fi

