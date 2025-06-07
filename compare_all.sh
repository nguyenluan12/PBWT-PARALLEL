#!/bin/bash

# Compile both implementations
g++-14 -O3 -std=c++17 -I/opt/homebrew/include -L/opt/homebrew/lib -o pbwt_single pbwt-single.cpp util.cpp -lboost_iostreams -lz
g++-14 -O3 -fopenmp -std=c++17 -I/opt/homebrew/include -L/opt/homebrew/lib -o pbwt_haps pbwt-haps.cpp util.cpp -lboost_iostreams -lz

# Define test files
test_files=(
    "10SitesX10Haps.hap"
    "20SitesX10Haps.hap"
    "30SitesX20Haps.hap"
    "33SitesX22Haps.hap"
    "100SitesX50Haps.hap"
)

# Compare outputs for each test file
for file in "${test_files[@]}"; do
    echo "===== Testing $file ====="
    
    # Run both implementations
    ./pbwt_single "$file" > "single_output_$file.txt"
    ./pbwt_haps "$file" > "parallel_output_$file.txt"
    
    # Extract just the maximal matches section from both outputs
    grep -A1000 "========== Maximal Matches ==========" "single_output_$file.txt" > "single_maximal_$file.txt"
    grep -A1000 "========== Maximal Matches ==========" "parallel_output_$file.txt" > "parallel_maximal_$file.txt"
    
    # Compare the outputs
    echo "Comparing outputs for $file..."
    diff "single_maximal_$file.txt" "parallel_maximal_$file.txt" > "diff_$file.txt"
    
    # Check if there are differences
    if [ -s "diff_$file.txt" ]; then
        echo "DIFFERENCES FOUND in $file"
        echo "First 10 differences:"
        head -n 10 "diff_$file.txt"
    else
        echo "SUCCESS: Outputs match for $file"
    fi
    
    echo ""
done

echo "All tests completed!"
