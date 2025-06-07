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
    
    # Compare the files
    if diff -q "single_maximal_$file.txt" "parallel_maximal_$file.txt" > /dev/null; then
        echo "MATCH: Outputs for $file match exactly!"
    else
        echo "DIFFERENCES FOUND in $file"
        # Count the number of matches in each file
        single_count=$(grep -c "Position" "single_maximal_$file.txt")
        parallel_count=$(grep -c "Position" "parallel_maximal_$file.txt")
        echo "Single implementation: $single_count matches"
        echo "Parallel implementation: $parallel_count matches"
        
        # Show the first few differences
        echo "First 5 differences:"
        diff -y --suppress-common-lines "single_maximal_$file.txt" "parallel_maximal_$file.txt" | head -5
    fi
    
    echo ""
done

echo "All tests completed!"
