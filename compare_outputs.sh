#!/bin/bash

# Compile both implementations
g++-14 -O3 -fopenmp -std=c++17 -I/opt/homebrew/include -L/opt/homebrew/lib -o pbwt_haps pbwt-haps.cpp util.cpp -lboost_iostreams -lz
g++-14 -O3 -std=c++17 -I/opt/homebrew/include -L/opt/homebrew/lib -o pbwt_single pbwt-single.cpp util.cpp -lboost_iostreams -lz

# Test files
test_files=("10SitesX10Haps.hap" "20SitesX10Haps.hap" "30SitesX20Haps.hap" "33SitesX22Haps.hap" "100SitesX50Haps.hap")

# Run tests and compare outputs
for file in "${test_files[@]}"; do
    echo "Testing $file..."
    
    # Run both implementations and save outputs
    ./pbwt_haps "$file" > haps_output.txt
    ./pbwt_single "$file" > single_output.txt
    
    # Extract just the maximal matches section from both outputs
    grep -A 1000 "========== Maximal Matches ==========" haps_output.txt > haps_maximal.txt
    grep -A 1000 "========== Maximal Matches ==========" single_output.txt > single_maximal.txt
    
    # Compare the maximal matches sections
    diff_result=$(diff haps_maximal.txt single_maximal.txt)
    
    if [ -z "$diff_result" ]; then
        echo "✅ MATCH: Outputs for $file match exactly!"
    else
        echo "❌ MISMATCH: Outputs for $file differ!"
        echo "Differences:"
        echo "$diff_result"
    fi
    
    echo "-----------------------------------------"
done

# Clean up temporary files
rm haps_output.txt single_output.txt haps_maximal.txt single_maximal.txt
