#!/bin/bash

# Compile the parallel implementation
g++-14 -O3 -fopenmp -std=c++17 -I/opt/homebrew/include -L/opt/homebrew/lib -o pbwt_haps pbwt-haps.cpp util.cpp -lboost_iostreams -lz

# Define test files
test_files=(
    "10SitesX10Haps.hap"
    "20SitesX10Haps.hap"
    "30SitesX20Haps.hap"
    "33SitesX22Haps.hap"
    "100SitesX50Haps.hap"
)

# Run tests for each file
for file in "${test_files[@]}"; do
    echo "===== Testing $file ====="
    
    # Run the parallel implementation
    ./pbwt_haps "$file" > "parallel_output_$file.txt"
    
    # Extract just the maximal matches section
    grep -A1000 "========== Maximal Matches ==========" "parallel_output_$file.txt" > "maximal_$file.txt"
    
    # Count the number of matches
    match_count=$(grep -c "Position" "maximal_$file.txt")
    
    # Check if there are any matches
    if [ $match_count -gt 0 ]; then
        echo "SUCCESS: Found $match_count maximal matches for $file"
        
        # Show the first 5 matches
        echo "First 5 matches:"
        grep "Position" "maximal_$file.txt" | head -n 5
        
        # Show the last 5 matches
        echo "Last 5 matches:"
        grep "Position" "maximal_$file.txt" | tail -n 5
    else
        echo "ERROR: No maximal matches found for $file"
    fi
    
    echo ""
done

echo "All tests completed!"
