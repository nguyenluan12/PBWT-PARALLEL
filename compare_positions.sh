#!/bin/bash

# Compile both implementations
g++-14 -O3 -std=c++17 -I/opt/homebrew/include -L/opt/homebrew/lib -o pbwt_single pbwt-single.cpp util.cpp -lboost_iostreams -lz
g++-14 -O3 -fopenmp -std=c++17 -I/opt/homebrew/include -L/opt/homebrew/lib -o pbwt_haps pbwt-haps.cpp util.cpp -lboost_iostreams -lz

# Define test files
test_files=(
    "33SitesX22Haps.hap"
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
    
    # Count the number of matches for each position in single output
    echo "Position counts in single output:"
    grep "Position" "single_maximal_$file.txt" | cut -d' ' -f2 | cut -d':' -f1 | sort | uniq -c
    
    # Count the number of matches for each position in parallel output
    echo "Position counts in parallel output:"
    grep "Position" "parallel_maximal_$file.txt" | cut -d' ' -f2 | cut -d':' -f1 | sort | uniq -c
    
    # Compare specific positions
    echo "Comparing position 0 matches:"
    echo "Single output position 0:"
    grep "Position 0:" "single_maximal_$file.txt"
    echo "Parallel output position 0:"
    grep "Position 0:" "parallel_maximal_$file.txt"
    
    echo ""
done

echo "All tests completed!"
