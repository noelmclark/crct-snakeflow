#!/bin/bash

# Usage: ./calculate_coverage_percentage.sh input.bam output.txt

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.bam> <output.txt>"
    exit 1
fi

# Input arguments
input_bam=$1       # First argument: input BAM file
output_file=$2     # Second argument: output file for percent coverage

# Calculate bpCountZero (number of base pairs with 0 coverage)
zero=$(bedtools genomecov -ibam "$input_bam" -bga | \
       awk '$4==0 {bpCountZero+=($3-$2)} END {print bpCountZero}')

# Calculate bpCountNonZero (number of base pairs with non-zero coverage)
nonzero=$(bedtools genomecov -ibam "$input_bam" -bga | \
          awk '$4>0 {bpCountNonZero+=($3-$2)} END {print bpCountNonZero}')

# Check if both counts are available to avoid division by zero
if [ -z "$zero" ] || [ -z "$nonzero" ]; then
    echo "Error: Unable to calculate coverage."
    exit 1
fi

# Calculate the percentage of non-zero coverage
percent=$(bc <<< "scale=6; ($nonzero / ($zero + $nonzero)) * 100")

# Write the percentage to the output file
echo "$percent" > "$output_file"

# Output a success message
echo "Coverage percentage calculated and written to $output_file"
