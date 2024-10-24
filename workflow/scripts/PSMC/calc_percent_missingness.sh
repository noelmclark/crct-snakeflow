#!/bin/bash

# Usage: ./calculate_coverage_percentage.sh output_file input1.bam input2.bam ... inputN.bam

# Check if at least two arguments are provided (one for the output file, and at least one BAM file)
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <output_file> <input1.bam> [input2.bam ... inputN.bam]"
    exit 1
fi

# Output file (first argument)
output_file=$1

# Shift the first argument (output_file) so that "$@" contains only BAM files
shift

# Initialize the output file (optional, this clears any existing content)
> "$output_file"

# Loop over each input BAM file
for input_bam in "$@"; do

    # Get the base name of the input BAM file (without extension)
    base_name=$(basename "$input_bam" .bam)

    # Calculate bpCountZero (number of base pairs with 0 coverage)
    zero=$(bedtools genomecov -ibam "$input_bam" -bga | \
           awk '$4==0 {bpCountZero+=($3-$2)} END {print bpCountZero}')

    # Calculate bpCountNonZero (number of base pairs with non-zero coverage)
    nonzero=$(bedtools genomecov -ibam "$input_bam" -bga | \
              awk '$4>0 {bpCountNonZero+=($3-$2)} END {print bpCountNonZero}')

    # Check if both counts are available to avoid division by zero
    if [ -z "$zero" ] || [ -z "$nonzero" ]; then
        echo "Error: Unable to calculate coverage for $input_bam" >> "$output_file"
        continue
    fi

    # Calculate the percentage of non-zero coverage
    percent=$(bc <<< "scale=6; ($nonzero / ($zero + $nonzero)) * 100")

    # Append the result to the output file with the BAM file name
    echo "$base_name: $percent%" >> "$output_file"

    # Output a success message for each file
    echo "Coverage percentage calculated for $input_bam"

done

# Final success message
echo "All results written to $output_file"
