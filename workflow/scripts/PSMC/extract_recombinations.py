## This is a script I asked chat gpt to write to extract the n_sum information from the header of each of
# the .psmc files in order to compare total recombinations for QC 

# run the script with
# python extract_recombinations.py <input_directory> <output_file>

import os
import re
import argparse

def extract_values_from_file(file_path):
    """
    Extract values from a single file: n_seqs, sum_L, and sum_n.
    """
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Extract n_seqs, sum_L, and sum_n values using regular expressions
    n_seqs = re.search(r'n_seqs:(\d+)', content)
    sum_L = re.search(r'sum_L:(\d+)', content)
    sum_n = re.search(r'sum_n:(\d+)', content)
    
    # Convert extracted values to integers or default to 0 if not found
    n_seqs = int(n_seqs.group(1)) if n_seqs else 0
    sum_L = int(sum_L.group(1)) if sum_L else 0
    sum_n = int(sum_n.group(1)) if sum_n else 0
    
    return n_seqs, sum_L, sum_n

def process_files(input_directory, output_file):
    """
    Process all files in the input directory and write results to a TSV file.
    """
    # Open the output file for writing
    with open(output_file, 'w') as out_f:
        # Write the header
        out_f.write("Base_Name\tn_seqs\tsum_L\tsum_n\n")
        
        # Iterate over all files in the input directory
        for filename in os.listdir(input_directory):
            file_path = os.path.join(input_directory, filename)
            
            # Skip if it's not a file
            if not os.path.isfile(file_path):
                continue
            
            # Extract values from the current file
            n_seqs, sum_L, sum_n = extract_values_from_file(file_path)
            
            # Write the results to the output file
            base_name = os.path.splitext(filename)[0]
            out_f.write(f"{base_name}\t{n_seqs}\t{sum_L}\t{sum_n}\n")

def main():
    """
    Main function to handle command-line arguments and execute the script.
    """
    parser = argparse.ArgumentParser(description="Extract and summarize data from files.")
    parser.add_argument("input_directory", help="Path to the input directory containing files.")
    parser.add_argument("output_file", help="Path to the output TSV file.")
    
    args = parser.parse_args()
    
    # Process the files using provided arguments
    process_files(args.input_directory, args.output_file)

if __name__ == "__main__":
    main()
