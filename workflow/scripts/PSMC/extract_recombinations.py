## This is a script I asked chat gpt to write to extract the n_sum information from the header of each of
# the .psmc files in order to compare total recombinations for QC 

import os

def parse_file(file_path):
    # Initialize variables to store the extracted values
    n_seqs = None
    sum_L = None
    sum_n = None

    # Open the file and extract the required values
    with open(file_path, 'r') as file:
        for line in file:
            if "n_seqs:" in line:
                n_seqs = int(line.split("n_seqs:")[1].split(",")[0].strip())
            if "sum_L:" in line:
                sum_L = int(line.split("sum_L:")[1].split(",")[0].strip())
            if "sum_n:" in line:
                sum_n = int(line.split("sum_n:")[1].strip())
            # Stop reading further once all values are found
            if n_seqs is not None and sum_L is not None and sum_n is not None:
                break

    return n_seqs, sum_L, sum_n


def process_files(input_dir, output_tsv):
    # List all files in the input directory
    files = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]

    # Open the output TSV file for writing
    with open(output_tsv, 'w') as out_file:
        # Write the header
        out_file.write("Filename\tN_Seqs\tSum_L\tSum_N\n")

        # Process each file
        for file in files:
            file_path = os.path.join(input_dir, file)
            base_name = os.path.splitext(file)[0]  # Get the base filename
            n_seqs, sum_L, sum_n = parse_file(file_path)
            # Write the extracted values to the TSV file
            out_file.write(f"{base_name}\t{n_seqs}\t{sum_L}\t{sum_n}\n")


# Example usage:
input_directory = "path/to/your/input/directory"
output_file = "output.tsv"
process_files(input_directory, output_file)
