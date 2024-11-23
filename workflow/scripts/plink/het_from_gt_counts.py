# thank you chat gpt
# run with python het_from_gt_counts.py <input_directory> <output_file>
import os
import glob
import sys

def process_files(input_dir, output_file):
    # Initialize a list to hold the parsed data
    data = []

    # Iterate over all input files in the directory
    for file_path in glob.glob(os.path.join(input_dir, "*.txt")):
        # Extract the sample name from the file name
        full_name = os.path.splitext(os.path.basename(file_path))[0]
        sample_name = full_name.split("-")[0]  # Take only the first part before the first '-'

        # Initialize values
        sum_het_ref_alt = total_variants = average_het_ref_alt = None

        # Read the file line by line
        with open(file_path, "r") as f:
            for line in f:
                if line.startswith("Sum of HET_REF_ALT_CTS:"):
                    sum_het_ref_alt = line.split(":")[1].strip()
                elif line.startswith("Total variants:"):
                    total_variants = line.split(":")[1].strip()
                elif line.startswith("Average HET_REF_ALT_CTS:"):
                    average_het_ref_alt = line.split(":")[1].strip()

        # Append the extracted values to the data list
        if sum_het_ref_alt is not None and total_variants is not None and average_het_ref_alt is not None:
            data.append([sample_name, sum_het_ref_alt, total_variants, average_het_ref_alt])

    # Write the output to a TSV file
    with open(output_file, "w") as out:
        # Write the header with the updated column names
        out.write("sample\tsum_het\ttotal_var\tavg_het\n")

        # Write the data rows
        for row in data:
            out.write("\t".join(row) + "\n")

    print(f"Summary file created: {output_file}")

if __name__ == "__main__":
    # Check for correct number of arguments
    if len(sys.argv) != 3:
        print("Usage: python het_from_gt_counts.py <input_directory> <output_file>")
        sys.exit(1)

    # Get the input directory and output file from command-line arguments
    input_dir = sys.argv[1]
    output_file = sys.argv[2]

    # Run the processing function
    process_files(input_dir, output_file)
