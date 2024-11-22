# thank you chat gpt
# run with: python process_variants.py <input_file>
import pandas as pd
import sys

def main():
    # Check for command-line arguments
    if len(sys.argv) != 2:
        print("Usage: python process_variants.py <input_file>")
        sys.exit(1)

    # Get the input file path from the command-line argument
    input_file = sys.argv[1]

    try:
        # Read the file into a DataFrame
        df = pd.read_csv(input_file, sep="\t", comment='#')

        # Filter rows where MISSING_CT is greater than 0
        filtered_df = df[df['MISSING_CT'] == 0]

        # Sum the HET_REF_ALT_CTS column
        het_ref_alt_sum = filtered_df['HET_REF_ALT_CTS'].sum()

        # Count the total number of remaining rows
        total_variants = len(filtered_df)

        # Calculate the average
        average_het_ref_alt = het_ref_alt_sum / total_variants if total_variants > 0 else 0

        # Print the results
        print(f"Sum of HET_REF_ALT_CTS: {het_ref_alt_sum}")
        print(f"Total variants: {total_variants}")
        print(f"Average HET_REF_ALT_CTS: {average_het_ref_alt}")
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
