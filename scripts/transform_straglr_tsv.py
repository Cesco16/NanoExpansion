import csv
import argparse

parser = argparse.ArgumentParser(description='missing_data',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input')
parser.add_argument('--output')
args = parser.parse_args()

input_tsv = args.input
output_tsv = args.output


def transform_tsv(input_file, output_file):
    # Define the columns to keep and the column renaming mapping
    columns_to_keep = [
        "#chrom", "start", "end", "target_repeat", "genotype", "read_name",
        "copy_number", "size", "read_start", "strand", "allele"
    ]
    renamed_columns = {
        "target_repeat": "repeat_unit",
        "read_name": "read"
    }

    # Open the input and output TSV files
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        # Read the input file
        reader = csv.DictReader(infile, delimiter='\t')
        # Determine the new column order
        new_columns = [renamed_columns.get(col, col) for col in columns_to_keep]
        # Create the writer with the new column headers
        writer = csv.DictWriter(outfile, fieldnames=new_columns, delimiter='\t')
        writer.writeheader()

        # Process each row
        for row in reader:
            # Keep only the necessary columns and rename them
            new_row = {renamed_columns.get(col, col): row[col] for col in columns_to_keep if col in row}
            writer.writerow(new_row)

# Transform the TSV file
transform_tsv(input_tsv, output_tsv)

print(f"Transformed TSV file saved as: {output_tsv}")
