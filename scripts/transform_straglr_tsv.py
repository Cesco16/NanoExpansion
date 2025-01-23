import csv
import argparse

parser = argparse.ArgumentParser(description='missing_data',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--input')
parser.add_argument('--output')
args = parser.parse_args()

input_file = args.input
output_file = args.output


# Define the columns to keep and their new names
columns_to_keep = {
    'chrom': 'chrom',
    'start': 'start',
    'end': 'end',
    'target_repeat': 'repeat_unit',
    'genotype': 'genotype',
    'read_name': 'read',
    'copy_number': 'copy_number',
    'size': 'size',
    'read_start': 'read_start',
    'strand': 'strand',
    'allele': 'allele'
}

# Read the input file and write the transformed data to the output file
with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.DictReader(infile, delimiter='\t')
    writer = csv.DictWriter(outfile, fieldnames=columns_to_keep.values(), delimiter='\t')

    # Write the header
    outfile.write(next(infile))  # Write the first line as is
    writer.writeheader()

    # Write the transformed rows
    for row in reader:
        transformed_row = {columns_to_keep[key]: row[key] for key in columns_to_keep}
        writer.writerow(transformed_row)


#print(f"Transformed TSV file saved as: {output_tsv}")
