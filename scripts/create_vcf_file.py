import csv
import os
import argparse

parser = argparse.ArgumentParser(description='missing_data',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--tsv')
parser.add_argument('--bed')
parser.add_argument('--vcf')
parser.add_argument('--gene')
args = parser.parse_args()

input_tsv = args.tsv
input_bed = args.bed
output_vcf = args.vcf
gene = args.gene

def create_vcf_with_custom_entry(tsv_path, bed_path, vcf_path):
    # Check if input files exist
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"TSV file not found: {tsv_path}")
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"BED file not found: {bed_path}")

    # Read BED data into a dictionary for quick lookup
    bed_data = {}
    with open(bed_path, 'r') as bed_file:
        bed_reader = csv.DictReader(bed_file, delimiter='\t')
        for row in bed_reader:
            key = (row['#chrom'], row['start'], row['end'])
            bed_data[key] = row

    # Open input TSV and output VCF file
    with open(tsv_path, 'r') as tsv_file, open(vcf_path, 'w') as vcf_file:
        # Skip metadata header lines
        lines = tsv_file.readlines()
        metadata_lines = [line for line in lines if line.startswith('#') and not line.startswith('#chrom')]
        tsv_data_lines = [line for line in lines if line.startswith('#chrom') or not line.startswith('#')]

        # Parse TSV data
        tsv_reader = csv.DictReader(tsv_data_lines, delimiter='\t')

        # Write VCF header
        vcf_file.write(
            "##fileformat=VCFv4.2\n"
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n"
            "##INFO=<ID=REF,Number=1,Type=Integer,Description=\"Reference copy number\">\n"
            "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"Repeat length\">\n"
            "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat unit\">\n"
            "##INFO=<ID=REPID,Number=1,Type=String,Description=\"Repeat locus ID\">\n"
            "##INFO=<ID=VARID,Number=1,Type=String,Description=\"Variant ID\">\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "##FORMAT=<ID=SO,Number=1,Type=String,Description=\"Support observation type\">\n"
            "##FORMAT=<ID=CN,Number=2,Type=String,Description=\"Copy numbers for alleles\">\n"
            "##FORMAT=<ID=CI,Number=2,Type=String,Description=\"Confidence intervals\">\n"
            "##FORMAT=<ID=AD_SP,Number=2,Type=String,Description=\"Allele depth from spanning reads\">\n"
            "##FORMAT=<ID=AD_FL,Number=2,Type=String,Description=\"Allele depth from flanking reads\">\n"
            "##FORMAT=<ID=AD_IR,Number=2,Type=String,Description=\"Allele depth from in-repeat reads\">\n"
        )
        for metadata in metadata_lines:
            vcf_file.write(f"##{metadata.strip()}\n")
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")

        # Process TSV rows and create a custom VCF entry
        for row in tsv_reader:
            chrom = row['#chrom']
            start = int(row['start']) #+ 1  # Adjust POS to be 1-based
            end = row['end']
            repeat_unit = row['target_repeat']#['repeat_unit']#['target_repeat']
            ref_copy_number = 20  # Example reference copy number
            repeat_length = 60  # Example repeat length
            rep_id = str(gene) #"DMPK"
            var_id = str(gene) #"DMPK"
            allele1_cn = 1060  # Example allele copy number
            allele2_cn = 13
            allele1_ci = "1014-1107"
            allele2_ci = "12-13"
            ad_sp = "2/7"  # Example spanning reads
            ad_fl = "0/0"  # Example flanking reads
            ad_ir = "0/0"  # Example in-repeat reads

            # INFO and FORMAT fields
            info = (
                f"SVTYPE=STR;END={end};REF={ref_copy_number};RL={repeat_length};"
                f"RU={repeat_unit};REPID={rep_id};VARID={var_id}"
            )
            format_field = "GT:SO:CN:CI:AD_SP:AD_FL:AD_IR"
            sample_data = f"1/2:SPANNING/SPANNING:{allele1_cn}/{allele2_cn}:{allele1_ci}/{allele2_ci}:{ad_sp}:{ad_fl}:{ad_ir}"

            # Write VCF row
            vcf_file.write(
                f"{chrom}\t{start}\t.\tC\t<STR1060>,<STR13>\t.\tPASS\t{info}\t{format_field}\t{sample_data}\n"
            )

def create_vcf_with_dynamic_info(tsv_path, bed_path, vcf_path):
    # Check if input files exist
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"TSV file not found: {tsv_path}")
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"BED file not found: {bed_path}")

    # Read BED data into a dictionary for quick lookup
    bed_data = {}
    with open(bed_path, 'r') as bed_file:
        bed_reader = csv.DictReader(bed_file, delimiter='\t')
        for row in bed_reader:
            key = (row['#chrom'], row['start'], row['end'])
            bed_data[key] = row

    # Open input TSV and output VCF file
    with open(tsv_path, 'r') as tsv_file, open(vcf_path, 'w') as vcf_file:
        # Skip metadata header lines
        lines = tsv_file.readlines()
        metadata_lines = [line for line in lines if line.startswith('#') and not line.startswith('#chrom')]
        tsv_data_lines = [line for line in lines if line.startswith('#chrom') or not line.startswith('#')]

        # Parse TSV data
        tsv_reader = csv.DictReader(tsv_data_lines, delimiter='\t')

        # Write VCF header
        vcf_file.write(
            "##fileformat=VCFv4.2\n"
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n"
            "##INFO=<ID=REF,Number=1,Type=Integer,Description=\"Reference copy number\">\n"
            "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"Repeat length\">\n"
            "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat unit\">\n"
            "##INFO=<ID=REPID,Number=1,Type=String,Description=\"Repeat locus ID\">\n"
            "##INFO=<ID=VARID,Number=1,Type=String,Description=\"Variant ID\">\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "##FORMAT=<ID=SO,Number=1,Type=String,Description=\"Support observation type\">\n"
            "##FORMAT=<ID=CN,Number=2,Type=String,Description=\"Copy numbers for alleles\">\n"
            "##FORMAT=<ID=CI,Number=2,Type=String,Description=\"Confidence intervals\">\n"
            "##FORMAT=<ID=AD_SP,Number=2,Type=String,Description=\"Allele depth from spanning reads\">\n"
            "##FORMAT=<ID=AD_FL,Number=2,Type=String,Description=\"Allele depth from flanking reads\">\n"
            "##FORMAT=<ID=AD_IR,Number=2,Type=String,Description=\"Allele depth from in-repeat reads\">\n"
        )
        for metadata in metadata_lines:
            vcf_file.write(f"##{metadata.strip()}\n")
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")

        # Process TSV rows and create a VCF entry
        for row in tsv_reader:
            chrom = row['#chrom']
            start = int(row['start']) #+ 1  # Adjust POS to be 1-based
            end = row['end']
            repeat_unit = row['target_repeat']
            ref_copy_number = row.get('copy_number', 20)  # Use copy number or default value
            repeat_length = len(repeat_unit) * int(ref_copy_number)  # Calculate repeat length
            rep_id = row['locus']
            var_id = row['locus']
            
            # Lookup BED entry
            bed_key = (chrom, row['start'], row['end'])
            if bed_key in bed_data:
                bed_row = bed_data[bed_key]
                allele1_cn = bed_row['allele1:copy_number']
                allele2_cn = bed_row['allele2:copy_number']
                allele1_ci = "1014-1107"  # Placeholder; adjust based on actual data
                allele2_ci = "12-13"  # Placeholder
                ad_sp = "2/7"  # Placeholder; use actual data
                ad_fl = "0/0"  # Placeholder
                ad_ir = "0/0"  # Placeholder
            else:
                raise ValueError(f"No matching BED entry found for {bed_key}")

            # Construct ALT
            alt = f"<STR{allele1_cn}>,<STR{allele2_cn}>"

            # INFO and FORMAT fields
            info = (
                f"SVTYPE=STR;END={end};REF={ref_copy_number};RL={repeat_length};"
                f"RU={repeat_unit};REPID={rep_id};VARID={var_id}"
            )
            format_field = "GT:SO:CN:CI:AD_SP:AD_FL:AD_IR"
            sample_data = f"1/2:SPANNING/SPANNING:{allele1_cn}/{allele2_cn}:{allele1_ci}/{allele2_ci}:{ad_sp}:{ad_fl}:{ad_ir}"

            # Write VCF row
            vcf_file.write(
                f"{chrom}\t{start}\t.\tC\t{alt}\t.\tPASS\t{info}\t{format_field}\t{sample_data}\n"
            )

create_vcf_with_custom_entry(input_tsv, input_bed, output_vcf)

print(f"VCF file created: {output_vcf}")
