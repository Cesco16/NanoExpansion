# NanoExpansion

Tool for extraction and characterization of Short Tandem Repeats (STRs) data from nanopore sequencing.

## How to use NanoExpansion

1. Download the repository

2. Create and activate the conda environment 

```bash
conda env create -f requirements.yaml

```

```bash
conda activate nanoexpansion
```

3. Index .bam STR file

```bash
samtools view -b -h -o native9411_str_regions.bam -L ../native13204/bed_filter.bed native9411_sort.bam
```
```bash
samtools index native9411_str_regions.bam
```
```bash
tail -n +3 native9411_straglr.tsv | cut -f 6 > reads_to_filter.txt
```
```bash
samtools view --write-index -N reads_to_filter.txt -o native9411_str_reads.bam
```

4. Annotate .vcf using Stranger

```bash
stranger -f ../variant_catalog_hg38.json native9411_straglr.vcf > native9411_straglr_annot.vcf | sed 's/\\ / _/g'
```
```bash
bgzip native9411_straglr_annot.vcf
```
```bash
tabix native9411_straglr_annot.vcf.gz
```

5. Extract fields of interest
```bash
SnpSift extractFields native9411_straglr_annot.vcf.gz CHROM POS ALT FILTER REF RL RU REPID VARID STR_STATUS > native9411_rep_annot.tsv
```
```bash
SnpSift extractFields native9411_straglr_annot.vcf.gz CHROM POS DisplayRU STR_NORMAL_MAX STR_PATHOLOGIC_MIN VARID Disease > native9411_rep_plot.tsv
```
6. Execute NanoExpansion

NanoExpansion.py --sample 9411native --repeat CAG --interruption CAA --path /path/to/file/GridIon

## Example of usage

![Example of STR with interruption pattern in gene DMPK](images/STR_eg.PNG)

![Example of wild-type allele in gene DMPK](images/str_eg2.PNG)
