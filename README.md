# NanoExpansion: a tool for the characterization of Repeat Expansion Pattern in Nanopore sequencing samples

NanoExpansion is a python software for the extraction and characterization of Short Tandem Repeats (STRs) data from nanopore sequencing.
It exploits the result from straglr to generate plots of the expansion site of the region of interest (e.g. gene DMPK for DM1) and to return the compact expansion pattern string.

## Requirements

Some files are needed in order to run NanoExpansion:

* a sorted and indexed .bam file of the sample of interest
* .tsv and .vcf output files from straglr
* the catalog for STR annotation with Stranger
* a .bed file with the region and the motif of expansion

Moreover, the folder structure must be the following

```
sample/
│
├── nanoexpansion/
    ├── sample-straglr.tsv
    ├── sample-straglr.vcf    
    ├── sample.sort.bam    
    ├── sample.sort.bam.bai    
    ├── variant_catalog_hg38.json    
    ├── bed_filter.bed    
    └── wf_str_repeats.bed
```

and the required files must be inside nanoexpansion folder.

Depending on the straglr version used, you would need to transform the output .tsv file in order to have only the following columns:

'chrom', 'start', 'end', 'repeat_unit', 'genotype', 'read', 'copy_number', 'size', 'read_start', 'strand', 'allele'

If your .tsv file does not satisfy this requirement, you should first run

```bash
python transform_straglr_tsv.py --input sample-straglr_old.tsv --output sample-straglr.tsv
```

If your version of straglr does not output the .vcf file, you can create it starting from the .tsv and the .bed files, by running:

```bash
python create_vcf_file.py --tsv sample-straglr.tsv --bed sample-straglr.bed --vcf sample-straglr.vcf
```

## How to use NanoExpansion

1. Download the repository

```bash
git clone https://github.com/Cesco16/NanoExpansion.git
cd NanoExpansion
```

2. Create and activate the conda environment 

```bash
conda env create -f requirements.yaml
```

```bash
conda activate nanoexpansion
```

3. Index .bam STR file and keep only reads with STR of interest

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
samtools view --write-index -N reads_to_filter.txt -o native9411_str_reads.bam native9411_str_regions.bam
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
```bash
NanoExpansion.py --sample 9411native --repeat CAG --interruption CAA --path /path/to/file/GridIon
```

N.B. Please, do not change the filenames created in steps 3-5.

## Example of usage

Here an example of NanoExpansion applied to a patient affected by Mytonic Dystrophy type 1 (DM1), which is characterized by an expansion of the CTG triplet in gene *DMPK*.
Thanks to NanoExpansion, it is possible to characterize the wild-type allele.
The numbers in the plots represents the number of nucleotides in each region. The number of repeats is obtained dividing those numbers by the length of the repeat motif (in this case, 3).

![Example of wild-type allele in gene DMPK](images/str_eg2.PNG)

and also the mutated reads. Here an example of an expanded read, that shows a TTG interruption pattern:

![Example of STR with interruption pattern in gene DMPK](images/STR_eg.PNG)

Finally, NanoExpansion returns the complete characterization of repeat patterns in all the available reads:

```
1145b1e2-58bb-433c-afa1-939a27d713f3 :  (CTG)8
fddbf9a9-73b1-4bfb-bab8-4e386dad1720 :  (CTG)12
15234493-05c1-45e8-a844-6a8f88846125 :  (CTG)12
9ae4ba29-c4f2-493e-b67a-74254b9bd9a5 :  (CTG)11
551d2b3a-7a47-4dc9-bd14-d6e227cffab3 :  (CTG)648(TTG)1(CTG)132
8b3c7dcb-3e01-4642-b1c0-7fa506faf26c :  (CTG)114(TTG)757(CTG)91
17c4a40a-8861-4141-99aa-f5a9440e5166 :  (CTG)12
7cae6e57-bd34-4410-aac5-1bb2024430be :  (CTG)87(TTG)317(GTC)(TTG)145(TGTCG)(TTG)21(TCG)(TTG)234(TC)(TTG)237(CTG)38
503ccf68-7c1c-4350-a7cf-83d1ae02b101 :  (CTG)12
8ca91fb6-e90c-4fde-828d-4d9df868ae6a :  (CTG)12
fa9331e7-441a-4d13-bace-b6b2c5e11a40 :  (CTG)69(CCGCCG)(CTG)35(CCGCCGCCG)(CTG)22(TTG)118(TCG)(TTG)31(TCG)(TTG)123(CTG)118(CCGCCCG)(CTG)198
```

## License

This project is licensed under the [MIT License](LICENSE).  
You are free to use, modify, and distribute this software under the terms of the license.

## Citation

If you use NanoExpansion in your research or work, please cite the GitHub repository:

```
@misc{NanoExpansion
author = {Francesco Casadei},
title = {NanoExpansion: a tool for the characterization of Repeat Expansion Pattern in Nanopore sequencing samples},
year = {2024},
publisher = {GitHub},
journal = {GitHub repository},
howpublished = {\url{https://github.com/Cesco16/NanoExpansion)}
}
```
