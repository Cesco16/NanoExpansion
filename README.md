# NanoExpansion

**NanoExpansion** is a Python software package for the comprehensive characterization of pathogenic repeat expansions from long-read sequencing data.

Unlike conventional repeat genotyping tools, NanoExpansion reconstructs the internal architecture of expanded alleles, identifies interruption motifs, extracts locus-specific DNA methylation information, and produces publication-ready graphical representations of repeat structures.

---

## Features

NanoExpansion provides:

* reconstruction of repeat architectures from long-read sequencing data;
* automatic extraction of reads spanning repeat expansion loci;
* identification of canonical repeat motifs;
* identification of interruption motifs (user-defined or automatically detected);
* compact representation of repeat architectures (e.g. `(CTG)120(CCG)3(CTG)42`);
* allele clustering;
* repeat length statistics;
* locus-specific DNA methylation analysis;
* graphical visualization of repeat organization;
* publication-ready figures and summary tables.

---

## Workflow

```
Aligned BAM
      │
      ▼
Read extraction
      │
      ▼
Repeat parsing
      │
      ▼
Motif identification
      │
      ▼
Architecture reconstruction
      │
      ├────────► Methylation extraction
      │
      ▼
Allele clustering
      │
      ▼
Plots + Reports
```

---

## Installation

Clone the repository

```bash
git clone https://github.com/<username>/NanoExpansion.git

cd NanoExpansion
```

Create a conda environment

```bash
conda create -n nanoexpansion python=3.11

conda activate nanoexpansion
```

Install dependencies

```bash
pip install -r requirements.txt
```

or

```bash
pip install .
```

---

## Requirements

Python ≥ 3.10

Main dependencies

* pysam
* numpy
* pandas
* scipy
* matplotlib
* scikit-learn
* edlib
* tqdm

---

## Input

NanoExpansion requires:

* aligned BAM file
* indexed BAM (.bai)
* reference genome
* genomic coordinates of the repeat locus
* canonical repeat motif

Optionally:

* interruption motifs
* methylation tags
* sample name

---

## Usage

Basic analysis

```bash
nanoexpansion \
    --bam sample.bam \
    --reference hg38.fa \
    --gene DMPK
```

Specify interruption motifs

```bash
nanoexpansion \
    --bam sample.bam \
    --reference hg38.fa \
    --gene DMPK \
    --interruptions CCG CAG
```

Automatic interruption discovery

```bash
nanoexpansion \
    --bam sample.bam \
    --reference hg38.fa \
    --gene DMPK \
    --discover-interruptions
```

---

## Command line arguments

| Argument                 | Description                               |
| ------------------------ | ----------------------------------------- |
| --bam                    | aligned BAM file                          |
| --reference              | reference genome                          |
| --gene                   | target locus                              |
| --sample                 | sample name                               |
| --output                 | output directory                          |
| --motif                  | canonical repeat motif                    |
| --interruptions          | interruption motifs                       |
| --discover-interruptions | automatically search for secondary motifs |
| --min-support            | minimum number of supporting reads        |
| --threads                | number of CPU threads                     |
| --plots                  | generate figures                          |
| --methylation            | extract methylation information           |

---

## Output

NanoExpansion generates

```
results/

summary.tsv

alleles.tsv

repeat_architecture.tsv

interruptions.tsv

methylation.tsv

plots/

report.html
```

---

## Repeat architecture

Example

```
(CTG)138

(CTG)94(CCG)3(CTG)41

(GAA)812(GGA)2(GAA)174
```

---

## Methylation analysis

When methylation tags are available in the BAM file (e.g. Oxford Nanopore), NanoExpansion extracts

* CpG methylation

* allele-specific methylation

* locus average methylation

* methylation plots

---

## Automatic interruption discovery

If interruption motifs are not supplied, NanoExpansion can automatically search for recurrent secondary motifs within the repeat tract.

Detected motifs are ranked according to

* frequency
* number of supporting reads
* repeat consistency

allowing the identification of candidate interruption motifs without prior knowledge.

---

## Supported repeat expansion disorders

NanoExpansion is locus-independent and can be applied to any repeat expansion disorder by specifying the genomic coordinates and repeat motif.

Examples include

* Myotonic dystrophy type 1 (*DMPK*)
* Huntington disease (*HTT*)
* Fragile X syndrome (*FMR1*)
* CANVAS (*RFC1*)
* Friedreich ataxia (*FXN*)
* Spinocerebellar ataxias
* C9orf72-related ALS/FTD

---

## Citation

If you use NanoExpansion in your work, please cite

> *Citation will be added after publication.*

---

## License

MIT License

---

## Acknowledgements

NanoExpansion was developed at the Department of Medical and Surgical Sciences, University of Bologna.
