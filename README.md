# Salmon pseudo-alignment WTS pipeline

Snakemake workflow for paired-end whole transcriptome sequencing (WTS) preprocessing and gene expression quantification using **fastp** and **Salmon**.

This pipeline performs:
- paired-end FASTQ quality trimming with **fastp**
- transcript-level quantification with **Salmon**
- gene-level quantification from Salmon output using an **R script**

The current workflow is designed for **stranded paired-end RNA-seq** and runs Salmon with library type `ISR`.

---

## Overview

The pipeline processes paired-end FASTQ files listed in a sample sheet and generates trimmed reads, QC reports, transcript-level quantification, gene-level quantification, and sample metadata summaries.

### Main steps

1. **Read sample sheet**
   - loads sample names and FASTQ locations from a CSV file

2. **Quality trimming with fastp**
   - adapter detection for paired-end reads
   - trimming of low-quality bases from read ends
   - minimum read length filtering
   - HTML and JSON QC reports

3. **Pseudo-alignment and transcript quantification with Salmon**
   - quantifies transcript abundance against a prebuilt Salmon transcriptome index
   - uses `--validateMappings`, `--seqBias`, and `--gcBias`

4. **Gene-level quantification**
   - converts transcript-level Salmon output into gene-level expression using an R helper script

5. **Sample metadata output**
   - writes a simple text file per sample with the original FASTQ paths and folder information

---

## Workflow outputs

For each sample, the pipeline generates:

- trimmed FASTQ files  
  `trimmed/{sample}_R1.trimmed.fastq.gz`  
  `trimmed/{sample}_R2.trimmed.fastq.gz`

- fastp QC reports  
  `qc/fastp/{sample}.html`  
  `qc/fastp/{sample}.json`

- Salmon transcript quantification  
  `quants/{sample}_ISR/quant.sf`

- gene-level quantification  
  `quants/{sample}_ISR/quant.gene.sf`

- sample metadata  
  `quants/{sample}_ISR/{sample}.txt`

---

## Repository structure

A recommended structure for this repository is:

```text
salmon-pseudo-alignment-wts-pipeline/
├── Snakefile
├── README.md
├── config/
│   └── config.example.yaml
├── samples/
│   └── samples.example.csv
├── helper_scripts/
│   └── gene_quant_generation.R
├── envs/
│   └── salmon.yaml
└── .gitignore

##Notes on libry type:
This workflow currently runs Salmon with **-l ISR** (read 1 maps to the reverse strand, read 2 maps to the forward strand).
Please verify that this matches your library preparation protocol before running the workflow on new datasets


##Runing the pipeline: 
snakemake --cores 20 --use-conda --configfile config/config.yaml

