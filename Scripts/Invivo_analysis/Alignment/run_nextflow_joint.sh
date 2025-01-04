#!/bin/bash

# this contains featureCounts and SNPsplit
source ~/miniconda3/bin/activate sequencing_tools

module load fastqc
module load trim-galore/0.6.6
module load STAR
module load samtools/1.9
module load picard
module load R
module load deeptools/3.5.1

CONFIG_FILE=/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/nextflow.config.exp2
CONFIG_FILE=nextflow.config.exp2

nextflow run ../asRNA_pipelines/allele-specific_RNA-seq/allelic_RNA-seq.nf -c $CONFIG_FILE -resume
