#!/bin/bash

module load bowtie2/2.3.5.1

sample_id=$1

DATA_DIR=/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/data/20231213_invivo_exp4/2023-12-11-AAFCTG2M5/

fastq_1=$DATA_DIR/${sample_id}_1_sequence.txt.gz
fastq_2=$DATA_DIR/${sample_id}_2_sequence.txt.gz

bowtie2 -x /omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/rtta_alignment/rtta_index \
    -1 $fastq_1 \
    -2 $fastq_2 \
    -S /omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/rtta_alignment/output/exp4/${sample_id}.rtta.sam
