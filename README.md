# Escape from X inactivation is directly modulated by levels of Xist non-coding RNA

This repository contains all code to recreate the figures for the publication "Escape from X inactivation is directly modulated by levels of Xist non-coding RNA" (https://www.biorxiv.org/content/10.1101/2024.02.22.581559v1).

# Structure of the repository

## Data availability
All data (RNA-Seq, Cut and Run, WGBS and capture HIC) is available from GEO at XXX. 

## Preprocessing
- RNA-Seq: Preprocessing of RNA-Seq data into (allele-specific) count matrices can be done using the workflow at https://github.com/yuviaapr/allele-specific_RNA-seq. These matrices are then read into SCE objects using the scripts in ./Preprocessing/, which are the base for the rest of the analysis.
- Cut-and-Run: Preprocessing of Cut-and-Run data into (allele-specific) alignment files can be done using the workflow at https://github.com/yuviaapr/allele-specific_CUTandRUN, from which one can generate (allele-specific) count matrices for windows, genes and other genomic features (./Preprocessing).
- WGBS: WGBS data is aligned using the methylseq pipeline (https://nf-co.re/methylseq/2.6.0/) and then preprocessed into (allele-specific) per-base methylation value files (file names?) using SNP-split and (which tool?)
- Hi-C data: XXX

## Analysis
- 
