# Escape from X inactivation is directly modulated by levels of Xist non-coding RNA

This repository contains all code to recreate the figures for the publication "Escape from X inactivation is directly modulated by levels of Xist non-coding RNA" (https://www.biorxiv.org/content/10.1101/2024.02.22.581559v1).

# Structure of the repository

## Data availability
All data (RNA-Seq, Cut and Run, WGBS and capture HIC) is available from GEO at GSE259400. 

## Preprocessing
- RNA-Seq: Preprocessing of RNA-Seq data into (allele-specific) count matrices can be done using the workflow at https://github.com/yuviaapr/allele-specific_RNA-seq. These matrices are then read into SCE objects using the scripts in ./Preprocessing/, which are the base for the rest of the analysis (generate_sce_XXX.R).
- Cut-and-Run: Preprocessing of Cut-and-Run data into (allele-specific) alignment files can be done using the workflow at https://github.com/yuviaapr/allele-specific_CUTandRUN, from which one can generate (allele-specific) count matrices for windows, genes and other genomic features (./Preprocessing).
- WGBS: WGBS data is aligned using the methylseq pipeline (https://nf-co.re/methylseq/2.6.0/) and then preprocessed into (allele-specific) per-base methylation value files (file names?) using SNP-split and (which tool?)
- Hi-C data: XXX

## Analysis
- Figure 1 and Figures S1-4 shows the analysis of escape across clonal cell lines, as well as the analysis of escape after Xist-overexpression, the analysis can be reproduced using figure_1_escapee_silencing.Rmd and figure_1_multiclone_analysis.Rmd
- Figure 2 and Figures SX show the analysis of escape in astrocytes after Xist overexpression (figure_2_astrocytes.RmD)
- Figure 3 and Figures SX show the dependency of Xist-mediated silenced to SPEN (figure_3_spen_knockdown.Rmd)
- Figure 4 ...
- Figure 6 shows the reversibility of escape after stopping Xist-overexpression (figure_6_final_script.Rmd)
- Figure 7 shows the analysis of H3K4me3 and H2AK119ub1 as well as DNA methylation marks after Xist overexpression (figure_7_cnr_met_coverage.Rmd, figure_7_cnr_met.Rmd)
- Figure 8 shows the analysis of escape after Xist overexpression in embryos (...)
