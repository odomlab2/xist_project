### create dataset for Binomial factor analysis

setwd("~/Desktop/Projects/XChromosome_Antonia/")

library(tidyverse)
library(ggplot2)
library(scran)
library(scater)
library(DESeq2)

source("./Scripts/auxiliary.R")

data <- readRDS("./ProcessedData/merged_dataset.rds")
data <- data[rowMeans(counts(data)) > 100, ]
# data <- data[,data$Clone == "E6"]

## 
counts_total <- counts(data)

counts_total_deseq <- DESeqDataSetFromMatrix(countData = counts(data), colData = colData(data), design = ~1)
counts_total_deseq <- logNormCounts(counts_total_deseq)

counts_total_lognorm <- assays(counts_total_deseq)$logcounts
counts_total_lognorm_hvg <- counts_total_lognorm[order(rowVars(counts_total_lognorm), decreasing = T)[1:1000], ]

pca_res <- runPCA(t(counts_total_lognorm_hvg), rank = 10)
pca_res_rotation <- pca_res$rotation

data.frame(
  PC1 = pca_res$x[,1], 
  PC2 = pca_res$x[,2], 
  col = data$Clone
) %>%
  ggplot(aes(PC1, PC2, col = col)) + geom_point()

## 

counts_reference <- counts_active(data)
counts_alternative <- counts_inactive(data)

counts_reference_hvg <- counts_reference[order(rowVars(as.matrix(counts_reference / (counts_reference + counts_alternative))), decreasing = T)[1:1000], ]
counts_alternative_hvg <- counts_alternative[order(rowVars(as.matrix(counts_reference / (counts_reference + counts_alternative))), decreasing = T)[1:1000], ]
table(as.character(seqnames(rowRanges(data)[rownames(counts_reference_hvg), ])))

pca_res_bino <- runPCA(t(counts_reference_hvg / (counts_reference_hvg + counts_alternative_hvg)), rank = 10)
pca_res_bino_rotations <- pca_res_bino$rotation

pca_res_bino_rotations[order(abs(pca_res_bino_rotations[,2]), decreasing = T), ] %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  add_column("chromosome" = as.character(seqnames(rowRanges(data)[.$gene, ]))) %>%
  head(n = 30)

data.frame(
  PC1 = pca_res_bino$x[,1], 
  PC2 = pca_res_bino$x[,2], 
  col = data$Clone, 
  label = data$CorrectedSampleName
) %>%
  ggplot(aes(PC1, PC2, col = col)) + geom_point()


## write data
write.csv(counts(data)[rownames(counts_total_lognorm_hvg), ], "./ProcessedData/total_expression_bulk_new.csv")

write.csv(counts_reference_hvg, "./ProcessedData/reference_expression_bulk.csv")
write.csv(counts_alternative_hvg, "./ProcessedData/alternative_expression_bulk.csv")


