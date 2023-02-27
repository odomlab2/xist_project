### Preprocessing and generation of expression objects

library(tidyverse)
library(ggplot2)

source("./Scripts/auxiliary.R")

setwd("~/Desktop/Projects/XChromosome_Antonia/")

data <- readRDS("./ProcessedData/merged_dataset.rds")

# Run some QC: 
# How many total reads / sample
data.frame(counts = colSums(counts(data))) %>%
  rownames_to_column("Sample") %>%
  ggplot(aes(x = Sample, y = counts)) + geom_bar(stat = "identity") + 
  theme_paper() + coord_flip() + 
  ggtitle("Total Counts (All)")

# How many reads b6 / sample
data.frame(counts = colSums(counts_inactive(data))) %>%
  rownames_to_column("Sample") %>%
  ggplot(aes(x = Sample, y = counts)) + geom_bar(stat = "identity") + 
  theme_paper() + coord_flip() + 
  ggtitle("Total Counts (B6)")

# How many reads cast / sample
data.frame(counts = colSums(counts_active(data))) %>%
  rownames_to_column("Sample") %>%
  ggplot(aes(x = Sample, y = counts)) + geom_bar(stat = "identity") + 
  theme_paper() + coord_flip() + 
  ggtitle("Total Counts (Cast)")

# Ratios per sample?
data.frame(ratio =  colSums(counts_inactive(data)) / (colSums(counts_inactive(data)) + colSums(counts_active(data))) ) %>%
  rownames_to_column("Sample") %>%
  ggplot(aes(x = Sample, y = ratio)) + geom_bar(stat = "identity") + 
  theme_paper() + coord_flip() + 
  ggtitle("Allelic ratios (Autosomes + X)") + ylim(0, 1) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "red")

# look at counts per gene
data.frame(counts = rowSums(counts(data))) %>%
  ggplot(aes(x = counts)) + geom_histogram(bins = 100) + 
  scale_x_log10() + theme_paper()

# 
data.frame(
  B6_Counts = rowSums(counts_inactive(data)), 
  CAST_Counts = rowSums(counts_active(data)), 
  chromosome = seqnames(rowRanges(data))
) %>%
  ggplot(aes(x = B6_Counts + 1, y = CAST_Counts + 1, col = chromosome == "X")) + geom_point() + 
  theme_paper() + scale_x_log10() + scale_y_log10()

# Check Xist expression
data.frame(
  xist_b6 = as.numeric(counts_inactive(data["Xist", ])), 
  xist_cast = as.numeric(counts_active(data["Xist", ])), 
  samples = data$SampleName
) %>% ggplot(aes(x = samples, y = xist_b6 / (xist_b6 + xist_cast))) + geom_point() + coord_flip() + 
  ylim(c(0, 1))

# Look at autosomal gene expression patterns
library(DESeq2)

des <- DESeqDataSetFromMatrix(counts(data), 
                              colData = colData(data), 
                              design = ~ 1)
des.vst <- vst(des)
rowData(des.vst)$symbol <- rowData(data)$gene_name

var_data <- rowData(des.vst)
data.frame(head(var_data[order(var_data$baseVar, decreasing = T), ], n = 50))

plotPCA(des.vst, intgroup = "Condition")
# 


