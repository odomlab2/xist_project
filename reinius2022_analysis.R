library(tidyverse)
library(SingleCellExperiment)
library(scran)
library(scater)
library(scuttle)

setwd("~/Desktop/Projects/XChromosome_Antonia/")

data.b6 <- read.csv("~/Desktop/Projects/XChromosome_Antonia/Data/Reinius2022/processed_diff/mESC_diff_clone5.BL6_reads.txt", sep = "\t", row.names = 1)
data.cast <- read.csv("~/Desktop/Projects/XChromosome_Antonia/Data/Reinius2022/processed_diff/mESC_diff_clone5.CAST_reads.txt", sep = "\t", row.names = 1)
data.total <- data.b6 + data.cast

testy <- read.csv("~/Desktop/Projects/XChromosome_Antonia/Data/Reinius2022/processed_diff/mESC_diff_clone5.fract_CAST_reads.txt", sep = "\t", row.names = 1)

metadata <- read.csv("./Data/Reinius2022/processed_diff/metadata_ss3_clone5diff.tsv", sep = "\t") %>%
  column_to_rownames("sample_bc")
metadata <- metadata[colnames(data.total), ]

# isolate genes information
genes <- rownames(data.total)

# add chromosome information to genes
library("EnsDb.Mmusculus.v79")
gene_info <- data.frame(genes(EnsDb.Mmusculus.v79))
rownames(gene_info) <- gene_info$gene_id
gene_info <- gene_info[!duplicated(gene_info$symbol), ]

genes_keep <- intersect(genes, gene_info$gene_id)
gene_info <- gene_info[genes_keep, ]
genes <- genes[genes %in% genes_keep]

genes <- gene_info[genes, ]
genes <- makeGRangesFromDataFrame(genes, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = T)

# parse into SCE object 
sce_dataset1 <- SingleCellExperiment(assays = list("counts" = data.total[names(genes), ], 
                                                   "counts_b6" = data.b6[names(genes), ], 
                                                   "counts_cast" = data.cast[names(genes), ]), 
                                     colData = DataFrame(metadata), rowData = genes)

# exclude cells with few reads
hist(log10(colSums(counts(sce_dataset1))))

sce_dataset1 <- sce_dataset1[,colSums(counts(sce_dataset1)) > 1000]
sce_dataset1 <- logNormCounts(sce_dataset1)
set.seed(1234)
sce_dataset1 <- runPCA(sce_dataset1)
plotPCA(sce_dataset1, colour_by = "day")

# subset on the X chromosome
sce_dataset1.x <- sce_dataset1[seqnames(rowRanges(sce_dataset1)) == "X", ]

x_ratios <- colSums(assays(sce_dataset1.x)[["counts_b6"]]) / ( colSums(assays(sce_dataset1.x)[["counts_b6"]]) +  colSums(assays(sce_dataset1.x)[["counts_cast"]]))
xist_ratios <- colSums(assays(sce_dataset1.x["ENSMUSG00000086503", ])[["counts_b6"]]) / 
  ( colSums(assays(sce_dataset1.x["ENSMUSG00000086503", ])[["counts_b6"]]) +  colSums(assays(sce_dataset1.x["ENSMUSG00000086503", ])[["counts_cast"]]))

data.frame(
  PC1 = reducedDims(sce_dataset1)[["PCA"]][,1], 
  PC2 = reducedDims(sce_dataset1)[["PCA"]][,2], 
  x_ratio = x_ratios
) %>%
  ggplot(aes(x = PC1, y = PC2, col = x_ratio)) + geom_point() + 
  scale_color_viridis_c()

data.frame(
  day = sce_dataset1$day, 
  x_ratio = x_ratios
) %>%
  ggplot(aes(x = factor(day), y = x_ratio)) + ggbeeswarm::geom_quasirandom()

data.frame(
  day = sce_dataset1$day, 
  x_ratio = x_ratios
) %>%
  ggplot(aes(x =  x_ratio)) + geom_histogram(bins = 100) + 
  geom_vline(xintercept = 0.2, linetype = "dashed") + 
  geom_vline(xintercept = 0.4, linetype = "dashed") + 
  geom_vline(xintercept = 0.6, linetype = "dashed") + 
  geom_vline(xintercept = 0.8, linetype = "dashed")

x_ratios_intervals <- cut(x_ratios, breaks = 5)
convert_intervals <- unlist(setNames(c("CAST", "unclear", "dual", "unclear", "B6"), levels(x_ratios_intervals))[x_ratios_intervals])

data.frame(
  PC1 = reducedDims(sce_dataset1)[["PCA"]][,1], 
  PC2 = reducedDims(sce_dataset1)[["PCA"]][,2], 
  x_ratio = convert_intervals, 
  days = factor(sce_dataset1$day)
) %>%
  ggplot(aes(x = PC1, y = PC2, col = x_ratio, shape = days)) + geom_point()

# 
sce_dataset1$assigned_x <- convert_intervals
sce_dataset1_subset <- sce_dataset1[,sce_dataset1$assigned_x %in% c("CAST", "B6")]

active_counts <- cbind(assays(sce_dataset1_subset[,sce_dataset1_subset$assigned_x == "B6"])[["counts_b6"]],
                       assays(sce_dataset1_subset[,sce_dataset1_subset$assigned_x == "CAST"])[["counts_cast"]])

inactive_counts <- cbind(assays(sce_dataset1_subset[,sce_dataset1_subset$assigned_x == "B6"])[["counts_cast"]], 
                         assays(sce_dataset1_subset[,sce_dataset1_subset$assigned_x == "CAST"])[["counts_b6"]])

assays(sce_dataset1_subset)[["active_counts"]] <- active_counts[,colnames(sce_dataset1_subset)]
assays(sce_dataset1_subset)[["inactive_counts"]] <- inactive_counts[,colnames(sce_dataset1_subset)]

sce_dataset1_subset.x <- sce_dataset1_subset[seqnames(rowRanges(sce_dataset1_subset)) == "X", ]

# look at escapees

escape_data <- data.frame(
  per_gene_active = rowSums(assays(sce_dataset1_subset.x)[["active_counts"]]), 
  per_gene_inactive = rowSums(assays(sce_dataset1_subset.x)[["inactive_counts"]]), 
  gene =  rowData(sce_dataset1_subset.x)$gene_name
) %>%
  add_column(ratio = .$per_gene_inactive / (.$per_gene_inactive + .$per_gene_active))

escape_data %>%
  dplyr::filter(per_gene_active + per_gene_inactive > 100) %>%
  dplyr::arrange(ratio) %>%
  tail()

data.frame(
  per_gene_active = rowSums(assays(sce_dataset1_subset.x)[["active_counts"]]), 
  per_gene_inactive = rowSums(assays(sce_dataset1_subset.x)[["inactive_counts"]]), 
  gene =  rowData(sce_dataset1_subset.x)$gene_name
) %>%
  ggplot(aes(x = per_gene_active + per_gene_inactive, y = per_gene_inactive / (per_gene_active + per_gene_inactive))) + 
  geom_point() + scale_x_log10()

data.frame(
  per_gene_active = rowSums(assays(sce_dataset1_subset.x)[["active_counts"]]), 
  per_gene_inactive = rowSums(assays(sce_dataset1_subset.x)[["inactive_counts"]]), 
  gene =  rowData(sce_dataset1_subset.x)$gene_name
) %>%
  ggplot(aes(x = per_gene_active + per_gene_inactive, y = per_gene_inactive / (per_gene_active + per_gene_inactive))) + 
  geom_point() + scale_x_log10() + 
  ggrepel::geom_text_repel(aes(label = gene))
  






