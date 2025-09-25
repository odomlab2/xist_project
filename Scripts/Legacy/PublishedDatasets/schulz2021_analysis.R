library(tidyverse)
library(SingleCellExperiment)
library(scran)
library(scater)
library(scuttle)

setwd("~/Desktop/Projects/XChromosome_Antonia/")

data.b6 <- read_tsv("~/Desktop/Projects/XChromosome_Antonia/Data/Schulz2021/GSE151009_B6_UMICountMatrix.txt", col_names = F, skip = 1) %>%
  column_to_rownames("X1")
colnames(data.b6) <- colnames(read_tsv("~/Desktop/Projects/XChromosome_Antonia/Data/Schulz2021/GSE151009_B6_UMICountMatrix.txt", n_max = 0))

data.cast <- read_tsv("~/Desktop/Projects/XChromosome_Antonia/Data/Schulz2021/GSE151009_Cast_UMICountMatrix.txt", col_names = F, skip = 1) %>%
  column_to_rownames("X1")
colnames(data.cast) <- colnames(read_tsv("~/Desktop/Projects/XChromosome_Antonia/Data/Schulz2021/GSE151009_Cast_UMICountMatrix.txt", n_max = 0))

data.total <- data.b6 + data.cast


### 
# testy <- read.csv("~/Desktop/Projects/XChromosome_Antonia/Data/Reinius2022/processed_diff/mESC_diff_clone5.fract_CAST_reads.txt", sep = "\t", row.names = 1)

metadata <- data.frame(
  row.names = colnames(data.total),
  differentiation_day = unlist(lapply(str_split(colnames(data.total), "_"), function(x){x[[1]]}))
)

# isolate genes information
genes <- rownames(data.total)

# add chromosome information to genes
library("EnsDb.Mmusculus.v79")
gene_info <- data.frame(genes(EnsDb.Mmusculus.v79))
rownames(gene_info) <- gene_info$gene_id
gene_info <- gene_info[!duplicated(gene_info$symbol), ]

genes_keep <- intersect(genes, gene_info$symbol)
gene_info <- gene_info[gene_info$symbol %in% genes_keep, ]
genes <- genes[genes %in% genes_keep]

genes <- gene_info[match(genes, gene_info$symbol), ]
genes <- makeGRangesFromDataFrame(genes, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = T)

# parse into SCE object 
sce_dataset1 <- SingleCellExperiment(assays = list("counts" = data.total[genes$symbol, ], 
                                                   "counts_b6" = data.b6[genes$symbol, ], 
                                                   "counts_cast" = data.cast[genes$symbol, ]), 
                                     colData = DataFrame(metadata), rowData = genes)

# exclude cells with few reads
hist(log10(colSums(counts(sce_dataset1))))

sce_dataset1 <- sce_dataset1[,colSums(counts(sce_dataset1)) > 10000]
sce_dataset1 <- logNormCounts(sce_dataset1)
set.seed(1234)
sce_dataset1 <- runPCA(sce_dataset1)
plotPCA(sce_dataset1, colour_by = "differentiation_day")
plotPCA(sce_dataset1, colour_by = "Pitrm1")
rownames(sce_dataset1) <- rowData(sce_dataset1)$gene_name

pca_here <- prcomp(t(logcounts(sce_dataset1)))

genes_rel = names(sort(rowSums(pca_here$rotation[,1:20] ** 2 ), decreasing = T)[1:100])
rowData(sce_dataset1[genes_rel, ])$gene_name
plotPCA(sce_dataset1, colour_by = "Krt18")

# look at allelic ratios
rowData(sce_dataset1)$total_counts <- rowSums(counts(sce_dataset1))
rowData(sce_dataset1)$chromosome <- seqnames(rowRanges(sce_dataset1))

rowData(sce_dataset1) %>%
  data.frame() %>%
  dplyr::filter(!is.na(chromosome)) %>%
  dplyr::filter(chromosome == "X") %>%
  dplyr::arrange(-total_counts, decreasing = T) %>%
  head()

data.frame(
  PC1 = reducedDims(sce_dataset1)[["PCA"]][,1], 
  PC2 = reducedDims(sce_dataset1)[["PCA"]][,2], 
  cov = as.numeric(allelic_ratio)
) %>%
  ggplot(aes(x = PC1, y = PC2, col = cov)) + geom_point() + viridis::scale_color_viridis()

gene = "Car2"
allelic_ratio = assays(sce_dataset1[gene, ])[["counts_b6"]] / (assays(sce_dataset1[gene, ])[["counts_b6"]] + assays(sce_dataset1[gene, ])[["counts_cast"]])
data.frame(
  cond = sce_dataset1$differentiation_day,
  cov = as.numeric(allelic_ratio)
) %>%
  ggplot(aes(x = cond, y = cov)) + stat_summary(size = 3) + geom_jitter() + viridis::scale_color_viridis()

# subset on the X chromosome
sce_dataset1.x <- sce_dataset1[seqnames(rowRanges(sce_dataset1)) == "X", ]

x_ratios <- colSums(assays(sce_dataset1.x)[["counts_b6"]]) / ( colSums(assays(sce_dataset1.x)[["counts_b6"]]) +  colSums(assays(sce_dataset1.x)[["counts_cast"]]))
xist_ratios <- colSums(assays(sce_dataset1.x["Xist", ])[["counts_b6"]]) / 
  ( colSums(assays(sce_dataset1.x["Xist", ])[["counts_b6"]]) +  colSums(assays(sce_dataset1.x["Xist", ])[["counts_cast"]]))

plot(x_ratios, xist_ratios)

data.frame(
  PC1 = reducedDims(sce_dataset1)[["PCA"]][,1], 
  PC2 = reducedDims(sce_dataset1)[["PCA"]][,2], 
  x_ratio = x_ratios
) %>%
  ggplot(aes(x = PC1, y = PC2, col = x_ratio)) + geom_point() + 
  scale_color_viridis_c()

data.frame(
  day = sce_dataset1$differentiation_day, 
  x_ratio = x_ratios
) %>%
  ggplot(aes(x = factor(day), y = x_ratio)) + ggbeeswarm::geom_quasirandom()

data.frame(
  day = sce_dataset1$differentiation_day, 
  x_ratio = x_ratios
) %>%
  ggplot(aes(x =  x_ratio)) + geom_histogram(bins = 100) + 
  geom_vline(xintercept = 0.2, linetype = "dashed") + 
  geom_vline(xintercept = 0.4, linetype = "dashed") + 
  geom_vline(xintercept = 0.6, linetype = "dashed") + 
  geom_vline(xintercept = 0.8, linetype = "dashed")

x_ratios_intervals <- cut(x_ratios, breaks = 4)
convert_intervals <- unlist(setNames(c("CAST", "unclear", "unclear", "B6"), levels(x_ratios_intervals))[x_ratios_intervals])

data.frame(
  PC1 = reducedDims(sce_dataset1)[["PCA"]][,1], 
  PC2 = reducedDims(sce_dataset1)[["PCA"]][,2], 
  x_ratio = convert_intervals, 
  days = factor(sce_dataset1$differentiation_day)
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
escapee_annotation <- readxl::read_excel("~/Desktop/Projects/XChromosome_Project/ProcessedData/ListOfEscapeeFromEdithLab.xlsx") %>%
  dplyr::select(c("ENSEMBL_v102", "final status"))
library("EnsDb.Mmusculus.v79") 
symbols <- mapIds(EnsDb.Mmusculus.v79, keys = escapee_annotation$ENSEMBL_v102, keytype = "GENEID", column="SYMBOL")
escapee_annotation$symbol <- symbols
escapee_annotation <- escapee_annotation %>% column_to_rownames("ENSEMBL_v102")

escape_data <- data.frame(
  per_gene_active = rowSums(assays(sce_dataset1_subset.x)[["active_counts"]]), 
  per_gene_inactive = rowSums(assays(sce_dataset1_subset.x)[["inactive_counts"]]), 
  gene =  rowData(sce_dataset1_subset.x)$gene_name
) %>%
  add_column(ratio = .$per_gene_inactive / (.$per_gene_inactive + .$per_gene_active)) %>%
  add_column(escape_annotation = escapee_annotation[rownames(.), 1])

escape_data %>%
  dplyr::filter(per_gene_active + per_gene_inactive > 100) %>%
  dplyr::filter(ratio > 0.1 & ratio < 0.5) %>%
  dplyr::arrange(ratio)

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

## 
genes_use <- escape_data %>%
  dplyr::filter(per_gene_active + per_gene_inactive > 500) %>%
  dplyr::filter(ratio > 0.1 & ratio < 0.5) %>%
  dplyr::arrange(ratio)

sce_dataset1_subset.x <- sce_dataset1_subset.x[,sce_dataset1_subset.x$differentiation_day %in% c("d3", "d4")]
ratios_here <- assays(sce_dataset1_subset.x)[["inactive_counts"]] / (assays(sce_dataset1_subset.x)[["inactive_counts"]] + assays(sce_dataset1_subset.x)[["active_counts"]])
ratios_here[assays(sce_dataset1_subset.x)[["inactive_counts"]] + assays(sce_dataset1_subset.x)[["active_counts"]] <= 2] <- 0
ratios_here <- ratios_here[rownames(genes_use), ]
ratios_here <- ratios_here[rownames(genes_use), ]
rownames(ratios_here) <- genes_use$gene

pheatmap::pheatmap(ratios_here, cluster_rows = F)

pca_here <- prcomp(t(ratios_here))
pca_here_rotations <- data.frame(pca_here$rotation[,1:5])
head(pca_here_rotations[order(pca_here_rotations$PC2, decreasing = T), ])

data.frame(
  PC1 = pca_here$x[,1],
  PC2 = pca_here$x[,2], 
  col = sce_dataset1_subset.x$assigned_x
) %>%
  ggplot(aes(x = PC1, y = PC2, col = col)) + geom_point()



