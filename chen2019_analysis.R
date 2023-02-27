library(tidyverse)
library(SingleCellExperiment)
library(scran)
library(scater)
library(scuttle)

setwd("~/Desktop/Projects/XChromosome_Antonia/")

data <- readRDS("./Data/Chen2019/raw_data_matrices.rds")
metadata <- read_csv("./Data/Chen2019/SraRunTable.txt")

# subset on mice with F1 genotype
metadata$Genotype <- gsub("C57BL/6J \\(F\\)x CAST/EiJ \\(M\\)", "C57BL/6J \\(F\\) x CAST/EiJ \\(M\\)", metadata$Genotype)

data_b6 <- data[[1]]
data_cast <- data[[2]]
data_undet <- data[[3]]

colnames(data_b6) <- gsub("_R", "", colnames(data_b6))
colnames(data_cast) <- gsub("_A", "", colnames(data_cast))
colnames(data_undet) <- gsub("_U", "", colnames(data_undet))

data_total <- data_b6 + data_cast + data_undet

metadata <- metadata %>%
  dplyr::filter(Run %in% colnames(data_total)) %>%
  column_to_rownames("Run")

# add chromosome information to genes
library("EnsDb.Mmusculus.v79")
gene_info <- data.frame(genes(EnsDb.Mmusculus.v79))
rownames(gene_info) <- gene_info$gene_id
gene_info <- gene_info[!duplicated(gene_info$symbol), ]

genes_keep <- intersect(rownames(data_total), gene_info$symbol)
gene_info <- gene_info[gene_info$symbol %in% genes_keep, ]
genes <- rownames(data_total)[rownames(data_total) %in% genes_keep]

genes <- gene_info[match(genes, gene_info$symbol), ]
genes <- makeGRangesFromDataFrame(genes, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = T)

# parse into SCE object 
sce_dataset1 <- SingleCellExperiment(assays = list("counts" = data_total[genes$symbol, rownames(metadata)], 
                                                   "counts_b6" = data_b6[genes$symbol, rownames(metadata)], 
                                                   "counts_cast" = data_cast[genes$symbol, rownames(metadata)]), 
                                     colData = DataFrame(metadata), rowData = genes)

rownames(sce_dataset1) <- genes$symbol

# exclude cells with few reads
hist(log10(colSums(counts(sce_dataset1))))
sce_dataset1 <- sce_dataset1[,colSums(counts(sce_dataset1)) > 10000]

# exclude non-F1 cells
sce_dataset1 <- sce_dataset1[,!is.na(sce_dataset1$Genotype)]
sce_dataset1 <- sce_dataset1[,!sce_dataset1$Genotype == "C57BL/6J (F) x C57BL/6J (M)"]

sce_dataset1 <- logNormCounts(sce_dataset1)
set.seed(1234)
sce_dataset1 <- runPCA(sce_dataset1)
sce_dataset1 <- runTSNE(sce_dataset1)

sce_dataset1_hvgs <- getTopHVGs(sce_dataset1)

# plotPCA(sce_dataset1, colour_by = "Age")
# plotPCA(sce_dataset1, colour_by = "embryo")
# plotPCA(sce_dataset1, colour_by = "Genotype")
# 
# plotTSNE(sce_dataset1, colour_by = "Age")
# plotTSNE(sce_dataset1, colour_by = "embryo")
# plotTSNE(sce_dataset1, colour_by = "Genotype")
# plotTSNE(sce_dataset1, colour_by = "clusters")

sce.reference <- MouseGastrulationData::EmbryoAtlasData(samples = 21)
sce.reference <- logNormCounts(sce.reference)
rownames(sce.reference) <- rowData(sce.reference)$SYMBOL

library(SingleR)
sce.all.annotated <- SingleR(test=sce_dataset1[sce_dataset1_hvgs, ], ref=sce.reference, labels=colData(sce.reference)$celltype)
sce_dataset1$SingleR_celltype <- sce.all.annotated$labels

plotPCA(sce_dataset1, colour_by = "SingleR_celltype")
plotTSNE(sce_dataset1, colour_by = "SingleR_celltype")

# annotate sex
data.frame(
  FemaleSpecific = logcounts(sce_dataset1)["Xist", ], 
  MaleSpecific = logcounts(sce_dataset1)["Eif2s3y", ], 
  cov = sce_dataset1$Genotype
) %>%
  ggplot(aes(x = FemaleSpecific, MaleSpecific, col = cov)) + geom_point()

df_here <- data.frame(
  FemaleSpecific = logcounts(sce_dataset1)["Xist", ], 
  MaleSpecific = logcounts(sce_dataset1)["Eif2s3y", ]
) %>% 
  add_column("male" = .$MaleSpecific > 2) %>%
  add_column("female" = .$FemaleSpecific > 2)

sce_dataset1$female_signal <- df_here$female
sce_dataset1$male_signal <- df_here$male

sce_dataset1 <- sce_dataset1[,!sce_dataset1$male_signal & sce_dataset1$female_signal]

# annotate celltypes 

# subset on the X chromosome
sce_dataset1.x <- sce_dataset1[seqnames(rowRanges(sce_dataset1)) == "X", ]
rownames(sce_dataset1.x) <- rowData(sce_dataset1.x)$gene_name

x_ratios <- colSums(assays(sce_dataset1.x)[["counts_b6"]]) / ( colSums(assays(sce_dataset1.x)[["counts_b6"]]) +  colSums(assays(sce_dataset1.x)[["counts_cast"]]))
xist_ratios <- colSums(assays(sce_dataset1.x["Xist", ])[["counts_b6"]]) / 
  ( colSums(assays(sce_dataset1.x["Xist", ])[["counts_b6"]]) +  colSums(assays(sce_dataset1.x["Xist", ])[["counts_cast"]]))
x_ratios_noxist <- 
  colSums(assays(sce_dataset1.x[rownames(sce_dataset1.x) != "Xist", ])[["counts_b6"]]) / 
  ( colSums(assays(sce_dataset1.x[rownames(sce_dataset1.x) != "Xist", ])[["counts_b6"]]) +  
      colSums(assays(sce_dataset1.x[rownames(sce_dataset1.x) != "Xist", ])[["counts_cast"]]))

data.frame(
  x_ratios = x_ratios_noxist, 
  xist_ratios = xist_ratios,
  celltype = sce_dataset1.x$SingleR_celltype
) %>%
  ggplot(aes(x = x_ratios, y = xist_ratios, col = celltype)) + geom_point()

data.frame(
  PC1 = reducedDims(sce_dataset1)[["PCA"]][,1], 
  PC2 = reducedDims(sce_dataset1)[["PCA"]][,2], 
  x_ratio = x_ratios
) %>%
  ggplot(aes(x = PC1, y = PC2, col = x_ratio)) + geom_point() + 
  scale_color_viridis_c()

data.frame(
  day = sce_dataset1$Age, 
  x_ratio = x_ratios_noxist
) %>%
  ggplot(aes(x = factor(day), y = x_ratio)) + ggbeeswarm::geom_quasirandom()

data.frame(
  day = sce_dataset1$SingleR_celltype, 
  x_ratio = x_ratios_noxist, 
  cross = sce_dataset1$Genotype
) %>%
  ggplot(aes(x = factor(day), y = x_ratio, col = cross)) + 
  ggbeeswarm::geom_quasirandom() +
  coord_flip() + facet_wrap(~cross)

data.frame(
  day = sce_dataset1$Age, 
  x_ratio = x_ratios
) %>%
  ggplot(aes(x =  x_ratio)) + geom_histogram(bins = 100) + 
  geom_vline(xintercept = 0.33, linetype = "dashed") + 
  geom_vline(xintercept = 0.4, linetype = "dashed") + 
  geom_vline(xintercept = 0.6, linetype = "dashed") + 
  geom_vline(xintercept = 0.66, linetype = "dashed")

x_ratios_intervals <- cut(x_ratios_noxist, breaks = 3)
convert_intervals <- unlist(setNames(c("CAST", "unclear", "B6"), levels(x_ratios_intervals))[x_ratios_intervals])

data.frame(
  PC1 = reducedDims(sce_dataset1)[["PCA"]][,1], 
  PC2 = reducedDims(sce_dataset1)[["PCA"]][,2], 
  x_ratio = convert_intervals, 
  days = factor(sce_dataset1$Age)
) %>%
  ggplot(aes(x = PC1, y = PC2, col = x_ratio, shape = days)) + geom_point()

# 
sce_dataset1$assigned_x <- convert_intervals
sce_dataset1_subset <- sce_dataset1[,sce_dataset1$assigned_x %in% c("CAST", "B6")]

sce_dataset1_subset <- sce_dataset1_subset[,!sce_dataset1_subset$SingleR_celltype %in% c("ExE ectoderm", "ExE endoderm", "Veisceral endoderm", "Parietal endoderm")]

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

#sce_dataset1_subset.x <- sce_dataset1_subset.x[,!sce_dataset1_subset.x$SingleR_celltype %in% c("ExE ectoderm", "ExE endoderm", "Parietal endoderm", "Visceral endoderm")]
sce_dataset1_subset.x <- sce_dataset1_subset.x[,sce_dataset1_subset.x$SingleR_celltype == 'Epiblast']
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
  col = sce_dataset1_subset.x$Genotype
) %>%
  ggplot(aes(x = PC1, y = PC2, col = col)) + geom_point()

data.frame(
  average_escape = colMeans(ratios_here), 
  xist_expression = logcounts(sce_dataset1_subset.x)["Xist", ]
) %>% 
  ggplot(aes(x = xist_expression, y = average_escape)) + geom_point() + 
    geom_smooth(method = "lm")

testy <- cor(t(ratios_here), logcounts(sce_dataset1_subset.x)["Xist", ])
testy[order(testy[,1]), ]


