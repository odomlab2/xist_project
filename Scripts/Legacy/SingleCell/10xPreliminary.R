library(DropletUtils)
library(scran)
library(reshape2)
library(scater)
library(ggplot2)
library(Rtsne)
library(irlba)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(tidyverse)
source("~/Desktop/PhD/RCode/Multipurpose/auxiliary.R")

setwd("~/Desktop/Projects/XChromosome_Antonia/")

#### Read individual libraries

sce.sample1 <- read10xCounts("~/Desktop/Projects/XChromosome_Antonia/Data/ALNPCs092022/")
colData(sce.sample1)$Sample <- rep("B6xCast", ncol(sce.sample1))
colData(sce.sample1)$Library <- rep("Sample1", ncol(sce.sample1))

#### merge dataset 
sce.all <- sce.sample1
rm("sce.sample1")

colnames(sce.all) <- paste0(sce.all$Library, "_", sce.all$Barcode)

# calculate number of detected cells per library
cur_stats <- melt(table(colData(sce.all)$Sample, colData(sce.all)$Library))
cur_stats <- cur_stats[cur_stats$value > 0,]
cur_stats <- cur_stats[order(cur_stats$Var1),]
stats.df <- data.frame(row.names = cur_stats$Var2,
                       Sample = cur_stats$Var1,
                       Library = cur_stats$Var2,
                       No_cells = cur_stats$value)

# ---------------------------------------------------------------------------------------------
ggplot(stats.df, aes(x = Library, y = No_cells, fill = Sample)) + 
  geom_bar(stat = "identity") + theme_classic() + 
  coord_flip() + ylim(0, 20000) + ggtitle("Number of cells before filtering")
# ---------------------------------------------------------------------------------------------

QC.metrics <- perCellQCMetrics(sce.all)
colData(sce.all) <- cbind(colData(sce.all), QC.metrics)

QC.metrics %>% data.frame() %>%
  ggplot(aes(x = 1, y = total)) + geom_violin() + scale_y_log10() # need to filter for cells > ~10k

QC.metrics %>% data.frame() %>%
  ggplot(aes(x = 1, y = detected)) + geom_violin() + scale_y_log10() # need to filter for cells > ~3k

sce.all <- sce.all[,colData(sce.all)$detected > 2500]
sce.all <- sce.all[,colData(sce.all)$total > 5000]
# Remove genes that are not expressed
sce.all <- sce.all[Matrix::rowSums(counts(sce.all)) > 0,]
# Add to stats data frame
cur_stats <- melt(table(colData(sce.all)$Sample, colData(sce.all)$Library))
cur_stats <- cur_stats[cur_stats$value > 0,]
cur_stats <- cur_stats[order(cur_stats$Var1),]
stats.df$AfterFiltering <- cur_stats$value

# ---------------------------------------------------------------------------------------------
ggplot(stats.df, aes(x = Library, y = AfterFiltering, fill = Sample)) + 
  geom_bar(stat = "identity") + theme_classic() + 
  coord_flip() + ylim(0, 10000) + ggtitle("Number of cells after filtering")
# ---------------------------------------------------------------------------------------------

rownames(sce.all) <- rowData(sce.all)$Symbol

clusters <- quickCluster(sce.all, method = "igraph",use.ranks=FALSE, min.size = 100)
sce.all <- computeSumFactors(sce.all, clusters=clusters)
sce.all <- logNormCounts(sce.all)

HVgenes <- HVG(sce = sce.all, n = 1000)
HVgenes[1:50]

pca <- prcomp_irlba(t(logcounts(sce.all[HVgenes,])), n = 50)
tsne <- Rtsne(pca$x, pca = FALSE)
reducedDims(sce.all)$PCA <- pca$x
reducedDims(sce.all)$TSNE <- tsne$Y

# ---------------------------------------------------------------------------------------------
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  batch = paste(colData(sce.all)$Sample, 
                                colData(sce.all)$Library))) +
  geom_point(aes(tsne1, tsne2), alpha = 0.6)
# ---------------------------------------------------------------------------------------------

# Number genes expressed
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = colData(sce.all)$detected)) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle("Number of genes expressed") + theme_bw()

# Plot gene expression
gene = "Nes"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene) + theme_bw()

gene = "Pax6"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene) + theme_bw()

gene = "Vim"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene) + theme_bw()

gene = "Olig2"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene) + theme_bw()


# Plot gene expression
gene = "Hist1h1b"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene) + theme_bw()

gene = "Hmgb2"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene) + theme_bw()

gene = "Top2a"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene) + theme_bw()

gene = "Id3"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene) + theme_bw()

gene = "Pou3f1"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene) + theme_bw()

gene = "Nav2"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene) + theme_bw()

gene = "Mki67"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene) + theme_bw()

gene = "Xist"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene) + theme_bw()

gene = "Kdm5c"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene) + theme_bw()

# annotate cell cycle stage
library(scran)
library(AnnotationDbi)
library(org.Mm.eg.db)

hs.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
hs.pairs.symbol <- lapply(hs.pairs, function(pp){
  pp.1 = mapIds(org.Mm.eg.db, keys = pp$first, column="SYMBOL",keytype="ENSEMBL",multiVals="first")
  pp.2 = mapIds(org.Mm.eg.db, keys = pp$second, column="SYMBOL",keytype="ENSEMBL",multiVals="first")
  pp.df = data.frame(
    first = pp.1,
    second = pp.2
  )
  pp.df[!duplicated(pp.df),]
})

cc.assignment <- cyclone(as.matrix(logcounts(sce.all)), hs.pairs.symbol)
sce.all$cycle_stage <- cc.assignment$phases

ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = sce.all$cycle_stage))  +
  geom_point(aes(tsne1, tsne2, colour = gene)) + 
  ggtitle("Cell cycle phase") + theme_bw()

# ---------------------------------------------------------------------------------------------

# clusters batch corrected data
set.seed(1234)
clusters <- quickCluster(sce.all)
sce.all$clusters <- clusters

##### annotate clusters
marker.genes <- findMarkers(sce.all, clusters)
marker.genes <- lapply(marker.genes, function(x){x[x$summary.logFC > 0,]})

mito.proportion <- colCounts(as.matrix(counts(sce.all[grepl("^mt-", rownames(sce.all)),]))) / 
  colCounts(as.matrix(counts(sce.all)))

ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  mito_prop = mito.proportion)) +
  geom_point(aes(tsne1, tsne2, colour = mito_prop)) + 
  scale_color_viridis()

head(data.frame(marker.genes$`13`), n = 30)

# ---------------------------------------------------------------------------------------------
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  cluster = as.factor(clusters))) +
  geom_point(aes(tsne1, tsne2, colour = cluster))
# ---------------------------------------------------------------------------------------------

### add allelic data

### question: allelic coverage per cell and gene

source("~/Desktop/PhD/Projects/ASE_Spermatogenesis_Paper/Scripts/General/auxiliary.R")
source("~/Desktop/PhD/Projects/ASE_Spermatogenesis_Paper/Scripts/General/ase_functions.R")
source("~/Desktop/PhD/Projects/ASE_Spermatogenesis_Paper/Scripts/General/reuse_functions.R")

### read ase-data
data.sample1 <- read_10x_ase("./Data/ALNPCs092022/ase_feature_matrix/")

add_empty_rows <- function(data, vector){
  data = data[rownames(data) %in% vector,]
  new_rows = setdiff(vector, rownames(data))
  add_zeros = matrix(rep(0, length(new_rows) * ncol(data)), 
                     nrow = length(new_rows),
                     ncol = ncol(data))
  rownames(add_zeros) = new_rows
  dd = rbind(data, add_zeros)
  #rownames(dd) = c(rownames(data), new_rows)
  dd[vector, ]
  # dd_out <- do.call("rbind", lapply(vector, function(x){
  #   print(x)
  #   dd[x,]
  # }))
  #dd
}
make_ase_sce <- function(sce_full, data_reference, data_alternative){
  colnames(sce_full) <- sce_full$Barcode
  measured.genes = union(rownames(data_reference), 
                         rownames(data_alternative))
  measured.genes = union(rownames(sce_full), measured.genes)
  data_full = add_empty_rows(as.matrix(counts(sce_full)), measured.genes)
  data_reference = add_empty_rows(data_reference, measured.genes)
  data_alternative = add_empty_rows(data_alternative, measured.genes)
  
  #return(list(data_full, data_reference, data_alternative))
  
  barcodes_use <- sce_full$Barcode
  barcodes_use <- intersect(barcodes_use, 
                            intersect(colnames(data_reference), 
                                      colnames(data_alternative)))
  
  sce <- SingleCellExperiment(
    assays = list("counts" = data_full,
                  "counts_reference" = data_reference[,barcodes_use],
                  "counts_alternative" = data_alternative[,barcodes_use])
  )
  sce
}

# Sample1
sce.sample1 <- make_ase_sce(sce.all[ ,sce.all$Library == "Sample1"], 
                            data.sample1$reference,
                            data.sample1$alternative)
sce.sample1$Library <- "Sample1"
sce.sample1$Species <- "B6"
rm(data.sample1)

reducedDims(sce.sample1)$PCA <- reducedDims(sce.all)$PCA
reducedDims(sce.sample1)$TSNE <- reducedDims(sce.all)$TSNE

assays(sce.sample1)[["counts"]] <- as(counts(sce.sample1), "dgTMatrix")
assays(sce.sample1)[["counts_reference"]] <- as(assays(sce.sample1)[["counts_reference"]], "dgTMatrix")
assays(sce.sample1)[["counts_alternative"]] <- as(assays(sce.sample1)[["counts_alternative"]], "dgTMatrix")

saveRDS(sce.sample1, "./ProcessedData/NPC_10x.rds")

###

annotate_chromosome_ens <- function(to_annotate, by_identifier = "gene_symbol"){
  library(EnsDb.Mmusculus.v79)
  all_genes <- genes(EnsDb.Mmusculus.v79)
  if ( by_identifier == "gene_symbol" ){
    df_here <- data.frame(
      chromosome_name = seqnames(all_genes),
      genes = all_genes$symbol
    )
    df_here <- df_here[!duplicated(df_here$genes), ]
    rownames(df_here) <- df_here$genes
  }
  df_here[to_annotate, ]$chromosome_name
}

library("EnsDb.Mmusculus.v79")

new_rowdata <- data.frame(
  row.names = rownames(sce.sample1), 
  gene_ids = rownames(sce.sample1), 
  symbol = rownames(sce.sample1)
) %>%
  add_column("chromosome" = annotate_chromosome_ens(.$symbol))

rowData(sce.sample1) <- new_rowdata

sce.sample1 <- sce.sample1[!is.na(rowData(sce.sample1)$chromosome), ]

## check that indeed, all cells have the same Xi status: 
cell_stats <- data.frame(
  b6_counts_total = colSums(counts_reference(sce.sample1)), 
  cast_counts_total = colSums(counts_alternative(sce.sample1)), 
  b6_counts_X =  colSums(counts_reference(sce.sample1[rowData(sce.sample1)$chromosome == "X" & rownames(sce.sample1) != "Xist", ])), 
  cast_counts_X =  colSums(counts_alternative(sce.sample1[rowData(sce.sample1)$chromosome == "X" & rownames(sce.sample1) != "Xist", ])), 
  b6_counts_Xist =  colSums(counts_reference(sce.sample1[rownames(sce.sample1) == "Xist", ])), 
  cast_counts_Xist =  colSums(counts_alternative(sce.sample1[rownames(sce.sample1) == "Xist", ]))
) %>% 
  mutate(ratio_total = cast_counts_total / (b6_counts_total + cast_counts_total)) %>%
  mutate(ratio_X = cast_counts_X / (b6_counts_X + cast_counts_X)) %>%
  mutate(ratio_Xist = cast_counts_Xist / (b6_counts_Xist + cast_counts_Xist))

# look at distribution of X and Xist d-scores
cell_stats %>%
  ggplot(aes(x = ratio_X)) + geom_density() + xlim(0, 1) + theme_bw()
cell_stats %>%
  ggplot(aes(x = ratio_Xist)) + geom_density() + xlim(0, 1) + theme_bw()
cell_stats %>%
  ggplot(aes(x = ratio_X, y = ratio_Xist)) + geom_point() + xlim(0, 1) + ylim(0, 1) + theme_bw()

gene_stats <- data.frame(
  gene = rownames(sce.sample1), 
  chromosome = rowData(sce.sample1)$chromosome, 
  b6_per_gene = rowSums(counts_reference(sce.sample1)), 
  cast_per_gene = rowSums(counts_alternative(sce.sample1))
) %>%
  mutate(total = b6_per_gene + cast_per_gene) %>%
  mutate(ase = cast_per_gene / total) # cast is inactive, so in the denominator here

gene_stats %>%
  dplyr::filter(total > 500) %>%
  dplyr::filter(!is.na(chromosome)) %>%
  mutate(chromosome_color = ifelse(!(chromosome %in% c("MT", "X")), "Autosomes", chromosome)) %>% 
  ggplot(aes(x = total, y = ase, col = chromosome_color)) + geom_point() + scale_x_log10() + 
  theme_bw() # looks like mom in the cross was B6, no reads on cast MT

gene_stats %>%
  dplyr::filter(total > 500) %>%
  arrange(-ase) %>%
  head(n = 50)

gene_stats %>%
  #dplyr::filter(total > 100) %>%
  dplyr::filter(chromosome == "X") %>%
  ggplot(aes(x = total, y = ase)) + geom_point() + scale_x_log10() + 
    ggrepel::geom_text_repel(aes(label = gene)) + theme_bw()

# plot individual gene
escapees <- gene_stats %>%
  dplyr::filter(chromosome == "X") %>%
  dplyr::filter(total > 300) %>%
  dplyr::filter(ase < 0.8 & ase > 0.1)

### Now look at variability of individual escapees
ratios = as.matrix(counts_alternative(sce.sample1[rownames(escapees), ]) / 
              ( counts_reference(sce.sample1[rownames(escapees), ]) + counts_alternative(sce.sample1[rownames(escapees), ])  ))

ratios_nona <- ratios
ratios_nona[is.na(ratios_nona)] <- 0

correlations_ratios = cor(t(ratios_nona))
correlations_ratios[correlations_ratios > 0.999] <- 0

correlations_ratios %>%
  data.frame() %>%
  rownames_to_column("gene1") %>%
  pivot_longer(-c("gene1")) %>%
  arrange(-value)

gene_stats %>%
  dplyr::filter(total > 500) %>% rownames() -> genes_here
genes_here <- intersect(genes_here, rownames(sce.sample1))

counts_inactive <- counts_alternative(sce.sample1[genes_here, ])
counts_active <- counts_reference(sce.sample1[genes_here, ])

# all_fits <- do.call("rbind", lapply(rownames(counts_active), function(i){
#   print(i)
#   tryCatch({
#     data_in <- data.frame(y = as.numeric(counts_inactive[i, ]), N = as.numeric(counts_active[i, ] + counts_inactive[i, ]))
#     data_in <- data_in[rowSums(data_in) > 0, ]
#     fit <- vglm(cbind(y, N - y) ~ 1, betabinomial, data = data_in)
#     coefs <- Coef(fit)
#     mean = coefs[[1]]
#     disp = coefs[[2]]
#     c(mean, disp)
#   }, 
#   error = function(cond){
#     return(c(NA, NA))
#   })
# }))
# 
# #### 
# 
# all_fits <- data.frame(all_fits)
# all_fits$symbol <- rownames(counts_active)
# colnames(all_fits) <- c("mean", "var", "symbol")
# #all_fits$theta_new <- 1 / (1 + all_fits$var)
# #all_fits$mu <- all_fits$X1 / (all_fits$X1 + all_fits$X2)
# #all_fits$theta <- 1 / (all_fits$X1 + all_fits$X2 + 1)
# 
# hist(all_fits$mean, breaks = 100)
# 
# saveRDS(all_fits, "~/Desktop/Projects/XChromosome_Project/ProcessedData/NPC_variability_model_coefficients.rds")
all_fits <- readRDS("~/Desktop/Projects/XChromosome_Project/ProcessedData/NPC_variability_model_coefficients.rds")

# all_fits <- all_fits %>%
#   column_to_rownames("symbol")

#gene_stats_here <- gene_stats[gene_stats$counts_active + gene_stats$counts_inactive > 100, ]

#gene_stats_here$ml_mu <- all_fits$mean
#gene_stats_here$ml_theta <- all_fits$var

all_fits <- all_fits %>%
  add_column(total = rowSums(counts_active) + rowSums(counts_inactive)) %>%
  add_column(empirical_mu = rowSums(counts_inactive) / .$total) %>%
  add_column(chromosome = rowData(sce.sample1[all_fits$symbol, ])$chromosome)

all_fits %>%
  ggplot(aes(x = total, y = empirical_mu)) + geom_point() + scale_x_log10()

all_fits %>%
  ggplot(aes(x = total, y = mean)) + geom_point() + scale_x_log10()

all_fits %>%
  ggplot(aes(x = empirical_mu, y = mean)) + geom_point()

all_fits %>%
  ggplot(aes(x = mean, y = var)) + geom_point()

# remove the extremely imbalanced genes

all_fits %>% data.frame() %>%
  dplyr::filter(!(mean < 0.1 | mean > 0.9)) %>%
  ggplot(aes(x = mean, y = var)) + geom_point()

all_fits_here <- all_fits %>% data.frame() %>%
  dplyr::filter(!(mean < 0.1 | mean > 0.9)) 

all_fits_here %>%
  arrange(-var) %>%
  head(n = 50)

all_fits_here %>%
  ggplot(aes(x = chromosome == "X", y = var)) +  geom_jitter() + geom_boxplot(outlier.color = NA, width = 0.3)  + theme_bw()

all_fits_here %>%
  ggplot(aes(x = mean, var)) + geom_point() + geom_smooth() + theme_bw()

all_fits_here %>%
  ggplot(aes(x = total, var)) + 
  geom_point(data = all_fits_here[all_fits_here$chromosome != "X", ], col = "grey") + 
  geom_point(data = all_fits_here[all_fits_here$chromosome == "X", ], col = "red") + 
  ggrepel::geom_text_repel(data = all_fits_here[all_fits_here$chromosome == "X", ], col = "red", aes(label = symbol)) + 
  geom_smooth() + scale_x_log10() + theme_bw()

all_fits_here %>%
  ggplot(aes(x = mean, var)) + 
  geom_point(data = all_fits_here[all_fits_here$chromosome != "X", ], col = "grey") + 
  geom_point(data = all_fits_here[all_fits_here$chromosome == "X", ], col = "red") + 
  ggrepel::geom_text_repel(data = all_fits_here[all_fits_here$chromosome == "X", ], col = "red", aes(label = symbol)) + 
  geom_smooth() + xlim(0, 1) + theme_bw()

all_fits_here %>%
  dplyr::filter(chromosome == "X") %>%
  dplyr::arrange(var) %>%
  head()

all_fits_here %>%
  dplyr::arrange(-var) %>%
  dplyr::filter(total > 2000) %>%
  dplyr::filter(chromosome == "X") %>%
  head(n = 30)

gene_here <- "Gfod2"
data.frame("total" = counts_active[gene_here, ] + counts_inactive[gene_here, ], 
           "ase" = counts_inactive[gene_here, ] / (counts_active[gene_here, ] + counts_inactive[gene_here, ])) %>%
  ggplot(aes(x = total, y = ase)) + geom_jitter() + stat_summary(aes(x = 10), col = "red") + scale_x_log10() + ggtitle(gene_here)

## Check that the inferred parameters make sense
library(VGAM)
data_in <- data.frame(y = sample(c(0, 1), 100, replace = T) * 5, N = 1 * 5)
data_in <- data.frame(y = rbinom(n = 100, size = 5, prob = 0.5), N = 1 * 5)
fit <- vglm(cbind(y, N - y) ~ 1, betabinomial, data = data_in)
coefs <- Coef(fit)
mean = coefs[[1]]
disp = coefs[[2]]
c(mean, disp)

# Alternative analysis: Which genes have higher than expected autosomal overdispersion? RAME?
all_fits_here %>%
  ggplot(aes(x = total, var)) + 
    geom_point(data = all_fits_here[all_fits_here$chromosome != "X", ], col = "grey") + 
    geom_point(data = all_fits_here[all_fits_here$chromosome == "X", ], col = "red") + 
    ggrepel::geom_text_repel(data = all_fits_here[all_fits_here$chromosome == "X", ], col = "red", aes(label = symbol)) + 
    geom_smooth() + scale_x_log10() + theme_bw()

all_fits_here %>%
  ggplot(aes(x = mean, var)) + 
    geom_point(data = all_fits_here[all_fits_here$chromosome != "X", ], col = "grey") + 
    geom_point(data = all_fits_here[all_fits_here$chromosome == "X", ], col = "red") + 
    ggrepel::geom_text_repel(data = all_fits_here[all_fits_here$chromosome == "X", ], col = "red", aes(label = symbol)) + 
    geom_smooth() + xlim(0, 1) + theme_bw()

# We first filter out genes with < 100 reads
#all_fits_here_here <- all_fits_here %>%
#  dplyr::filter(total > 1000)

# Fit overdispersion as a loess of log(expression) and abs(mean - 0.5): 
loess_fit <- loess(data = all_fits_here, formula = var ~ log(total) + abs(mean - 0.5))

plot(all_fits_here$var, predict(loess_fit))
cor(all_fits_here$var, predict(loess_fit))

all_fits_here$residuals <- loess_fit$residuals

all_fits_here <- all_fits_here %>% dplyr::arrange(-residuals) %>%
  mutate(highlight = residuals > 0.2)

all_fits_here %>%
  ggplot(aes(x = total, var)) + 
  geom_point() + 
  ggrepel::geom_text_repel(data = all_fits_here[all_fits_here$highlight, ], col = "red", aes(label = symbol)) + 
  geom_smooth() + scale_x_log10() + theme_bw()

gene_here <- "Nrg1"
data.frame("total" = counts_active[gene_here, ] + counts_inactive[gene_here, ], 
           "ase" = counts_inactive[gene_here, ] / (counts_active[gene_here, ] + counts_inactive[gene_here, ])) %>%
  ggplot(aes(x = total, y = ase)) + geom_jitter() + stat_summary(aes(x = 10), col = "red") + scale_x_log10() + ggtitle(gene_here)

rame_genes <- c("Eya3", "App", "Ptk2b", "Snca", "Cnrip1", "Eya1", "Eya2", "Kcnq2", "A2m", "Acyp2", "Bag3", "Eya4", "Grik2")

all_fits_here <- all_fits_here %>% 
  mutate(is_rame = symbol %in% rame_genes)

all_fits_here %>%
  ggplot(aes(x = total, var)) + 
  geom_point() + 
  ggrepel::geom_text_repel(data = all_fits_here[all_fits_here$is_rame, ], col = "red", aes(label = symbol)) + 
  geom_smooth() + scale_x_log10() + theme_bw()

### compare the escapee profile to the bulk RNA-Seq data (which clone is it?)

bulk_data <- readRDS("~/Desktop/Projects/XChromosome_Antonia/ProcessedData/merged_dataset.rds")
bulk_data_X <- bulk_data[seqnames(rowRanges(bulk_data)) == "X", ]
bulk_data_X <- bulk_data_X[,is.na(bulk_data_X$Condition) | (bulk_data_X$Condition == "Aux_0_Dox_0_WO_NO_WOAuxNO")]
bulk_data_X <- bulk_data_X[rowMeans((assays(bulk_data_X)[["counts_active"]] + assays(bulk_data_X)[["counts_inactive"]])) > 50, ]
colData(bulk_data)

ratios_bulk <- assays(bulk_data_X)[["counts_active"]] / (assays(bulk_data_X)[["counts_active"]] + assays(bulk_data_X)[["counts_inactive"]])

gene_stats_X <- gene_stats %>%
  dplyr::filter(chromosome == "X") %>%
  dplyr::filter(total > 50) %>%
  dplyr::select(c("gene", "ase")) %>%
  add_column("sample" = "bulk") %>%
  tibble() %>%
  relocate(c("gene", "sample", "ase"))

### merge 
ratios_bulk_merged <- ratios_bulk %>%
  rownames_to_column("gene") %>%
  pivot_longer(-c("gene")) %>%
  rename("ase" = value, "sample" = name) %>%
  rbind(., gene_stats_X)

ratios_bulk_merged %>%
  pivot_wider(values_from = ase, names_from = sample) %>%
  tidyr::gather(value, variable, -c(bulk, gene)) %>%
  mutate(is_na_bulk = is.na(bulk)) %>%
  mutate(is_na_samples = is.na(variable)) %>%
  mutate(is_na_color = paste0(is_na_bulk, "_", is_na_samples)) %>%
  mutate(bulk = replace_na(bulk, 1)) %>%
  mutate(variable = replace_na(variable, 0)) %>%
  ggplot(aes(x = bulk, y = variable, col = is_na_color)) + geom_point() + facet_wrap(~value) + theme_bw()

ratios_bulk_merged %>%
  dplyr::filter(sample %in% c("RNA_NPC_C5_C57BL.6J_rep1", "RNA_NPC_C5_C57BL.6J_rep2", "bulk")) %>%
  pivot_wider(values_from = ase, names_from = sample) %>%
  tidyr::gather(value, variable, -c(bulk, gene)) %>%
  mutate(is_na_bulk = is.na(bulk)) %>%
  mutate(is_na_samples = is.na(variable)) %>%
  mutate(is_na_color = paste0(is_na_bulk, "_", is_na_samples)) %>%
  mutate(bulk = replace_na(bulk, 1)) %>%
  mutate(variable = replace_na(variable, 0)) %>%
  ggplot(aes(y = bulk, x = 1 - variable, col = is_na_color)) + geom_point() + facet_wrap(~value) + 
  geom_abline(linetype = 'dashed') + theme_paper() + 
    xlab("Single-cell RNA-Seq") + ylab("Bulk RNA-Seq") + theme_bw()

ratios_bulk_merged %>%
  dplyr::filter(sample %in% c("RNA_NPC_C5_C57BL.6J_rep1", "RNA_NPC_C5_C57BL.6J_rep2", "bulk")) %>%
  pivot_wider(values_from = ase, names_from = sample) %>%
  tidyr::gather(value, variable, -c(bulk, gene)) %>%
  mutate(is_na_bulk = is.na(bulk)) %>%
  mutate(is_na_samples = is.na(variable)) %>%
  mutate(is_na_color = paste0(is_na_bulk, "_", is_na_samples)) -> test

# It looks like most of the genes are either < 0.1, > 0.6 or pseudogenes
test %>%
  dplyr::filter(is_na_bulk) %>%
  dplyr::filter(value != "RNA_NPC_C5_C57BL.6J_rep1") %>%
  View() # 10 genes between 0.1 / 0.7, most pseudogenes

test %>%
  dplyr::filter(is_na_samples) %>%
  dplyr::filter(value != "RNA_NPC_C5_C57BL.6J_rep1") %>%
  View() #17 genes between 0.1 / 0.7, most pseudogenes

## look at cell cycle dependencies of genes

## ---------------------- dali functions ---------------------- 

# initialize python path
library(reticulate)
use_python("/Users/jasper/Library/r-miniconda/bin/python", required = T)
dali_path = "~/Desktop/PhD/Projects/ASE_Spermatogenesis/Scripts/dali/"

# define R wrappers around scDALI functions
test_regions_R <- function(A, D, 
                           cell_state, 
                           n_cores = 1L)
{
  source_python(paste0(dali_path, "/dali/my_test_regions.py"))
  test_regions(np_array(A), 
               np_array(D), 
               np_array(cell_state), 
               n_cores = n_cores)
}

test_mean_R <- function(a, d, mean_mu = 0.5)
{
  source_python(paste0(dali_path, "/dali/my_test_mean.py"))
  res = test_mean(np_array(a), 
                  np_array(d), 
                  mean_mu = mean_mu)
  names(res) = c("Mean", "Theta", "NLL_null", "NLL_alt", "pval")
  res
}

run_gp_R <- function(A, D, 
                     cell_state, 
                     kernel = "Linear",
                     n_cores = 1L)
{
  source_python(paste0(dali_path, "/dali/my_run_gp.py"))
  res = run_gp(np_array(A), 
               np_array(D), 
               np_array(cell_state), 
               kernel = kernel,
               n_cores = n_cores)
  res
}

genes_test <- gene_stats %>%
  dplyr::filter(chromosome == "X") %>%
  dplyr::filter(b6_per_gene + cast_per_gene > 100) %>%
  dplyr::filter(ase > 0.1 & ase < 0.7)

test_pca <- calculatePCA(sce.all[rowData(sce.sample1)$chromosome != "X", ])

A <- as.matrix(t(counts_reference(sce.sample1[rownames(genes_test),])))
D <- as.matrix(t(counts_reference(sce.sample1[rownames(genes_test),]) + 
         counts_alternative(sce.sample1[rownames(genes_test),])))
pseudotime <- as.numeric(as.factor(sce.all$cycle_stage))
pseudotime <- as.matrix(test_pca[,1:10])

results_linear_kernel <- test_regions_R(A, D, pseudotime, n_cores = 4L)
names(results_linear_kernel) <- rownames(genes_test)
sort(results_linear_kernel)
sort(p.adjust(results_linear_kernel))
sort(p.adjust(results_linear_kernel)[p.adjust(results_linear_kernel) < 0.1])
tt <- sort(p.adjust(results_linear_kernel))

ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  cluster = as.numeric(logcounts(sce.all)["Gpm6b", ]))) +
  geom_point(aes(tsne1, tsne2, colour = cluster)) + scale_color_viridis()

rotations <- attr(test_pca, "rotation")
sort(rotations[,1], decreasing = F)

genes_mouse <- genes(EnsDb.Mmusculus.v79)

genes_show <- names(tt[tt < 0.1])
genes_show_reference <- names(tt)
genes_show_total <- c(genes_show, genes_show_reference)

pos <- position_jitter(height = 0.2, width = 0,  seed = 2)
data.frame(
  starts = start(genes_mouse[match(genes_show_total, genes_mouse$symbol), ]), 
  genes = genes_show_total, 
  cat = c(rep("variable", length(genes_show)), rep( "reference", length(genes_show_reference)))
) %>%
  ggplot(aes(x = starts, y = cat)) + geom_jitter(position = pos) + ggrepel::geom_text_repel(aes(label = genes), position = pos)

# run test against every PC to see what they associate with

all_results <- lapply(1:10, function(i){
  
  pseudotime <- as.matrix(test_pca[,i])
  results_linear_kernel <- test_regions_R(A, D, pseudotime, n_cores = 4L)
  names(results_linear_kernel) <- rownames(genes_test)
  
  results_linear_kernel

})

all_results_df <- lapply(1:10, function(i){
  data.frame(
    PC = paste0("PC_", i),
    genes = names(all_results[[i]]),
    pvals = as.numeric(all_results[[i]])
  )
}) %>% do.call("rbind", .)

all_results_df %>%
  dplyr::filter(genes %in% names(tt[tt < 0.1])) %>%
  mutate(padj = p.adjust(pvals)) %>%
  add_column(position = start(genes_mouse[match(.$genes, genes_mouse$symbol), ])) %>%
  ggplot(aes(x = PC, y = reorder(genes, position), fill = padj < 0.1)) + geom_tile(col = "black") + theme_classic()

# 

rotations_df <- lapply(1:10, function(i){
  data.frame(
    Gene = rownames(rotations), 
    Index = 1:nrow(rotations), 
    PC = paste0("PC_", i), 
    loading = rotations[,i]
  )
}) %>% do.call("rbind", .)

rotations_df %>%
  group_by(PC) %>%
  slice_max(order_by = abs(loading), n = 10) %>%
  mutate(Index = 1:10) %>%
  ungroup() %>%
  ggplot(aes(x = reorder(Index, loading), y = loading)) + geom_point() + facet_wrap(~PC) + 
    ggrepel::geom_text_repel(aes(label = Gene))

# 

gene = "Gpm6b"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  cluster = as.numeric(ratios[gene, ]))) +
  geom_point(aes(tsne1, tsne2, colour = cluster)) + scale_color_gradient2(low = "blue", high = "red", midpoint = 0.5)

data.frame(PC5 = reducedDims(sce.all)$PCA[,5],
           PC6 = reducedDims(sce.all)$PCA[,6],
           ase = as.numeric(ratios[gene, ])) %>% ggplot(aes(x = PC5, y = PC6, col = ase)) + geom_point() + 
            scale_color_gradient2(low = "blue", high = "red", midpoint = 0.5) + theme_bw()

data.frame(tsne1 = reducedDims(sce.all)$PCA[,5],
           tsne2 = reducedDims(sce.all)$PCA[,6],
           cluster = as.numeric(ratios[gene, ])) %>% ggplot(aes(x = tsne1, y = cluster)) + geom_point() + geom_smooth()


data.frame(tsne1 = reducedDims(sce.all)$PCA[,5],
           tsne2 = reducedDims(sce.all)$PCA[,6],
           cluster = as.numeric(logcounts(sce.all)["Hmga2", ])) %>% ggplot(aes(x = tsne1, y = cluster)) + geom_point() + geom_smooth()

data.frame(tsne1 = reducedDims(sce.all)$PCA[,5],
           tsne2 = reducedDims(sce.all)$PCA[,6],
           cluster = sce.all$cycle_stage) %>% ggplot(aes(x = tsne1, y = tsne2, col = cluster)) + geom_point()

pval_cutoff <- data.frame(pval = tt, padj = p.adjust(tt)) %>%
  arrange(pval) %>% dplyr::filter(padj < 0.1)

data.frame(gene = names(tt), "pval" = tt) %>%
  ggplot(aes(x = reorder(gene, pval), y = -log10(pval), col = p.adjust(pval) < 0.1)) + geom_point() + coord_flip() + theme_bw() + 
  theme(text = element_text(size = 8))

data.frame(gene = names(tt), "pval" = tt) %>%
  dplyr::filter(pval < 0.2) %>%
  ggplot(aes(x = reorder(gene, pval), y = -log10(pval), col = p.adjust(pval) < 0.1)) + geom_point() + coord_flip() + theme_bw() + 
  theme(text = element_text(size = 8)) + ggrepel::geom_text_repel(aes(label = gene))

### Now look at variability of individual escapees
ratios = (as.matrix(counts_alternative(sce.sample1[rownames(escapees), ]))) / 
                     ( as.matrix( counts_reference(sce.sample1[rownames(escapees), ]) + counts_alternative(sce.sample1[rownames(escapees), ])  ) )

#ratios[as.matrix( counts_reference(sce.sample1[rownames(escapees), ]) + counts_alternative(sce.sample1[rownames(escapees), ])  ) < 2] <- NA

ratios_nona <- ratios
# ratios_nona[is.na(ratios_nona)] <- 0

test_numbers <- apply(ratios_nona, 1, function(x){unlist(apply(ratios_nona, 1, function(y){sum(!is.na(x) & !is.na(y))}))})
test <- apply(ratios_nona, 1, function(x){unlist(apply(ratios_nona, 1, function(y){cor(x, y, use = "complete.obs")}))})

# look at how number of obs determines correlations: 

df_here <- data.frame(
  number_of_nonas = as.numeric(test_numbers),
  correlations = as.numeric(test)
)
df_here[is.na(df_here$correlations), ]$correlations <- 0

df_here %>%
  ggplot(aes(x = number_of_nonas, y = abs(correlations))) + geom_point() + theme_bw() + scale_x_log10()

# loess fit and extract residuals

loess_here <- loess((correlations) ~ (number_of_nonas), df_here)
df_here$residuals <-  loess_here$residuals

df_here %>%
  ggplot(aes(x = number_of_nonas, y = residuals)) + geom_point() + theme_bw() + scale_x_log10()

df_here$Gene1 <- rep(rownames(ratios_nona), nrow(ratios_nona))
df_here$Gene2 <- unlist(lapply(rownames(ratios_nona), function(x){rep(x, nrow(ratios_nona))}))

df_here %>%
  dplyr::filter(number_of_nonas > 100) %>%
  dplyr::arrange(-residuals) %>%
  head(n = 20)

pheatmap(test)

data.frame(
  x = as.numeric(ratios_nona["Trappc2", ]), 
  y = as.numeric(ratios_nona["Gpm6b", ])
) %>%
  ggplot(aes(x = x, y = y)) + geom_jitter() + theme_bw() + geom_smooth(method = "lm")





library(MOFA2)
# 
# # get variable genes in expression space
expression_data <- logcounts(sce.all[HVgenes], )

MOFAobject <- create_mofa(data = list("View1" = as.matrix(expression_data)))
MOFAobject <- prepare_mofa(MOFAobject)
MOFAobject <- run_mofa(MOFAobject)

plot_variance_explained(MOFAobject)
plot_weights(MOFAobject, factors = 1:5, view = 1)

plot_factor(MOFAobject, factors = 1:10)
MOFAobject <- MOFA2::run_tsne(MOFAobject)
MOFA2::plot_dimred(MOFAobject)

all_factors <- get_factors(MOFAobject)$group1
all_weights <- get_weights(MOFAobject)$View1

head(all_weights[order((all_weights[,10]), decreasing = F), ], n = 20)

tsne_mofa <- Rtsne::Rtsne(all_factors)

data.frame(
  Dim1 = tsne_mofa$Y[,1], 
  Dim2 = tsne_mofa$Y[,2], 
  Factor = all_factors[,10]
) %>% ggplot(aes(Dim1, Dim2, col = Factor)) + geom_point() + theme_bw() + scale_color_viridis()



all_results_mofa <- lapply(1:10, function(i){
  
  pseudotime <- as.matrix(all_factors[,i])
  results_linear_kernel <- test_regions_R(A, D, pseudotime, n_cores = 4L)
  names(results_linear_kernel) <- rownames(genes_test)
  
  results_linear_kernel
  
})

all_results_df <- lapply(1:10, function(i){
  data.frame(
    Dim = paste0("Dim_", i),
    genes = names(all_results_mofa[[i]]),
    pvals = as.numeric(all_results_mofa[[i]])
  )
}) %>% do.call("rbind", .)

all_results_df %>%
  dplyr::filter(genes %in% names(tt[tt < 0.1])) %>%
  mutate(padj = p.adjust(pvals)) %>%
  add_column(position = start(genes_mouse[match(.$genes, genes_mouse$symbol), ])) %>%
  ggplot(aes(x = Dim, y = reorder(genes, position), fill = padj < 0.1)) + geom_tile(col = "black") + theme_classic()

# 

