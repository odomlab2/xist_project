# exploratory analysis of antonias data

library(ggplot2)
library(tidyverse)

source("~/Desktop/PhD/RCode/Multipurpose/auxiliary.R")
source("~/Desktop/Projects/XChromosome_Project/Scripts/functions.R")

data <- read.csv("~/Desktop/Projects/XChromosome_Antonia/Data/December_2021/table_raw_counts.txt", sep = "\t")

# isolate genes information

genes <- data[,1:2]
data <- data[,-c(1:2)]

# add chromosome information to genes
# ...
# first add gene symbols
library("EnsDb.Mmusculus.v79")
symbols <- mapIds(EnsDb.Mmusculus.v79, keys = genes$Geneid, keytype = "GENEID", column="SYMBOL")
genes$symbol <- symbols
# genes <- genes[!is.na(genes$symbol), ] # ???? 10k dont map to gene symbols

gene_info <- data.frame(genes(EnsDb.Mmusculus.v79))
rownames(gene_info) <- gene_info$gene_id
genes$chromosome <- gene_info[genes$Geneid, ]$seqnames
genes <- genes[,colnames(genes) != "Length"]

# read in metadata

metadata <- readxl::read_excel("~/Desktop/Projects/XChromosome_Antonia/Data/December_2021/RNAseq-samples-CL30_31.xlsx")
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Allele", "DayOfExperiment", "ndAux", "timeAux", "ndDox", "timeDox", 
                        "washout", "timeWO", "WOwithAux", "Replicates")
metadata <- metadata %>% 
  mutate(SampleName = gsub(('\"'), "", gsub("^[^CL]*", "", SampleName))) %>%
  add_column("Condition" = paste0("Aux_", .$ndAux, "_Dox_", .$ndDox, "_WO_", .$washout, "_WOAux", .$WOwithAux))

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "allelic reads + unassigned reads")

metadata_b6 <- metadata %>%
  dplyr::filter(Allele == "C57BL.6J")

metadata_cast <- metadata %>%
  dplyr::filter(Allele == "CAST.EiJ")

# reorder samples
data <- data[,metadata$SampleName]

# analysis
# First look at total counts
data.all <- data[, metadata$Allele == "allelic reads + unassigned reads"]
data.b6 <- data[, metadata$Allele == "C57BL.6J"]
data.cast <- data[, metadata$Allele == "CAST.EiJ"]

data.frame(counts = colSums(data.all)) %>%
  rownames_to_column("Sample") %>%
  ggplot(aes(x = Sample, y = counts)) + geom_bar(stat = "identity") + 
  theme_paper() + coord_flip() + 
  ggtitle("Total Counts (All)")

data.frame(counts = colSums(data.b6)) %>%
  rownames_to_column("Sample") %>%
  ggplot(aes(x = Sample, y = counts)) + geom_bar(stat = "identity") + 
  theme_paper() + coord_flip() + 
  ggtitle("Total Counts (B6)")

data.frame(counts = colSums(data.cast)) %>%
  rownames_to_column("Sample") %>%
  ggplot(aes(x = Sample, y = counts)) + geom_bar(stat = "identity") + 
  theme_paper() + coord_flip() + 
  ggtitle("Total Counts (Cast)")

data.frame(ratio =  colSums(data.b6) / ( colSums(data.b6) + colSums(data.cast)) ) %>%
  rownames_to_column("Sample") %>%
  ggplot(aes(x = Sample, y = ratio)) + geom_bar(stat = "identity") + 
  theme_paper() + coord_flip() + 
  ggtitle("Allelic ratios (Autosomes + X)") + ylim(0, 1) + 
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "red")

# look at counts per gene
data.frame(counts = rowSums(data.all)) %>%
  ggplot(aes(x = counts)) + geom_histogram(bins = 100) + 
  scale_x_log10() + theme_paper()

data.frame(
  B6_Counts = rowSums(data.b6), 
  CAST_Counts = rowSums(data.cast), 
  chromosome = genes$chromosome
) %>%
  ggplot(aes(x = B6_Counts + 1, y = CAST_Counts + 1)) + geom_point() + 
  theme_paper() + scale_x_log10() + scale_y_log10()

data.frame(
  B6_Counts = rowSums(data.b6), 
  CAST_Counts = rowSums(data.cast), 
  chromosome = genes$chromosome
) %>%
  ggplot(aes(x = B6_Counts + 1, y = CAST_Counts + 1, col = chromosome == "X")) + geom_point() + 
  theme_paper() + scale_x_log10() + scale_y_log10()

# Look at autosomal gene expression patterns
library(DESeq2)

des <- DESeqDataSetFromMatrix(data.all[,metadata_reduced$SampleName], colData = data.frame(metadata_reduced %>% column_to_rownames("SampleName")), design = ~ 1)
des.vst <- vst(des)
rowData(des.vst)$symbol <- genes$symbol

var_data <- rowData(des.vst)
data.frame(head(var_data[order(var_data$baseVar, decreasing = T), ], n = 50))

plotPCA(des.vst, intgroup = "Condition")
# 

genes_x <- genes[genes$chromosome == "X", ]
genes_x <- genes_x[!is.na(genes_x$chromosome), ]

rownames(data.b6) <- genes$Geneid
rownames(data.cast) <- genes$Geneid

data.b6.X <- data.b6[genes_x$Geneid, ]
rownames(data.b6.X) <- make.unique(genes_x$symbol)
data.cast.X <- data.cast[genes_x$Geneid, ]
rownames(data.cast.X) <- make.unique(genes_x$symbol)

genes_include <- rowSums(data.b6.X) > 50 | rowSums(data.cast.X) > 50
table(genes_include)

data.b6.X <- data.b6.X[genes_include, ]
data.cast.X <- data.cast.X[genes_include, ]

data.ratios <- data.b6.X / (data.b6.X + data.cast.X)
data.ratios[is.na(data.ratios)] <- 0

pca_ratios <- prcomp(t(data.ratios))

data.frame(
  PC1 = pca_ratios$x[,1], 
  PC2 = pca_ratios$x[,2], 
  Covariate = metadata_b6$Clone
) %>%
  ggplot(aes(PC1, PC2, col = Covariate)) + geom_point(size = 3)

data.frame(
  PC1 = pca_ratios$x[,1], 
  PC2 = pca_ratios$x[,2], 
  Covariate = metadata_b6$Condition, 
  Clone = metadata_b6$Clone
) %>%
  ggplot(aes(PC1, PC2, col = Covariate, shape = Clone)) + geom_point() + 
  ggrepel::geom_text_repel(aes(label = Covariate))

data.frame(
  PC1 = pca_ratios$rotation[,1], 
  PC2 = pca_ratios$rotation[,2], 
  Name = rownames(pca_ratios$rotation)
) %>%
  ggplot(aes(PC1, PC2)) + geom_point() + ggrepel::geom_text_repel(aes(label = Name))

# try the same thing with genes that show any kind of escape
data.ratios.with.escape <- data.ratios[apply(data.ratios, 1, function(x){any(x > 0.1)}), ]
data.ratios.with.escape[is.na(data.ratios.with.escape)] <- 0

pca_ratios_escapees <- prcomp(t(data.ratios.with.escape))

# basically indistinguishable
plot(pca_ratios_escapees$x[,1], pca_ratios$x[,1])
plot(pca_ratios_escapees$x[,2], pca_ratios$x[,2])

# look for genes with differential allelic balance between conditions
# first global

# first look at Xist
data.frame(
  Sample = metadata_reduced$Condition, 
  Expression = 2 ** des.vst@assays@data[[1]][which(genes$symbol == "Xist"), ], 
  Clone = metadata_reduced$Clone
) %>% 
  ggplot(aes(x = Sample, y = Expression, fill = Clone)) + geom_bar(stat = "identity", position = "dodge") + 
  coord_flip() + theme_paper()

# now look at distribution of escapees

data.ratios.with.escape %>% pivot_longer(-c()) %>%
  add_column(Condition = metadata_b6[match(.$name, metadata_b6$SampleName), ]$Condition) %>%
  add_column(Clone = metadata_b6[match(.$name, metadata_b6$SampleName), ]$Clone) %>%
  ggplot(aes(x = Condition, y = value, col = Clone)) + geom_boxplot(outlier.color = NA) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.05), alpha = 0.2) + coord_flip() + theme_paper()

clone1 <- data.ratios.with.escape[,grepl("CL30", colnames(data.ratios.with.escape))]
clone2 <- data.ratios.with.escape[,grepl("CL31", colnames(data.ratios.with.escape))]

results_clone1 <- do.call("rbind", lapply(1:ncol(clone1), function(i){
  median_diff <- mean(clone1[,i]) - mean(clone1[,1])
  pval_wilc <- wilcox.test(clone1[,i], clone1[,1], paired = T)$p.value
  return(c(colnames(clone1)[i], median_diff, pval_wilc))
}))
colnames(results_clone1) <- c("Sample", "mean_difference", "pval")

results_clone2 <- do.call("rbind", lapply(1:ncol(clone2), function(i){
  median_diff <- mean(clone2[,i]) - mean(clone2[,1])
  pval_wilc <- wilcox.test(clone2[,i], clone2[,1], paired = T)$p.value
  return(c(colnames(clone2)[i], median_diff, pval_wilc))
}))
colnames(results_clone2) <- c("Sample", "mean_difference", "pval")

results_to_baseline <- rbind(results_clone1, results_clone2) %>%
  data.frame() %>% mutate(pval = as.numeric(pval), mean_difference = as.numeric(mean_difference)) %>%
  add_column(p_adj = p.adjust(.$pval)) %>%
  add_column(Condition = metadata[match(.$Sample, metadata$SampleName), c("Condition")]$Condition) %>%
  add_column(Clone = metadata[match(.$Sample, metadata$SampleName), c("Clone")]$Clone)

results_to_baseline %>%
  ggplot(aes(Condition, mean_difference, shape = Clone)) + geom_point(aes(size = -log10(pval), col = p_adj < 0.1)) + 
  coord_flip() + geom_hline(yintercept = 0, linetype = "dashed") + theme_paper()

# check how escapee covaries across genomic space
pheatmap::pheatmap(cor(t(data.ratios), use = "everything"), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, 
                   cellheight=1, cellwidth = 1)

# reproduce heatmaps from antonias presentation
testy <- data.ratios
testy <- testy[order(match(rownames(testy), gene_info$gene_name)), ]

# pick samples to show
conds_show <- c("Aux_0_Dox_0_WO_NO_WOAuxNA", "Aux_0_Dox_3_WO_NO_WOAuxNA", "Aux_0_Dox_7_WO_NO_WOAuxNA", 
                "Aux_0_Dox_3_WO_YES_WOAuxNO", "Aux_0_Dox_7_WO_YES_WOAuxYES")
samples_show <- metadata_b6$Clone == "CL30" & metadata_b6$Condition %in% conds_show

testy_cl30 <- testy[,samples_show]
colnames(testy_cl30) <- metadata_b6[samples_show, ]$Condition

# order by chromosome position
testy_cl30 %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene) %>%
  mutate(gene = factor(gene, levels = rownames(testy))) %>%
  ggplot(aes(x = gene, y = factor(name, levels = conds_show), fill = value)) + geom_tile() + theme_paper() + 
  scale_fill_gradient2(low = "#aa7c01", mid = "#fe93a0", high = "#0b6d68", midpoint = 0.5) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 4)) + xlab("")
ggsave("~/Desktop/testy.pdf", width = 50, height = 5, limitsize = F)
ggsave("~/Desktop/testy.pdf", width = 15, height = 5, limitsize = F)

# zoom in 

ints <- c(which(rownames(testy_cl30) == "Cetn2"), which(rownames(testy_cl30) == "Brcc3"))
genes_show <- rownames(testy_cl30)[ints[[1]]:ints[[2]]]

testy_cl30[genes_show, ] %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene) %>%
  mutate(gene = factor(gene, levels = genes_show)) %>%
  ggplot(aes(x = gene, y = factor(name, levels = conds_show), fill = value)) + geom_tile() + theme_paper() + 
  scale_fill_gradient2(low = "#aa7c01", mid = "#fe93a0", high = "#0b6d68", midpoint = 0.5) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 4)) + xlab("")

###### 
# the same for clone 31
######

testy <- data.ratios
testy <- testy[order(match(rownames(testy), gene_info$gene_name)), ]

# pick samples to show
conds_show <- c("Aux_0_Dox_0_WO_NO_WOAuxNA", "Aux_0_Dox_3_WO_NO_WOAuxNA", "Aux_0_Dox_7_WO_NO_WOAuxNA", 
                "Aux_0_Dox_3_WO_YES_WOAuxNO", "Aux_0_Dox_7_WO_YES_WOAuxYES")
samples_show <- metadata_b6$Clone == "CL31" & metadata_b6$Condition %in% conds_show

testy_cl30 <- testy[,samples_show]
colnames(testy_cl30) <- metadata_b6[samples_show, ]$Condition

# order by chromosome position
testy_cl30 %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene) %>%
  mutate(gene = factor(gene, levels = rownames(testy))) %>%
  ggplot(aes(x = gene, y = factor(name, levels = conds_show), fill = value)) + geom_tile() + theme_paper() + 
  scale_fill_gradient2(low = "#aa7c01", mid = "#fe93a0", high = "#0b6d68", midpoint = 0.5) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 4)) + xlab("")

ggsave("~/Desktop/testy.pdf", width = 50, height = 5, limitsize = F)
ggsave("~/Desktop/testy.pdf", width = 15, height = 5, limitsize = F)

# zoom in 
ints <- c(which(rownames(testy_cl30) == "Cetn2"), which(rownames(testy_cl30) == "Brcc3"))
genes_show <- rownames(testy_cl30)[ints[[1]]:ints[[2]]]

testy_cl30[genes_show, ] %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene) %>%
  mutate(gene = factor(gene, levels = genes_show)) %>%
  ggplot(aes(x = gene, y = factor(name, levels = conds_show), fill = value)) + geom_tile() + theme_paper() + 
  scale_fill_gradient2(low = "#aa7c01", mid = "#fe93a0", high = "#0b6d68", midpoint = 0.5) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, size = 4)) + xlab("")

# ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ 
# save stuff... 
# ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ # ------ 

# show 

testy <- cor(t(data.ratios), use = "everything") %>%
  data.frame() %>%
  rownames_to_column("gene_2") %>%
  pivot_longer(-gene_2) %>%
  dplyr::filter(value < 0.99999999)

distance_between_genes <- abs(gene_info[match(testy$gene_2, gene_info$gene_name), ]$start - gene_info[match(testy$name, gene_info$gene_name), ]$start)

data.frame(
  distance = distance_between_genes, 
  correlation = testy$value
) %>%
  ggplot(aes(x = distance, y = correlation)) + geom_point() + 
  geom_smooth() +
  scale_x_log10()

# we might want annotation on which genes are known to be escapees

escapee_annotation <- readxl::read_excel("~/Desktop/Projects/XChromosome_Project/ProcessedData/ListOfEscapeeFromEdithLab.xlsx") %>%
  dplyr::select(c("ENSEMBL_v102", "final status"))

library("EnsDb.Mmusculus.v79") 
symbols <- mapIds(EnsDb.Mmusculus.v79, keys = escapee_annotation$ENSEMBL_v102, keytype = "GENEID", column="SYMBOL")

escapee_annotation$symbol <- symbols

# fit linear models 

# For simplicity exclude the aux late samples

#data.b6.X.here <- data.b6.X[,samples_in_b6]
#data.cast.X.here <- data.cast.X[,samples_in_cast]

metadata_b6_here <- metadata_b6[metadata_b6$WOwithAux == "NO" & metadata_b6$Condition != "Aux_9_Dox_7_WO_NO_WOAuxNO", ]
metadata_b6_here$Condition <- paste0("Aux_", metadata_b6_here$ndAux, "_Dox_", metadata_b6_here$ndDox, "_WO_", metadata_b6_here$washout)

metadata_cast_here <- metadata_cast[metadata_cast$WOwithAux == "NO" & metadata_cast$Condition != "Aux_9_Dox_7_WO_NO_WOAuxNO", ]
metadata_cast_here$Condition <- paste0("Aux_", metadata_cast_here$ndAux, "_Dox_", metadata_cast_here$ndDox, "_WO_", metadata_cast_here$washout)

data.b6.X.here <- data.b6.X[, metadata_b6_here$SampleName]
data.cast.X.here <- data.cast.X[, metadata_cast_here$SampleName]

fitting_metadata_here <- metadata_b6_here[,c("ndDox", "ndAux", "washout", "Clone")]
fitting_metadata_here$ndDox <- factor(fitting_metadata_here$ndDox, levels = unique(fitting_metadata_here$ndDox))
fitting_metadata_here$ndAux <- factor(fitting_metadata_here$ndAux, levels = unique(fitting_metadata_here$ndAux))
fitting_metadata_here$washout <- factor(fitting_metadata_here$washout, levels = unique(fitting_metadata_here$washout))
fitting_metadata_here$Clone <- factor(fitting_metadata_here$Clone, levels = unique(fitting_metadata_here$Clone))

library(VGAM)

# first check which genes we can actually explain any additional variance

data_test_inactive <- data.b6.X.here
data_test_active <- data.cast.X.here

all_coefs <- data.frame(do.call("rbind", lapply(1:nrow(data_test_active), function(i){
  tryCatch({
    y = as.numeric(data_test_inactive[i, ])
    N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
    # fit1 <- glm(cbind(y, N-y) ~ 1, family = "binomial")
    fit1 <- glm(cbind(y, N-y) ~ 1, family = "binomial", data = fitting_metadata_here)
    fit2 <- glm(cbind(y, N-y) ~ 1 + ndDox + ndAux + washout, family = "binomial", data = fitting_metadata_here)
    #kkk <- lrtest(fit2, fit1)
    return(c(logLik(fit1), logLik(fit2)))
  },error=function(cond) {
    return(c(NA))
  })
})))

all_coefs$lrr <- all_coefs$X2 - all_coefs$X1
all_coefs[all_coefs$lrr > 100, ]$lrr <- 100
all_coefs$av_escape <- rowSums(data_test_inactive) / rowSums(data_test_inactive + data_test_active)
all_coefs$av_counts <- rowMeans(data_test_inactive + data_test_active)
all_coefs$gene <- rownames(data_test_active)

all_coefs %>%
  ggplot(aes(av_counts, av_escape)) + geom_point() + scale_x_log10()

all_coefs %>%
  ggplot(aes(av_escape, lrr)) + geom_point()

all_coefs %>%
  ggplot(aes(av_counts, lrr)) + geom_point() + scale_x_log10()

all_coefs$cutoff <- all_coefs$av_escape > 0.05 & all_coefs$av_escape < 0.75
all_coefs$lrr_cutoff <- all_coefs$lrr > 10
all_coefs$combined_use <- all_coefs$av_counts > 10 & all_coefs$av_escape < 0.75 & all_coefs$av_escape > 0.05

genes_test <- all_coefs[all_coefs$lrr_cutoff, ]$gene
genes_test <- all_coefs[all_coefs$cutoff, ]$gene
genes_test <- all_coefs[all_coefs$combined_use, ]$gene

data_test_inactive <- data.b6.X.here[genes_test, ]
data_test_active <- data.cast.X.here[genes_test, ]

all_coefs <- data.frame(do.call("rbind", lapply(1:nrow(data_test_active), function(i){
  tryCatch({
    y = as.numeric(data_test_inactive[i, ])
    N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
    fit1 <- glm(cbind(y, N-y) ~ Clone + ndDox * ndAux * washout, family = "binomial", data = fitting_metadata_here)
    coefs <- coef(fit1)
    coefs
   return(coefs)
  }, error = function(cond) {
    return(c(NA))
  }, warning = function(cond){
    return(c(NA))
  })
})))

rownames(all_coefs) <- rownames(data_test_active)
rownames(all_coefs) <- rownames(data_test_inactive)
all_coefs <- all_coefs[!apply(all_coefs, 1, function(x){all(is.na(x))}), !apply(all_coefs, 2, function(x){all(is.na(x))})]
all_coefs[is.na(all_coefs)] <- 0

# all_coefs %>%
#   rownames_to_column("Gene") %>%
#   pivot_longer(-c(Gene)) %>%
#   ggplot(aes(x = name, y = value)) + geom_boxplot() + 
#   geom_line(aes(group = Gene)) + 
#   coord_flip()

# there are genes with crazy high / low coefficients across the board

all_coefs[abs(all_coefs) > 5] <- sign(all_coefs[abs(all_coefs) > 5]) * 5

row_annotation <- data.frame(
  Gene = rownames(all_coefs), 
  escapee_annotation = escapee_annotation[match(rownames(all_coefs), escapee_annotation$symbol), ][,2]
) %>% column_to_rownames("Gene")

pheatmap::pheatmap(all_coefs, cluster_cols = F, color=colorRampPalette(c("navy", "white", "red"))(50), annotation_row = row_annotation, fontsize = 7)

all_coefs %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-c(Gene)) %>%
  ggplot(aes(x = name, y = value)) + geom_boxplot() +
  coord_flip()

all_coefs %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-c(Gene)) %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ][,2])) %>% 
  # select(c("X.Intercept.", "ndDox3", "ndDox7", "ndDox7.washoutYES", "ndDox7.washoutYES")) %>%
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(x = name, y = value, col = escape_status)) + 
  geom_boxplot() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +
  coord_flip()

intercept_values <- rev_logit(all_coefs$X.Intercept.)
names(intercept_values) <- rownames(all_coefs)
ref_samples <- c("CL30.NoDoxNoAux_C57BL.6J_sorted.bam", "CL31.NoDoxNoAux_C57BL.6J_sorted.bam")
ref_samples_cast <- c("CL30.NoDoxNoAux_CAST.EiJ_sorted.bam", "CL31.NoDoxNoAux_CAST.EiJ_sorted.bam")
mean_references <- rowSums(data_test_inactive[,ref_samples]) / rowSums(data_test_inactive[,ref_samples] + data_test_active[,ref_samples_cast])

plot(mean_references[names(intercept_values)], intercept_values)

head(all_coefs[order(all_coefs$ndDox7, decreasing = F), ], n = 20)

# look at dox treatment vs washout
all_coefs[,c("ndDox7", "ndDox7.washoutYES")] %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(x = ndDox7, y = ndDox7.washoutYES)) + 
  geom_smooth(method = "lm") + 
  geom_abline(slope = -1, linetype = "dashed") +
  geom_point() + 
  ggrepel::geom_text_repel(aes(label = Gene))

# look at dox treatment vs spen
all_coefs[,c("ndDox7", "ndDox7.ndAux9")] %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(x = ndDox7, y = ndDox7.ndAux9)) + 
  geom_smooth(method = "lm") + 
  geom_abline(slope = -1, linetype = "dashed") + 
  geom_point() + 
  ggrepel::geom_text_repel(aes(label = Gene))

# example genes
i = "Kdm6a"
y = as.numeric(data_test_inactive[i, ])
N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
fit1 <- glm(cbind(y, N-y) ~ Clone + ndDox * ndAux * washout, family = "binomial", data = fitting_metadata_here)

summary(fit1)
coef(fit1)

data.frame(coef(fit1)) %>%
  rownames_to_column("Variable") %>%
  ggplot(aes(Variable, coef.fit1.)) + geom_bar(stat = "identity") +
  coord_flip()

data.frame(
  y = y,
  N = N,
  Cov = metadata_b6_here$Condition,
  Clone = metadata_b6_here$Clone
) %>% ggplot(aes(Cov, y / N, col = Clone)) + geom_point() + coord_flip() + ylim(0, 1) +
  ggtitle(i) +
  geom_hline(yintercept = rev_logit(coef(fit1)[[1]]), linetype = "dashed")

### given a vector of coefficients, check if they correlate with any genomic / epigenomic features 
### from https://genome.cshlp.org/content/29/7/1087.full.pdf+html

# their escapee-annotation
load("~/Desktop/Projects/XChromosome_Antonia/Data/Supplemental_Code/data/annotation_files/escapees/escapees.RData")

# their features # called data-set
load("~/Desktop/Projects/XChromosome_Antonia/Data/Supplemental_Code/data/feature_matrix/promoter_matrix_reannotated_normRAdjusted_all_chrX_genes.RData")

# overlap between our gene list and theirs
both_genes <- intersect(rownames(data_set), rownames(all_coefs))

data_set_fit <- data_set[both_genes, ]
data_set_fit_scaled <- apply(data_set_fit, 2, function(x){x <- as.numeric(x); return((x - mean(x)) / sd(x))})
coefs_fit_here <- all_coefs[both_genes, ]

# we can use multivariate linear regression instead of RFs
library(randomForest)
require(pROC)

rf_cv_auc <- function(data_fit, folds = 5){     
  n_chunks <- floor(nrow(data_fit) / folds)
  data_fit <- data_fit[sample(nrow(data_fit)), ]
  df_folds <- split(data_fit, (seq(nrow(data_fit))-1) %/% n_chunks)
  lapply(1:(length(df_folds) - 1), function(i){
    train_data <- do.call('rbind', df_folds[-i])
    test_data <- df_folds[[i]]
    rf <- randomForest(
      objective ~ .,
      data = train_data,
      importance = T, 
    )
    prediction <- stats::predict(rf, test_data, type = "vote")
    rf.roc <- roc(test_data$objective, prediction[,2])
    return(auc(rf.roc))
  })
}

# thanks, https://www.r-bloggers.com/2019/07/clean-consistent-column-names/
clean_names <- function(.data, unique = FALSE) {
  n <- if (is.data.frame(.data)) colnames(.data) else .data
  n <- gsub("%+", "_pct_", n)
  n <- gsub("\\$+", "_dollars_", n)
  n <- gsub("\\++", "_plus_", n)
  n <- gsub("-+", "_minus_", n)
  n <- gsub("\\*+", "_star_", n)
  n <- gsub("#+", "_cnt_", n)
  n <- gsub("&+", "_and_", n)
  n <- gsub("@+", "_at_", n)
  n <- gsub("[^a-zA-Z0-9_]+", "_", n)
  n <- gsub("([A-Z][a-z])", "_\\1", n)
  n <- tolower(trimws(n))
  
  n <- gsub("(^_+|_+$)", "", n)
  
  n <- gsub("_+", "_", n)
  
  if (unique) n <- make.unique(n, sep = "_")
  
  if (is.data.frame(.data)) {
    colnames(.data) <- n
    .data
  } else {
    n
  }
}

run_rf_per_vector <- function(coef_vector){
  
  data_fit <- data_set_fit
  
  correlations <- apply(data_fit, 2, function(x){cor(coef_vector, as.numeric(x))})
  
  names_save <- colnames(data_fit)
  colnames(data_fit) <- paste0("Feature", 1:ncol(data_fit))
  
  data_fit <- cbind(data_fit, "objective" = coef_vector)
  
  rf <- randomForest(
    objective ~ .,
    data=data_fit,
    importance = T, 
  )
  
  importance_df <- data.frame(rf$importance)
  importance_df$cor <- correlations
  rownames(importance_df) <- names_save
  
  list(importance_df, rf$rsq, rf$mse)
}

rf_results <- apply(coefs_fit_here, 2, run_rf_per_vector)

rsqs <- unlist(lapply(rf_results, function(x){mean(x[[2]])}))
mses <- unlist(lapply(rf_results, function(x){mean(x[[3]])}))

i = 1
colnames(coefs_fit_here)[i]
testy <- rf_results[[i]][[1]]
head(testy[order(abs(testy$cor), decreasing = T), ], n = 20)

# combine correlations

all_cors_heatmap <- do.call("rbind", lapply(rf_results, function(x){x[[1]][,"cor"]}))
colnames(all_cors_heatmap) <- colnames(data_set_fit)

pheatmap(all_cors_heatmap, cluster_rows = F, color=colorRampPalette(c("navy", "white", "red"))(50))

# 

# 
# # exclude the samples that aren't present in both clones
# 
# table(metadata_b6$Condition)
# 
# samples_in_b6 <- metadata_b6[!metadata_b6$Condition %in% c("Aux_9_Dox_7_WO_NO_WOAuxNA", "Aux_9_Dox_7_WO_YES_WOAuxYES"), ]$SampleName
# samples_in_cast <- metadata_cast[!metadata_cast$Condition %in% c("Aux_9_Dox_7_WO_NO_WOAuxNA", "Aux_9_Dox_7_WO_YES_WOAuxYES"), ]$SampleName
# 
# data.b6.X.here <- data.b6.X[,samples_in_b6]
# data.cast.X.here <- data.cast.X[,samples_in_cast]
# 
# fitting_metadata_here <- metadata_b6[!metadata_b6$Condition %in% c("Aux_9_Dox_7_WO_NO_WOAuxNA", "Aux_9_Dox_7_WO_YES_WOAuxYES"), 
#                                      c("ndDox", "ndAux", "washout", "WOwithAux")]
# fitting_metadata_here$ndDox <- round((fitting_metadata_here$ndDox) / 3)
# fitting_metadata_here$ndAux <- round((fitting_metadata_here$ndAux) / 3)
# 
# # fitting_metadata_here <- data.frame(apply(fitting_metadata_here, 2, function(x){factor(x)}))
# 
# genes_test <- rownames(data.b6.X.here)
# genes_test <- genes_test[!is.na(genes_test)]
# 
# data_test_inactive <- data.b6.X.here[genes_test, ]
# data_test_active <- data.cast.X.here[genes_test, ]
# 
# # first check which genes we can actually explain any additional variance
# 
# all_coefs <- data.frame(do.call("rbind", lapply(1:nrow(data_test_active), function(i){
#   tryCatch({
#     y = as.numeric(data_test_inactive[i, ])
#     N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
#     # fit1 <- glm(cbind(y, N-y) ~ 1, family = "binomial")
#     fit1 <- glm(cbind(y, N-y) ~ 1, family = "binomial", data = fitting_metadata_here)
#     fit2 <- glm(cbind(y, N-y) ~ 1 + ndDox + ndAux + washout + WOwithAux, family = "binomial", data = fitting_metadata_here)
#     # kkk <- lrtest(fit2, fit1)
#     return(c(logLik(fit1), logLik(fit2)))
#   },error=function(cond) {
#     return(c(NA))
#   })
# })))
# 
# all_coefs$genes <- rownames(data_test_active)
# colnames(all_coefs) <- c("LL_null", "LL_alt", "Gene")
# all_coefs$LLR <- all_coefs$LL_alt - all_coefs$LL_null
# 
# # example genes
# i = "Gm14639"
# y = as.numeric(data_test_inactive[i, ])
# N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
# fit1 <- glm(cbind(y, N-y) ~ 1 + ndDox + ndAux + washout + WOwithAux, family = "binomial", data = fitting_metadata_here)
# 
# summary(fit1)
# coef(fit1)
# 
# data.frame(
#   y = y, 
#   N = N, 
#   Cov = metadata_b6[!metadata_b6$Condition %in% c("Aux_9_Dox_7_WO_NO_WOAuxNA", "Aux_9_Dox_7_WO_YES_WOAuxYES"), ]$Condition, 
#   Clone = metadata_b6[!metadata_b6$Condition %in% c("Aux_9_Dox_7_WO_NO_WOAuxNA", "Aux_9_Dox_7_WO_YES_WOAuxYES"), ]$Clone
# ) %>% ggplot(aes(Cov, y / N, col = Clone)) + geom_point() + coord_flip() + ylim(0, 1) + 
#   ggtitle(i)
# 
# genes_test <- rowMeans(data.b6.X.here / (data.b6.X.here + data.cast.X.here)) > 0.05 & rowMeans(data.b6.X.here / (data.b6.X.here + data.cast.X.here)) < 0.8
# genes_test <- genes_test[!is.na(genes_test)]
# 
# all_coefs$gene_empirical <- all_coefs$Gene %in% names(which(genes_test))
# ggplot(all_coefs, aes(x = gene_empirical, y = LLR)) + geom_boxplot() + scale_y_log10(limits = c(10, 5000))
# 
# data_test_inactive <- data.b6.X.here[names(which(genes_test)), ]
# data_test_active <- data.cast.X.here[names(which(genes_test)), ]
# 
# fitting_metadata_here$ndDox <- factor(fitting_metadata_here$ndDox, levels = unique(fitting_metadata_here$ndDox))
# fitting_metadata_here$ndAux <- factor(fitting_metadata_here$ndAux, levels = unique(fitting_metadata_here$ndAux))
# fitting_metadata_here$washout <- factor(fitting_metadata_here$washout, levels = c("NO", "YES"))
# fitting_metadata_here$WOwithAux <- factor(fitting_metadata_here$washout, levels = c("NO", "YES"))
# 
# all_coefs <- data.frame(do.call("rbind", lapply(1:nrow(data_test_active), function(i){
#   tryCatch({
#     y = as.numeric(data_test_inactive[i, ])
#     N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
#     fit1 <- glm(cbind(y, N-y) ~ ndDox + ndAux + washout + WOwithAux, family = "binomial", data = fitting_metadata_here)
#     coefs <- coef(fit1)
#    return(coefs)
#   }, error = function(cond) {
#     return(c(NA))
#   }, warning = function(cond){
#     return(c(NA))
#   })
# })))
# 
# rownames(all_coefs) <- rownames(data_test_active)
# rownames(all_coefs) <- rownames(data_test_inactive)
# all_coefs <- all_coefs[!apply(all_coefs, 1, function(x){all(is.na(x))}), ]
# head(all_coefs[order(all_coefs$ndDox, decreasing = F), ])
# 
# # check what intercepts mean
# # rev_logit(intercept) is ~ average allelic imbalance
# plot(rowMeans(data_test_inactive / (data_test_inactive + data_test_active))[rownames(all_coefs)], rev_logit(all_coefs$X.Intercept.))
# 
# pheatmap::pheatmap(all_coefs[,-1], cluster_cols = F, color=colorRampPalette(c("navy", "white", "red"))(50))
# 
# # all_coefs <- all_coefs[!is.na(all_coefs$ndDox0), ]
# 
# all_coefs %>%
#   pivot_longer(-c()) %>%
#   ggplot(aes(x = name, y = value)) + geom_boxplot()
# 
# head(all_coefs[order(all_coefs$washout, decreasing = T), ], n = 20)
# 
# # example genes
# i = "Kdm6a"
# y = as.numeric(data_test_inactive[i, ])
# N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
# fit1 <- glm(cbind(y, N-y) ~ ndDox + ndAux + washout + WOwithAux, family = "binomial", data = fitting_metadata_here)
# # fit1 <- glm(cbind(y, N-y) ~ ndDox + ndAux + washout, family = "binomial", data = fitting_metadata_here)
# 
# summary(fit1)
# coef(fit1)
# 
# data.frame(coef(fit1)) %>%
#   rownames_to_column("Variable") %>%
#   ggplot(aes(Variable, coef.fit1.)) + geom_bar(stat = "identity") + 
#   coord_flip()
# 
# data.frame(
#   y = y, 
#   N = N, 
#   Cov = metadata_b6[!metadata_b6$Condition %in% c("Aux_9_Dox_7_WO_NO_WOAuxNA", "Aux_9_Dox_7_WO_YES_WOAuxYES"), ]$Condition, 
#   Clone = metadata_b6[!metadata_b6$Condition %in% c("Aux_9_Dox_7_WO_NO_WOAuxNA", "Aux_9_Dox_7_WO_YES_WOAuxYES"), ]$Clone
# ) %>% ggplot(aes(Cov, y / N, col = Clone)) + geom_point() + coord_flip() + ylim(0, 1) + 
#   ggtitle(i) + 
#   geom_hline(yintercept = rev_logit(coef(fit1)[[1]]), linetype = "dashed")

### ------------### ------------
### differential testying
### ------------### ------------


# head(metadata)
# 
# all_samples_b6 <- c("CL30.9dAUX_C57BL.6J_sorted.bam", "CL31.9dAUX_C57BL.6J_sorted.bam", "CL31.NoDoxNoAux_C57BL.6J_sorted.bam", "CL30.NoDoxNoAux_C57BL.6J_sorted.bam")
# all_samples_cast <- c("CL30.9dAUX_CAST.EiJ_sorted.bam", "CL31.9dAUX_CAST.EiJ_sorted.bam", "CL31.NoDoxNoAux_CAST.EiJ_sorted.bam", "CL30.NoDoxNoAux_CAST.EiJ_sorted.bam")
# 
# data_test_inactive <- data.b6.X[names(genes_include[which(genes_include)]), all_samples_b6]
# data_test_active <- data.cast.X[, all_samples_cast]
# rownames(data_test_active) <- rownames(data.b6.X)
# data_test_active <- data_test_active[names(genes_include[which(genes_include)]), ]
# 
# groups <- c("group1", "group1", "group2", "group2")
# 
# library(VGAM)
# all_coefs <- data.frame(do.call("rbind", lapply(1:nrow(data_test_active), function(i){
#   tryCatch({
#     y = as.numeric(data_test_inactive[i, ])
#     N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
#     fit1 <- vglm(cbind(y, N-y) ~ 1, betabinomialff, trace = TRUE)
#     fit2 <- vglm(cbind(y, N-y) ~ groups, betabinomialff, trace = TRUE)
#     kkk <- lrtest(fit2, fit1)
#     return(c(i, kkk@Body$`Pr(>Chisq)`[[2]]))
#   },error=function(cond) {
#     return(c(NA))
#   })
# })))
# 
# colnames(all_coefs) <- c("Index", "pval")
# all_coefs$Gene <- rownames(data_test_active)
# all_coefs$padj <- p.adjust(all_coefs$pval)
# all_coefs <- all_coefs[!is.na(all_coefs$pval), ]
# head(all_coefs[order(all_coefs$padj), ])
# 
# # use MOFA for dimensionality reduction
# library(MOFA2)
# 
# # get variable genes in expression space
# expression_data <- assays(des.vst)[[1]][HVG(assays(des.vst)[[1]], numberGenes = 1000), ]
# rownames(expression_data) <- make.unique(rowData(des.vst)$symbol[HVG(assays(des.vst)[[1]], numberGenes = 1000)])
# 
# # allelic x ratios
# colnames(data.ratios) <- colnames(expression_data)
# 
# MOFAobject <- create_mofa(data = list("View1" = as.matrix(expression_data), "View2" = as.matrix(data.ratios)))
# MOFAobject <- prepare_mofa(MOFAobject)
# MOFAobject <- run_mofa(MOFAobject)
# 
# plot_variance_explained(MOFAobject)
# plot_weights(MOFAobject, factors = 3, view = 2)
