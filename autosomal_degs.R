setwd("~/Desktop/Projects/XChromosome_Antonia/")

library(tidyverse)
library(ggplot2)
library(scran)
library(scater)
library(DESeq2)

theme_paper <- function(){
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        text = element_text(size = 20))
}

source("./Scripts/auxiliary.R")

# Read known information about escapee-status. 
# This is based on a literature survey and classifies genes into constitutive, facultative and non-escaping
# apart from Xist (which isn't really an escapee bc it's only from the Xi)
escapee_annotation <- readxl::read_excel("~/Desktop/Projects/XChromosome_Project/ProcessedData/ListOfEscapeeFromEdithLab.xlsx") %>%
  dplyr::select(c("ENSEMBL_v102", "final status"))
library("EnsDb.Mmusculus.v79") 
symbols <- mapIds(EnsDb.Mmusculus.v79, keys = escapee_annotation$ENSEMBL_v102, keytype = "GENEID", column="SYMBOL")
escapee_annotation$symbol <- symbols
convert_status_names <- setNames(c("constitutive", "facultative", "variable"), c("E", "V", "S"))
convert_status_names2 <- setNames(c("constitutive", "facultative", "silenced / variable"), c("E", "V", "S"))
escapee_annotation$our_status <- unlist(convert_status_names[escapee_annotation$`final status`])
escapee_annotation$our_status2 <- unlist(convert_status_names2[escapee_annotation$`final status`])

#### new analysis - autosomal + compensation

data <- readRDS("./ProcessedData/merged_dataset.rds")
data <- data[,!data$Guide %in% c("g13", "P3") & data$Experiment != "May2022"]

data[,data$Experiment == "September2022"]$Condition

# clean up names
name_conversion_condition <- list(
  "Aux_0_Dox_0_WO_0_WOAuxNO" = "Control", 
  
  "Aux_0_Dox_3_WO_0_WOAuxNO" = "Dox (3d)", 
  "Aux_0_Dox_3_WO_4_WOAuxNO" = "Dox (3d) - washout (4d)", 
  
  "Aux_0_Dox_7_WO_0_WOAuxNO" = "Dox (7d)", 
  "Aux_0_Dox_7_WO_4_WOAuxNO" = "Dox (7d) - washout (4d)", 
  "Aux_0_Dox_7_WO_7_WOAuxNO" = "Dox (7d) - washout (7d)",
  "Aux_0_Dox_7_WO_14_WOAuxNO" = "Dox (7d) - washout (14d)",
  "Aux_0_Dox_7_WO_4_WOAuxYES" = "Dox (7d) - washout (with Aux)",
  
  "Aux_0_Dox_14_WO_0_WOAuxNO" = "Dox (14d)",
  "Aux_0_Dox_14_WO_7_WOAuxNO" = "Dox (14d) - washout (7d)",
  "Aux_0_Dox_14_WO_14_WOAuxNO" = "Dox (14d) - washout (14d)",
  
  "Aux_0_Dox_21_WO_0_WOAuxNO" = "Dox (21d)",
  "Aux_0_Dox_21_WO_7_WOAuxNO" = "Dox (21d) - washout (7d)",
  "Aux_0_Dox_21_WO_14_WOAuxNO" = "Dox (21d) - washout (14d)",
  
  "Aux_2_Dox_0_WO_0_WOAuxNO" = "Aux (2d)", 
  "Aux_2_Dox_0_WO_4_WOAuxNO" = "Aux (2d) - washout (4d)", 
  "Aux_5_Dox_0_WO_0_WOAuxNO" = "Aux (5d)", 
  "Aux_5_Dox_3_WO_0_WOAuxNO" = "Dox (3d), Aux (5d)", 
  "Aux_5_Dox_3_WO_4_WOAuxNO" = "Dox (3d) - washout (4d), Aux (5d)", 
  "Aux_9_Dox_0_WO_0_WOAuxNO" = "Aux (9d)", 
  "Aux_9_Dox_7_WO_0_WOAuxNO" = "Dox (7d), Aux (9d)", 
  "Aux_9_Dox_7_WO_4_WOAuxNO" = "Dox (7d) - washout (4d), Aux (9d)", 
  "Aux_9_Dox_7_WO_4_WOAuxYES" = "Dox (7d) - washout (4d) (with Aux), Aux (9d)"
)


data$ConditionClean <- unlist(name_conversion_condition[data$Condition])

# remove genes with < 10 allelic reads per sample 
data <- data[(rowSums(counts_inactive(data)) + rowSums(counts_active(data))) / ncol(data) > 10, ]

# subset on X-linked genes
data <- data[seqnames(data) == "X", ]

## Verify that change in allelic balance comes from increase in Xi, not decrease in Xa: 
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_4_WOAuxNO", 
                                         "Aux_0_Dox_7_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_4_WOAuxNO", "Aux_0_Dox_14_WO_0_WOAuxNO", 
                                         "Aux_0_Dox_21_WO_4_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_4_WOAuxNO", 
                                         "Aux_0_Dox_14_WO_4_WOAuxNO", "Aux_0_Dox_21_WO_YES_14_WOAuxNO", "Aux_0_Dox_0_WO_0_WOAuxNO", 
                                         "Aux_0_Dox_21_WO_0_WOAuxNO", "Aux_0_Dox_14_WO_4_WOAuxNO")]

data_here <- data_here[,data_here$Experiment == "September2022"]
data_here$sizeFactor <- calculateSumFactors(data_here)

plot_gene_total <- function(gene, save = F, path = NULL, suffix = "_escape_plot.pdf"){
  y = as.numeric(assays(data_here)[["counts_active"]][gene, ])
  N = as.numeric(assays(data_here)[["counts_active"]][gene, ] + assays(data_here)[["counts_inactive"]][gene, ])
  size_factors = colData(data_here)$sizeFactor
  # fit1 <- glm(cbind(y, N-y) ~ Clone + ndDox + washout, family = "binomial", data = fitting_metadata_here)
  # fit0 <- glm(cbind(y, N-y) ~ 1, family = "binomial", data = fitting_metadata_here)
  
  pp <- data.frame(
    active = (y + 1) * size_factors,
    inactive = (N - y + 2) * size_factors,
    Dox = data_here$ndDox,
    Washout = data_here$washout, 
    Clone = paste0(data_here$Clone, " - ", data_here$Experiment)
  ) %>%
    add_column(TimeSeries = ifelse(.$Dox == 0, list(c("3", "7", "14", "21")), .$Dox)) %>%
    unnest(TimeSeries) %>%
    mutate(TimeSeries = paste0(unlist(TimeSeries), " days dox treatment")) %>%
    add_column(Condition = ifelse(.$Dox == 0, "Control", .$Washout)) %>%
    mutate(Condition = factor(Condition, levels = rev(c("Control", "NO", "YES", "YES_long")))) %>%
    pivot_longer(cols = c("inactive", "active")) %>%
    mutate(group = paste0(Clone, "_", name)) %>%
    ggplot(aes(x = Condition, y = value, col = Clone, shape = name)) + 
    geom_point(size = 4) + 
    geom_line(aes(group = group)) + 
    coord_flip() + 
    theme_paper() + 
    facet_wrap(~TimeSeries, nrow = 2) + 
    xlab("") + ylab("ASE (d-score)") + 
    #scale_x_discrete(labels = c("Dox + washout", "Dox only", "Control")) + 
    scale_y_log10()
  pp
  
  if (save){
    ggsave(paste0(path, "./", gene, suffix), pp)
  } 
  
  return(pp)
}

plot_gene_total("Idh3g")

# We need to do this systematically - use DESeq on both matrices with the 3 replicates between Control and 7d
data_here <- data[,data$Clone == "E6" & data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO")]

# check ratios: 

# get d-scores across genes + samples
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

# exclude genes with 0 counts in individual samples
rownames( ratios[!apply(ratios, 1, function(x){!any(is.na(x))}), ] )
ratios <- ratios[apply(ratios, 1, function(x){!any(is.na(x))}), ]

#colnames(ratios) <- c("CL30", "CL31", "E6 (rep1)", "E6 (rep2)")

genes_show <- c("Xist", "Mecp2", "Kdm6a", "Kdm5c")
genes_of_choice <- setNames(rep("", nrow(ratios)), rownames(ratios))
genes_of_choice[genes_show] <- genes_show

# Plot d-scores for genes after ordering the genes by chromosomale coordinate
ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  ggplot(aes(Clone, y = reorder(name, position), fill = value)) + geom_tile(width = 0.9) +
  theme_paper() + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
  theme(axis.ticks.y = element_blank()) + xlab("") + ylab("Chromosome position") + 
  scale_fill_gradientn(colors = c("#808000", "pink", "#008080")) + 
  scale_y_discrete(labels = genes_of_choice, position = "right") + ylab("") + 
  facet_wrap(~Condition)

# get sizefactors based on total counts of the libraries, as a control also get sizefactors from active + inactive x: 
data_for_sf <- readRDS("./ProcessedData/merged_dataset.rds")
data_for_sf <- data_for_sf[,data_for_sf$Clone == "E6" & data_for_sf$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_4_WOAuxNO")]
data_for_sf <- data_for_sf[(rowSums(counts_inactive(data_for_sf)) + rowSums(counts_active(data_for_sf))) / ncol(data_for_sf) > 10, ]

size_factors_global <- DESeq2::estimateSizeFactors(DESeqDataSetFromMatrix(counts(data_for_sf), colData = DataFrame(Cond = rep(1, 6)), design = ~0))
size_factors_x <- DESeq2::estimateSizeFactors(DESeqDataSetFromMatrix(counts_active(data_here) + counts_inactive(data_here), colData = DataFrame(Cond = rep(1, 6)), design = ~0))

deseq_inactive <- DESeqDataSetFromMatrix(countData = counts_inactive(data_here), colData = colData(data_here), design = ~ Experiment + ConditionClean)
deseq_inactive$sizeFactor <- size_factors_global$sizeFactor
deseq_inactive <- DESeq(deseq_inactive)

deseq_active <- DESeqDataSetFromMatrix(countData = counts_active(data_here), colData = colData(data_here), design = ~ Experiment + ConditionClean)
deseq_active$sizeFactor <- size_factors_global$sizeFactor
deseq_active <- DESeq(deseq_active)

# xist plot for sanity check: 

data.frame(
  xist_expression = counts(deseq_inactive)["Xist", ] / deseq_inactive$sizeFactor, 
  condition = deseq_inactive$Condition
) %>%
  ggplot(aes(x = condition, y = xist_expression)) + geom_boxplot() + geom_jitter(width = 0.1)

results(deseq_inactive) %>%
  data.frame() %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(x = baseMean, y = log2FoldChange, col = padj < 0.1)) + geom_point() + scale_x_log10() + 
  ggrepel::geom_text_repel(aes(label = Gene)) + ylim(c(-12, 3))

results(deseq_active) %>%
  data.frame() %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(x = baseMean, y = log2FoldChange, col = padj < 0.1)) + geom_point() + scale_x_log10() + 
  ggrepel::geom_text_repel(aes(label = Gene)) + ylim(c(-12, 3)) + 
  NULL

results(deseq_active) %>%
  data.frame() %>%
  dplyr::arrange(padj) %>% 
  head(n = 20)

data.frame(
  active = results(deseq_active)$log2FoldChange, 
  inactive = results(deseq_inactive)$log2FoldChange
) %>% ggplot(aes(x = active, y = inactive)) + geom_point()

deseq_active <- logNormCounts(deseq_active)
data.frame(
  Expression = as.numeric(assays(deseq_active)[["logcounts"]]["Firre", ]), 
  Condition = deseq_active$Condition, 
  Experiment = deseq_active$Experiment
) %>% 
  ggplot(aes(x = Condition, y = Expression, col = Experiment)) + geom_point()

# DESeq analysis across autosomes

# Now read our pre-processed dataset
data <- readRDS("./ProcessedData/merged_dataset.rds")
data <- data[rowMeans(counts(data)) > 50, ]
data <- computeSumFactors(data)

data_here <- data[,data$Clone == "E6" & data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO")]

deseq_autosomes <- DESeqDataSetFromMatrix(countData = counts(data_here), colData = colData(data_here), design = ~ Experiment + Condition)
deseq_autosomes <- DESeq(deseq_autosomes)
deseq_autosomes <- logNormCounts(deseq_autosomes)
deseq_autosomes_res <- results(deseq_autosomes) %>% data.frame()

head(deseq_autosomes_res[order(deseq_autosomes_res$padj), ], n = 20)

deseq_autosomes_res %>%
  dplyr::arrange(padj) %>%
  dplyr::filter(padj < 0.1)

results(deseq_autosomes) %>%
  data.frame() %>%
  rownames_to_column("Gene") -> df_plot

df_plot %>%
  ggplot(aes(x = baseMean, y = log2FoldChange, col = padj < 0.1)) + geom_point() + scale_x_log10() + 
    theme_paper() + xlab("Mean Expression") + ylab("log2( Xist / Control )") + scale_color_manual(values = c("grey", "red")) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    ggrepel::geom_text_repel(data = df_plot[df_plot$padj < 0.1, ], aes(label = Gene)) + 
    NULL
ggsave("./FiguresIllustrator/FigS1/DEGs_autosomal_Xist.pdf")

results(deseq_autosomes) %>%
  data.frame() %>%
  rownames_to_column("Gene") %>%
  dplyr::filter(padj < 0.1) %>%
  ggplot(aes(x = reorder(Gene, log2FoldChange), y = log2FoldChange)) + geom_point() + coord_flip()

data.frame(
  expression = assays(deseq_autosomes)[["logcounts"]]["Kdm6a", ], 
  condition = deseq_inactive$Condition, 
  experiment = deseq_inactive$Experiment
) %>% ggplot(aes(x = condition, y = expression)) + geom_boxplot(outlier.color = NA) + geom_jitter(width = 0.1, aes(col = experiment))

source("~/Desktop/PhD/RCode/Multipurpose/R_gsea_functions.R")

genesets_here <- get_gs_database_msigdb(which_set = "C5")
gsea.results <- run_gage_per_exp_preranked(exp.fc = setNames(deseq_autosomes_res$log2FoldChange, nm = rownames(deseq_autosomes_res)), genesets_here, is_human = F)

rbind(head(gsea.results$upregulated, n = 10), head(gsea.results$downregulated, n = 10)) %>%
  rownames_to_column("GO_term") %>%
  ggplot(aes(x = reorder(GO_term, q.val), size = -log10(q.val), y = stat.mean)) + geom_point() + coord_flip() + 
    xlab("") + ylab("Enrichment statistic") + geom_hline(linetype = "dashed", yintercept = 0) + 
    theme_paper(textsize = 20)
ggsave("./FiguresIllustrator/FigS1/GSEA_autosomal_xist.pdf")

# same including washout

data_here_washout <- data[,data$Clone == "E6" & data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_4_WOAuxNO")]
deseq_autosomes_washout <- DESeqDataSetFromMatrix(countData = counts(data_here_washout), colData = colData(data_here_washout), design = ~ Experiment +  Condition)
deseq_autosomes_washout <- DESeq(deseq_autosomes_washout)
deseq_autosomes_washout <- logNormCounts(deseq_autosomes_washout)
deseq_autosomes_washout_res <- results(deseq_autosomes_washout) %>% data.frame()

head(deseq_autosomes_washout_res[order(deseq_autosomes_washout_res$padj), ], n = 20)

deseq_autosomes_washout_res %>%
  dplyr::arrange(padj) %>%
  head(n = 20)

results(deseq_autosomes_washout) %>%
  data.frame() %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(x = baseMean, y = log2FoldChange, col = padj < 0.1)) + geom_point() + scale_x_log10() + 
  NULL

data.frame(
  expression = assays(deseq_autosomes_washout)[["logcounts"]]["Id4", ], 
  condition = deseq_autosomes_washout$Condition, 
  experiment = deseq_autosomes_washout$Experiment
) %>%
  ggplot(aes(x = condition, y = expression)) + geom_boxplot(outlier.color = NA) + geom_jitter(width = 0.1, aes(col = experiment), size = 3) + 
    theme_paper() + xlab("") + ylab("log2(Expression)") + ylim(c(0, 12))

# for degs in 7d vs cont, plot fcs 

deseq_autosomes_res <- deseq_autosomes_res[!is.na(deseq_autosomes_res$padj), ]

genes_regulated <- rownames(deseq_autosomes_res[deseq_autosomes_res$padj < 0.2, ])

data_here_washout <- logNormCounts(data_here_washout)

data.frame(
  gene = genes_regulated, 
  fcs_cont_xist = rowMeans(logcounts(data_here_washout[,data_here_washout$Condition == "Aux_0_Dox_7_WO_0_WOAuxNO"]))[genes_regulated] - 
    rowMeans(logcounts(data_here_washout[,data_here_washout$Condition == "Aux_0_Dox_0_WO_0_WOAuxNO"]))[genes_regulated], 
  fcs_xist_washout = rowMeans(logcounts(data_here_washout[,data_here_washout$Condition == "Aux_0_Dox_7_WO_0_WOAuxNO"]))[genes_regulated] - 
    rowMeans(logcounts(data_here_washout[,data_here_washout$Condition == "Aux_0_Dox_7_WO_4_WOAuxNO"]))[genes_regulated]
) %>%
  ggplot(aes(x = fcs_cont_xist, y = fcs_xist_washout)) + geom_point() + geom_abline(slope = 1, linetype = "dashed") + 
  theme_bw() + 
  ggrepel::geom_text_repel(aes(label = gene)) + 
  geom_hline(linetype = "dashed", yintercept = 0) + geom_vline(linetype = "dashed", xintercept = 0) + 
  xlab("Fold Change Xist / WT") + ylab("Fold Change Washout / Xist")
ggsave("./FiguresIllustrator/FigS1/DEG_autosomal_FC_comparison_all_genes.pdf")

genes_regulated <- rownames(deseq_autosomes_res)

data.frame(
  gene = genes_regulated, 
  fcs_cont_xist = rowMeans(logcounts(data_here_washout[,data_here_washout$Condition == "Aux_0_Dox_7_WO_0_WOAuxNO"]))[genes_regulated] - 
    rowMeans(logcounts(data_here_washout[,data_here_washout$Condition == "Aux_0_Dox_0_WO_0_WOAuxNO"]))[genes_regulated], 
  fcs_xist_washout = rowMeans(logcounts(data_here_washout[,data_here_washout$Condition == "Aux_0_Dox_7_WO_0_WOAuxNO"]))[genes_regulated] - 
    rowMeans(logcounts(data_here_washout[,data_here_washout$Condition == "Aux_0_Dox_7_WO_4_WOAuxNO"]))[genes_regulated]
) %>%
  ggplot(aes(x = fcs_cont_xist, y = fcs_xist_washout)) + geom_point() + geom_abline(slope = 1, linetype = "dashed") + 
  theme_bw() + 
  #ggrepel::geom_text_repel(aes(label = gene)) + 
  geom_hline(linetype = "dashed", yintercept = 0) + geom_vline(linetype = "dashed", xintercept = 0) + 
  xlab("Fold Change Xist / WT") + ylab("Fold Change Washout / Xist")
ggsave("./FiguresIllustrator/FigS1/DEG_autosomal_FC_comparison_all_genes.pdf")

data.frame(
  gene = genes_regulated, 
  fcs_cont_xist = rowMeans(logcounts(data_here_washout[,data_here_washout$Condition == "Aux_0_Dox_7_WO_0_WOAuxNO"][,1]))[genes_regulated] - 
    rowMeans(logcounts(data_here_washout[,data_here_washout$Condition == "Aux_0_Dox_0_WO_0_WOAuxNO"]))[genes_regulated], 
  fcs_xist_washout = rowMeans(logcounts(data_here_washout[,data_here_washout$Condition == "Aux_0_Dox_7_WO_4_WOAuxNO"][,3]))[genes_regulated] - 
    rowMeans(logcounts(data_here_washout[,data_here_washout$Condition == "Aux_0_Dox_7_WO_0_WOAuxNO"]))[genes_regulated]
) %>%
  ggplot(aes(x = fcs_cont_xist, y = fcs_xist_washout)) + geom_point() + geom_abline(slope = -1, linetype = "dashed") + 
  theme_bw() + ggrepel::geom_text_repel(aes(label = gene)) + geom_hline(linetype = "dashed", yintercept = 0) + geom_vline(linetype = "dashed", xintercept = 0) + 
  xlab("Fold Change Xist / WT") + ylab("Fold Change Washout / Xist")


### how many X-linked genes are when considering full counts?

# Now read our pre-processed dataset
data <- readRDS("./ProcessedData/merged_dataset.rds")
data <- data[rowMeans(counts(data)) > 50, ]
data <- computeSumFactors(data)

data_here <- data[,data$Clone == "E6" & data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO")]

deseq_autosomes <- DESeqDataSetFromMatrix(countData = counts(data_here), colData = colData(data_here), design = ~ Experiment + Condition)
deseq_autosomes <- DESeq(deseq_autosomes)
deseq_autosomes <- logNormCounts(deseq_autosomes)
deseq_autosomes_res <- results(deseq_autosomes) %>% data.frame()
deseq_autosomes_res$chromosome <- as.character(seqnames(rowRanges(data_here)[rownames(deseq_autosomes_res), ]))

deseq_autosomes_res %>% 
  dplyr::filter(padj < 0.1) %>% 
  add_column()


### what if we only test X genes?

deseq_autosomes <- DESeqDataSetFromMatrix(countData = counts(data_here), colData = colData(data_here), design = ~ Experiment + Condition)
deseq_autosomes <- deseq_autosomes[seqnames(rowRanges(data_here)) == "X", ]
deseq_autosomes <- DESeq(deseq_autosomes)
deseq_autosomes <- logNormCounts(deseq_autosomes)
deseq_autosomes_res <- results(deseq_autosomes) %>% data.frame()
deseq_autosomes_res$chromosome <- as.character(seqnames(rowRanges(data_here)[rownames(deseq_autosomes_res), ]))

deseq_autosomes_res %>% 
  dplyr::filter(padj < 0.1) %>% 
  add_column() %>%
  dplyr::arrange(padj)

deseq_autosomes_res %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(x = baseMean, y = log2FoldChange, col = padj < 0.1)) + geom_point() + scale_x_log10() + 
  NULL


### what if we only test escapees?

ratios <- counts_inactive(data_here[,data_here$Condition == "Aux_0_Dox_0_WO_0_WOAuxNO"]) / 
  (counts_active(data_here[,data_here$Condition == "Aux_0_Dox_0_WO_0_WOAuxNO"]) + counts_inactive(data_here[,data_here$Condition == "Aux_0_Dox_0_WO_0_WOAuxNO"]))

ratios <- ratios[as.character(seqnames(rowRanges(data_here))) == "X", ]

# Plot d-scores for genes that escape in any sample, excluding genes with > 0.7 d-score
ratios <- ratios[rownames(ratios) == "Xist" | (apply(ratios, 1, function(x){any(x > 0.1)}) & !apply(ratios, 1, function(x){any(x > 0.7)})), ]
ratios <- ratios[order(rowMeans(ratios), decreasing = T), ]

# Now we define the baseline escapees per sample as the genes with d-score > 0.1
escapees <- apply(ratios, 2, function(x){x > 0.1})
escapees_per_clone <- apply(escapees, 2, function(x){rownames(escapees)[x]})
escapees_per_clone <- lapply(escapees_per_clone, function(x){x[x != "Xist"]})
all_escapees <- intersect(escapees_per_clone$Xistlong_E6_NoDox_r1_sorted.bam, intersect(escapees_per_clone$C.E6.NoDox, escapees_per_clone$C.E6_NoDox_r2))
all_escapees <- all_escapees[!is.na(all_escapees)]

deseq_autosomes <- DESeqDataSetFromMatrix(countData = counts(data_here), colData = colData(data_here), design = ~ Experiment + Condition)
deseq_autosomes <- deseq_autosomes[all_escapees, ]
deseq_autosomes <- DESeq(deseq_autosomes)
deseq_autosomes <- logNormCounts(deseq_autosomes)
deseq_autosomes_res <- results(deseq_autosomes) %>% data.frame()
deseq_autosomes_res$chromosome <- as.character(seqnames(rowRanges(data_here)[rownames(deseq_autosomes_res), ]))

deseq_autosomes_res %>% 
  dplyr::filter(padj < 0.1) %>% 
  add_column() %>%
  dplyr::arrange(padj)

deseq_autosomes_res %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(x = baseMean, y = log2FoldChange, col = padj < 0.1)) + geom_point() + scale_x_log10() + 
  NULL




