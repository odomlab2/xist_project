setwd(paste0(dirname(rstudioapi::getSourceEditorContext()$path)))
setwd("../../../")

### joint analysis of DEG results
theme_paper <- function(){theme_classic(base_size = 20)}

## read deg lists and merge
data_e6_7d <- readRDS("./ProcessedData/autosomal_de_analysis_7d_dox.rds") %>% 
  add_column(Comparison = "E6_7d")

data_e6_21d <- readRDS("./ProcessedData/autosomal_de_analysis_21d_dox.rds") %>% 
  add_column(Comparison = "E6_21d")

data_cl30 <- readRDS("./ProcessedData/autosomal_de_analysis_cl30_7d_dox.rds") %>% 
  add_column(Comparison = "CL30_7d")

data_cl31 <- readRDS("./ProcessedData/autosomal_de_analysis_cl31_7d_dox.rds") %>% 
  add_column(Comparison = "CL31_7d")

data_astrocytes_3d <- readRDS("./ProcessedData/autosomal_de_analysis_astrocytes_3d_dox.rds") %>% 
  add_column(Comparison = "Astro_3d") %>% 
  dplyr::select(-c("Status"))

data_astrocytes_7d <- readRDS("./ProcessedData/autosomal_de_analysis_astrocytes_7d_dox.rds") %>% 
  add_column(Comparison = "Astro_7d") %>% 
  dplyr::select(-c("Status"))

datasets_include <- list(data_e6_7d, data_cl30, data_cl31, data_astrocytes_3d, data_astrocytes_7d)
reference_group <- c("E6_7d")

df_full <- do.call("rbind", datasets_include)

# only include genes with sufficient baseMean in all comparisons
df_full <- df_full %>%
  group_by(Gene) %>%
  mutate(coverage_threshold = sum(baseMean > 10)) %>%
  dplyr::filter(coverage_threshold == length(datasets_include)) %>%
  # mutate(ref_lfc = log2FoldChange[Comparison == reference_group]) %>%
  # mutate(ref_padj = padj[Comparison == reference_group]) %>%
  dplyr::select(-c("lfcSE", "stat", "coverage_threshold"))

# add chromosome information
genes <- genes(EnsDb.Mmusculus.v79) %>% data.frame() %>%
  dplyr::filter(gene_name %in% df_full$Gene) %>%
  dplyr::filter(!duplicated(gene_name)) %>%
  column_to_rownames("gene_name")

df_full <- df_full %>% 
  mutate(Chromosome = droplevels(genes[Gene, ]$seqnames))

df_full %>% 
  mutate(significant = padj < .1) %>%
  dplyr::filter(Chromosome != "X") %>%
  mutate(up_down = ifelse(log2FoldChange > 0, "up", "down")) %>%
  dplyr::filter(significant) %>%
  group_by(Comparison, up_down) %>%
  summarize(n = n())

## compare fold change correlations, exclude X-chromosome here
df_full %>%
  dplyr::filter(Chromosome != "X") %>%
  dplyr::select(-c(baseMean, pvalue, padj)) %>%
  pivot_wider(names_from = Comparison, values_from = log2FoldChange) -> df_ful_wide

full_join(
  pivot_longer(df_ful_wide, 
               cols = E6_7d:Astro_7d,
               names_to = "variable_x",
               values_to = "x"),
  pivot_longer(df_ful_wide, 
               cols = E6_7d:Astro_7d,
               names_to = "variable_y",
               values_to = "y"),
  by = c("Gene", "Chromosome"),
  relationship = "many-to-many"
) -> df_full_pairwise

df_full_pairwise %>%
  ggplot(aes(x = x, y = y)) + geom_point(size = .1) + facet_wrap(~variable_x+variable_y) + theme_paper() + 
    geom_abline(linetype = 'dashed')
ggsave("./Plots/FigsX_autosomal_joint_analysis/logfc_scatters.pdf", width = 14, height = 14)

df_full_pairwise %>%
  group_by(variable_x, variable_y) %>%
  summarize(cor = cor(x, y)) %>%
  dplyr::filter(cor < .99) %>%
  ggplot(aes(x = variable_x, y = variable_y)) + geom_point(aes(fill = cor), pch =21, size = 10) + 
    scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "white") + theme_paper() + xlab("") + ylab("") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("./Plots/FigsX_autosomal_joint_analysis/logfc_correlations.pdf", width = 7, height = 7)

df_full %>% 
  dplyr::filter(Chromosome != "X") %>%
  dplyr::filter(!is.na(padj)) %>%
  mutate(sig = ifelse(padj < .1, "sig", "not_sig")) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
  group_by(Comparison, sig, direction) %>%
  summarize(n = n()) %>% 
  dplyr::filter(sig != "not_sig")

df_full %>%
  dplyr::filter(Chromosome != "X") %>%
  dplyr::filter(padj < .1) %>%
  dplyr::filter(log2FoldChange > 0) %>%
  dplyr::filter(Comparison != "Astro_3d") %>%
  group_by(Comparison) %>%
  summarize(genes = list(Gene)) %>%
  pull(genes, name = "Comparison") %>%
  ggVennDiagram::ggVennDiagram(label = "count") + scale_fill_gradient(low = "white", high = "red")
ggsave("./Plots/FigsX_autosomal_joint_analysis/venn_genes_up.pdf")

df_full %>%
  dplyr::filter(Chromosome != "X") %>%
  dplyr::filter(padj < .1) %>%
  dplyr::filter(log2FoldChange < 0) %>%
  dplyr::filter(Comparison != "Astro_3d") %>%
  group_by(Comparison) %>%
  summarize(genes = list(Gene)) %>%
  pull(genes, name = "Comparison") %>%
  ggVennDiagram::ggVennDiagram(label = "count") + scale_fill_gradient(low = "white", high = "red")
ggsave("./Plots/FigsX_autosomal_joint_analysis/venn_genes_down.pdf")

df_full %>%
  dplyr::filter(Chromosome != "X") %>%
  dplyr::filter(padj < .1) %>%
  dplyr::filter(log2FoldChange < 0) %>%
  group_by(Gene) %>%
  summarize(n = n()) %>% 
  dplyr::filter(n > 2) %>%
  arrange(-n)

# run paired genes analysis (in astrocytes script)
df_full %>%
  mutate(chr = as.character(genes[Gene, ]$seqnames)) %>% 
  mutate(start = as.numeric(genes[Gene, ]$start)) %>%
  mutate(end = as.numeric(genes[Gene, ]$end)) -> df_full_position

calculate_gene_distances <- function(data, max_dist = 1e5){
  # calculate how often genes have another gene close by (atmost max_dist away)
  data %>% 
    group_by(chr) %>%
    mutate(
      close_gene = unlist(lapply(start, function(x){
        ifelse(any(start - x < max_dist & start - x > 0), "close", "far")
      }))
    ) %>%
    group_by(close_gene, .drop = F) %>% 
    summarize(n = n()) -> temp
  
  if (nrow(temp) == 1){
    temp <- rbind(c("close", 0), temp)
  }
  
  temp %>%
    pivot_wider(values_from = n, names_from = close_gene) %>% 
    mutate(close = as.numeric(close), far = as.numeric(far)) %>%
    mutate(ratio = close / (close + far))
  
}


compute_statistics <- function(res_test, cutoff = 0.2){
  
  lapply(1:10, function(i){
    res_test %>% 
      ungroup() %>%
      dplyr::filter(chr == "X") %>%
      mutate(padj = sample(padj)) %>%
      dplyr::filter(padj < cutoff) %>%
      calculate_gene_distances() %>% 
      add_column(it = i)
  }) %>% do.call("rbind", .) -> x_ratio_shuff
  
  lapply(1:10, function(i){
    res_test %>% 
      ungroup() %>%
      dplyr::filter(chr != "X") %>%
      mutate(padj = sample(padj)) %>%
      dplyr::filter(padj < cutoff) %>%
      calculate_gene_distances() %>% 
      add_column(it = i)
  }) %>% do.call("rbind", .) -> auto_ratio_shuff
  
  res_test %>% 
    dplyr::filter(chr == "X") %>%
    dplyr::filter(padj < cutoff) %>%
    calculate_gene_distances() %>% 
    add_column(it = 1) -> x_ratio
  
  res_test %>% 
    dplyr::filter(chr != "X") %>%
    dplyr::filter(padj < cutoff) %>%
    calculate_gene_distances() %>%
    add_column(it = 1) -> auto_ratio
  
  list(
    cbind("cond" = "x_ratio_shuff", x_ratio_shuff),
    cbind("cond" = "auto_ratio_shuff", auto_ratio_shuff),
    cbind("cond" = "x_ratio", x_ratio),
    cbind("cond" = "auto_ratio", auto_ratio)
  ) %>% do.call("rbind", .) %>%
    mutate(condition = factor(cond, levels = c("x_ratio", "x_ratio_shuff", "auto_ratio", "auto_ratio_shuff"))) %>%
    # group_by(condition) %>%
    # summarize(mean = mean(ratio), sd = sd(ratio)) %>% 
    mutate(comparison = unique(res_test$Comparison2)) -> results
  
  results
}

df_full_position %>%
  group_by(Comparison) %>%
  summarize(n_genes = sum(padj < .1, na.rm = T))

df_full_position %>%
  mutate(Comparison2 = Comparison) -> df_full_position_here
  
df_full_position_here %>% 
  dplyr::filter(Comparison %in% c("Astro_3d", "Astro_7d")) %>%
  group_by(Comparison) %>%
  group_map(~compute_statistics(.x)) -> statistics_results

df_full_position_here %>% 
  dplyr::filter(Comparison %in% c("Astro_3d", "Astro_7d", "CL31_7d")) %>%
  group_by(Comparison) %>%
  group_map(~compute_statistics(.x, .1)) -> statistics_results

do.call("rbind", statistics_results) %>%
  ggplot(aes(x = condition, y = ratio)) + 
    stat_summary(geom = "bar") +
    stat_summary(geom = "errorbar", fun.data = function(x) {mean_sdl(x, mult = 1)}) + 
    ggbeeswarm::geom_quasirandom(size = 5, pch = 21, fill = "grey") + 
    coord_flip() + 
    xlab("") + ylab("Fraction of genes with DE neighbour") + theme_paper() + facet_wrap(~comparison, scales = "free")
ggsave("./Plots/FigsX_autosomal_joint_analysis/gene_pairing_analysis.pdf", width = 12, height = 12)

# run go analysis
source("./xist_project/Scripts/Rscripts/R_gsea_functions.R")

run_gage_per_exp <- function(deseq.res, gs_database, is_human = T){
  print(unique(deseq.res$Comparison2))
  library(gage)
  exp.fc <- deseq.res$log2FoldChange
  names(exp.fc) <- deseq.res$Gene
  if (!is_human){
    #### this is just dummy to convert mouse to human, do thath later via gene ontology
    #### presumably gsea broad software does the same thing
    names(exp.fc) <- toupper(names(exp.fc))
  }
  fc.kegg.p <- gage(exp.fc, gsets = gs_database, ref = NULL, samp = NULL)
  res_great <- data.frame(fc.kegg.p$greater)
  res_lower <- data.frame(fc.kegg.p$less)
  list("upregulated" = res_great, "downregulated" = res_lower)
}

kegg.genesets <- get_gs_database_kegg()
hs.genesets <- get_gs_database_msigdb()
# hs.genesets.c1 <- get_gs_database_msigdb(which_set = "C1")
# hs.genesets.c5 <- get_gs_database_msigdb(which_set = "C5")

gs_results <- df_full %>%
  dplyr::filter(Chromosome != "X") %>%
  dplyr::select(c(Gene, log2FoldChange, Comparison)) %>%
  mutate(Comparison2 = Comparison) %>%
  group_by(Comparison) %>%
  group_map(~run_gage_per_exp(.x, is_human = F, gs_database = hs.genesets))

names <- c("Astro_3d", "Astro_7d", "CL30_7d", "CL31_7d", "E6_7d")
gs_results_proc <- lapply(1:length(gs_results), function(i){
  x = gs_results[[i]]
  rbind(
    cbind(x$upregulated, "direction" = "up"), 
    cbind(x$downregulated, "direction" = "down")
  ) %>% 
    add_column("Comparison" = names[[i]])
}) %>% do.call("rbind", .)

gs_results_proc %>% 
  rownames_to_column("geneset") %>%
  dplyr::filter(q.val < .1) %>%
  group_by(Comparison, direction) %>%
  group_map(~head(.x, n = 5), .keep = T) %>% do.call("rbind", .) %>%
  ggplot(aes(x = reorder(geneset, stat.mean), y = stat.mean, col = -log10(q.val))) + geom_point() + 
    facet_wrap(~Comparison, scales = "free_y", nrow = 1) + coord_flip() + 
    geom_hline(yintercept = 0)
ggsave("./Plots/FigsX_autosomal_joint_analysis/go_analysis.pdf", width = 20, height = 6)

## compare to cnr data
library(csaw)
library(DESeq2)
library(tidyverse)

theme_paper <- function(textsize = 20){
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        text = element_text(size = textsize, color = "black"), 
        axis.text = element_text(size = textsize, color = "black"), 
        axis.ticks = element_line(colour = 'black', size = 1), 
        axis.ticks.length = unit(.25, "cm"))
}

# setwd("/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/allele-specific_cut_and_run/")
test <- readRDS("./ProcessedData/cnr_processed/genebody_totalcounts_processed.rds")

metadata_corrected <- colData(test)

assays(test)[["counts_total"]] %>% 
  data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-c("Gene")) %>% 
  add_column(assay = metadata_corrected[.$name, ]$HistoneMark) %>%
  add_column(Condition = metadata_corrected[.$name, ]$Condition) %>%
  mutate(Condition = gsub("r2|r3", "", Condition)) %>%
  mutate(Condition = ifelse(Condition == "21dD", "21dDox", Condition)) %>%
  add_column(Replicate = metadata_corrected[.$name, ]$Replicate) %>%
  add_column(SizeFactor = metadata_corrected[.$name, ]$sizeFactor) -> long_data

long_data %>%
  dplyr::filter(Gene %in% genes$gene_id) %>%
  mutate(chr = genes[match(Gene, genes$gene_id), ]$seqnames) %>%
  mutate(Gene = genes[match(Gene, genes$gene_id), ]$symbol) -> long_data

long_data %>%
  # dplyr::filter(chr != "X") %>%
  group_by(Gene, assay) %>%
  mutate(average_coverage = mean(value)) %>%
  dplyr::filter(average_coverage > 20) %>%
  mutate(value_norm = value / SizeFactor) %>%
  ungroup() -> processed_data

processed_data %>%
  mutate(xchr = ifelse(chr == "X", "X", "Autosomes")) %>%
  group_by(assay, Condition, Gene, xchr) %>%
  summarize(mean = mean(value_norm)) %>%
  group_by(Gene, assay) %>%
  mutate(lfc = log2((mean + 1) / (mean[Condition == "NoDox"] + 1))) %>% 
  dplyr::filter(Condition != "NoDox") -> cnr_foldchanges

cnr_foldchanges %>%
  ggplot(aes(x = mean, y = lfc, col = xchr)) + geom_point(size = .1) + facet_wrap(~assay + Condition, nrow = 2) + scale_x_log10()

cnr_foldchanges %>%
  dplyr::select(-c("mean")) %>%
  dplyr::filter(Condition %in% c("7dDox")) %>%
  pivot_wider(names_from = assay, values_from = lfc) -> cnr_foldchanges_cut

df_full_position_here %>%
  right_join(cnr_foldchanges_cut) %>%
  pivot_longer(c("H2AK119ubi", "H3K27me3")) -> joint_data

joint_data %>%
  dplyr::filter(xchr == "Autosomes") %>%
  dplyr::filter(Comparison != "Astro_3d") %>%
  ggplot(aes(x = log2FoldChange, y = value)) +
  geom_point(size = .01) + 
    facet_wrap(~name+Comparison, nrow = 2) + theme_paper() + 
    xlab("log2Foldchange (RNA)") + ylab("log2Foldchange (CnR)")
ggsave("./Plots/FigsX_autosomal_joint_analysis/scatter_to_cnr.pdf", width = 10, height = 6)

joint_data %>%
  dplyr::filter(xchr == "Autosomes") %>%
  dplyr::filter(Comparison != "Astro_3d") %>%
  group_by(name, Comparison) %>%
  summarize(cor = cor(value, log2FoldChange, use = "complete.obs")) %>%
  ggplot(aes(x = Comparison, y = name)) + geom_point(aes(fill = cor), pch = 21, size = 10) + scale_fill_gradient2() + 
    theme_paper() + xlab("") + ylab("")
ggsave("./Plots/FigsX_autosomal_joint_analysis/cor_to_cnr.pdf", width = 10, height = 6)

joint_data %>%
  # dplyr::filter(xchr == "X") %>%
  dplyr::filter(Comparison != "Astro_3d") %>%
  ggplot(aes(x = log2FoldChange, y = value, col = xchr)) +
  geom_point(size = .01) + 
  facet_wrap(~name+Comparison, nrow = 2) + theme_paper() + 
  xlab("log2Foldchange (RNA)") + ylab("log2Foldchange (CnR)")
ggsave("./Plots/FigsX_autosomal_joint_analysis/scatter_to_cnr_xchr.pdf", width = 10, height = 6)

joint_data %>%
  # dplyr::filter(xchr == "X") %>%
  dplyr::filter(Comparison != "Astro_3d") %>%
  group_by(name, Comparison, xchr) %>%
  summarize(cor = cor(value, log2FoldChange, use = "complete.obs")) %>%
  ggplot(aes(x = Comparison, y = name, group = xchr)) + 
    geom_point(aes(fill = cor), pch = 21, size = 10, position = position_dodge(width = .4)) + scale_fill_gradient2() + 
    theme_paper() + xlab("") + ylab("")
ggsave("./Plots/FigsX_autosomal_joint_analysis/cor_to_cnr_xchr.pdf", width = 10, height = 6)

