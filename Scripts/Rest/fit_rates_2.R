setwd("~/Desktop/Projects/XChromosome_Antonia/")

library(tidyverse)
library(ggplot2)
library(scran)
library(scater)

theme_paper <- function(){
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        text = element_text(size = 20))
}

escape_colors <- setNames(c("grey", "darkgreen", "orange"), nm = c("silenced / variable", "facultative", "constitutive"))

source("./Scripts/auxiliary.R")

data <- readRDS("./ProcessedData/merged_dataset_paper.rds")
data_before_filtering <-  data[seqnames(data) == "X", ]

data <- data[(rowSums(counts_inactive(data)) + rowSums(counts_active(data))) / ncol(data) > 10, ]

data_with_autosomes <- data
data <- data[seqnames(data) == "X", ]

# define the genes which escape across experiments in baseline
data_here <- data[,data$ConditionClean == "Control"]
data_here <- data_here[,data_here$Clone == "E6"]

# get d-scores across genes + samples
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

# Genes with average ASE > 0.7 are likely mapping artefacts: 
genes_out <- names(rowMeans(ratios, na.rm = T)[rowMeans(ratios, na.rm = T) > 0.8 & rownames(ratios) != "Xist"])
genes_out_na <- rownames( ratios[!apply(ratios, 1, function(x){!any(is.na(x))}), ] )

genes_out <- c(genes_out, genes_out_na)

ratios <- ratios[!rownames(ratios) %in% genes_out, ]

escapees <- apply(ratios, 2, function(x){x > 0.1})
escapees_per_clone <- apply(escapees, 2, function(x){rownames(escapees)[x]})
escapees_per_clone <- lapply(escapees_per_clone, function(x){x[x != "Xist"]})

# define escapees we will look at
escapees_use <- as.character(do.call("c", escapees_per_clone)[!duplicated(do.call("c", escapees_per_clone))])





# 

### Now we look at Xist overexpression and how that silences escapees
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO", 
                                         "Aux_0_Dox_14_WO_0_WOAuxNO", "Aux_0_Dox_21_WO_0_WOAuxNO")]
data_here <- data_here[,data_here$Clone == "E6"]

exp_decay <- function(x, a, k) a * exp(-k*x) 
offset_exp_decay <- function(x, a, k, c) a * exp(-k*x) + c

fit_nls_curve_for_gene <- function(gene){
  print(gene)
  df_here <- data.frame(
    active = as.numeric(counts_active(data_here)[gene, ]) + 1, 
    inactive = as.numeric(counts_inactive(data_here)[gene, ]) + 1, 
    timepoint = as.numeric(data_here$ndDox), 
    experiment = data_here$Experiment
  ) %>% mutate(total = active + inactive) %>%
    mutate(rate = inactive / (active + inactive))
  tryCatch({
    rates.nls <- nls(rate ~ a * exp(- k * timepoint), data = df_here, algorithm = "port", start = c(a = 0.5, k = 0.1), lower = c(a = -Inf, k = 0))
    return(c(coef(rates.nls), modelr::rsquare(rates.nls, data = df_here), logLik(rates.nls)))
  },
  error=function(cond) {
    return(NA)
  })
}

fit_nls_curve_for_gene_offset <- function(gene){
  print(gene)
  df_here <- data.frame(
    active = as.numeric(counts_active(data_here)[gene, ]) + 1, 
    inactive = as.numeric(counts_inactive(data_here)[gene, ]) + 1, 
    timepoint = as.numeric(data_here$ndDox), 
    experiment = data_here$Experiment
  ) %>% mutate(total = active + inactive) %>%
    mutate(rate = inactive / (active + inactive))
  tryCatch({
    rates.nls.offset <- nls(rate ~ a * exp(- k * timepoint) + b, data = df_here, algorithm = "port", start = c(a = 0.5, k = 0.1, b = 0.1), lower = c(a = -Inf, k = 0, b = 0))
    return(c(coef(rates.nls.offset), modelr::rsquare(rates.nls.offset, data = df_here), logLik(rates.nls.offset)))
  }, 
  error=function(cond) {
    return(NA)
  })
}

general_escapees <- data_here[escapees_use, ]

all_coefs <- data.frame(do.call("rbind", lapply(rownames(general_escapees), fit_nls_curve_for_gene)))
all_coefs$gene <- rownames(general_escapees)
all_coefs$halflifes <- log(1/2) / ( - all_coefs$k )

# 
all_coefs_offset <- data.frame(do.call("rbind", lapply(rownames(general_escapees), fit_nls_curve_for_gene_offset)))
all_coefs_offset$gene <- rownames(general_escapees)
all_coefs_offset$halflifes <- log(1/2) / ( - all_coefs_offset$k )

# how many convergence errors?
table(is.na(all_coefs$a), is.na(all_coefs_offset$a)) # big issue here: genes with b = 0 run into convergence errors. not sure how to fix this.

all_coefs_offset$category <- rowData(data_here[all_coefs$gene, ])$EscapeAnnotation

all_coefs %>%
  ggplot(aes(x = halflifes, y = V3)) + geom_point() + xlab("Halflifes") + ylab("Regression R^2")

all_coefs_offset %>%
  ggplot(aes(x = halflifes, y = V4)) + geom_point() + xlab("Halflifes") + ylab("Regression R^2")

# Merge: 

merged_df <- data.frame(
  gene = all_coefs$gene, 
  category = all_coefs_offset$category, 
  a = all_coefs$a, 
  k = all_coefs$k, 
  r2 = all_coefs$V3, 
  loglik = all_coefs$V4, 
  halflifes = all_coefs$halflifes, 
  a_off = all_coefs_offset$a, 
  k_off = all_coefs_offset$k, 
  b_off = all_coefs_offset$b, 
  r2_off = all_coefs_offset$V4, 
  loglik_off = all_coefs_offset$V5,
  halflifes_off = all_coefs_offset$halflifes
) %>%
  mutate(BIC_native = 3 * log(ncol(data_here)) - 2 * loglik) %>%
  mutate(BIC_offset = 4 * log(ncol(data_here)) - 2 * loglik_off)

# Which genes are going out?

merged_df %>%
  replace(is.na(.), 0) %>%
  ggplot(aes(x = r2, y = r2_off)) + geom_jitter() + xlab("Regression R^2 native") + ylab("Regression R^2 offset")

# convergence error native model: 
merged_df %>%
  dplyr::filter(is.na(a))

merged_df %>%
  dplyr::filter(is.na(a_off))


# convergence error offset model: 
merged_df %>%
  dplyr::filter(is.na(a_off)) %>%
  dplyr::arrange(r2)


merged_df %>% 
  dplyr::filter(!(r2 > 0.5 & r2_off > 0.5)) %>% dim()

merged_df <- merged_df %>% 
  dplyr::filter(r2 > 0.5 & r2_off > 0.5)

merged_df %>%
  ggplot(aes(x = r2)) + geom_histogram()

merged_df %>%
  ggplot(aes(x = r2_off)) + geom_histogram()

merged_df %>%
  ggplot(aes(x = r2, r2_off, col = category)) + geom_point() + xlim(c(0.5, 1 )) + ylim(c(0.5, 1)) + geom_abline() +
    ggrepel::geom_text_repel(aes(label = gene))

saveRDS(merged_df, "./Data/processed/Fig1_model_fits.rds")


gene = "Firre"

df_here <- data.frame(
  active = as.numeric(counts_active(data_here)[gene, ]) + 1, 
  inactive = as.numeric(counts_inactive(data_here)[gene, ]) + 1, 
  timepoint = as.numeric(data_here$ndDox), 
  experiment = data_here$Experiment
) %>% mutate(total = active + inactive) %>% mutate(rate = inactive / (active + inactive))

df_native <- all_coefs %>% column_to_rownames("gene")
df_offset <- all_coefs_offset %>% column_to_rownames("gene")

data.frame(
  active = as.numeric(counts_active(data_here)[gene, ]) + 1, 
  inactive = as.numeric(counts_inactive(data_here)[gene, ]) + 1, 
  timepoint = as.numeric(data_here$ndDox), 
  experiment = data_here$Experiment
) %>% mutate(rates = inactive / (inactive + active)) %>%
  ggplot(aes(x = timepoint, y = rates, col = experiment)) + geom_point() + theme_paper() + 
  geom_function(fun = function(x, a = df_native[gene, ]$a, k = df_native[gene, ]$k) exp_decay(x, a, k), color = "red", linetype = "dashed") + ylab("allelic ratio") + 
  geom_function(fun = function(x, a = df_offset[gene, ]$a, k = df_offset[gene, ]$k, c = df_offset[gene, ]$b) offset_exp_decay(x, a, k, c), color = "darkred", linetype = "dashed") + ylab("allelic ratio") + 
  ylim(c(0, 1)) + ggtitle(gene) + 
  annotate("text", label = paste0("R^2: ", round(all_coefs[all_coefs$gene == gene, ]$V3, digits = 4)), x = 2, y = 1, size = 8) + 
  annotate("text", label = paste0("t1/2: ", round(all_coefs[all_coefs$gene == gene, ]$halflifes, digits = 4)), x = 2, y = 0.9, size = 8)

# look at BICs between models: 
merged_df %>%
  ggplot(aes(x = BIC_native, BIC_offset, col = category)) + geom_point() + geom_abline() +
  ggrepel::geom_text_repel(aes(label = gene))

merged_df %>%
  ggplot(aes(x = halflifes_off, b_off, col = BIC_native > BIC_offset)) + geom_point() + xlim(c(0, 15)) +  ylim(c(0, 0.3)) + ggrepel::geom_text_repel(aes(label = gene)) + 
  theme_paper()

merged_df %>%
  ggplot(aes(x = halflifes_off, b_off, col = category)) + geom_point() + xlim(c(0, 15)) +  ylim(c(0, 0.3)) + ggrepel::geom_text_repel(aes(label = gene)) + 
  theme_paper() + scale_color_manual(values = escape_colors)


merged_df %>%
  dplyr::filter(!is.na(category)) %>%
  ggplot(aes(x = category, y = halflifes_off, fill = category)) + geom_boxplot(width = 0.1, outlier.color = NA) + ggbeeswarm::geom_quasirandom() + theme_paper() + 
  scale_fill_manual(values = escape_colors) + xlab("") + ylab("half-life (non-linear regression) [days]") + ggtitle("Escape half-life")

merged_df %>%
  dplyr::filter(!is.na(category)) %>%
  ggplot(aes(x = category, y = b_off, fill = category)) + geom_boxplot(width = 0.1, outlier.color = NA) + ggbeeswarm::geom_quasirandom() + theme_paper() + 
  scale_fill_manual(values = escape_colors) + xlab("") + ylab("offset parameter (non-linear regression)") + ggtitle("Residual escape")

merged_df %>%
  dplyr::filter(!is.na(category)) %>%
  ggplot(aes(x = category, y = halflifes_off, fill = category)) + geom_boxplot(width = 0.1, outlier.color = NA) + ggbeeswarm::geom_quasirandom() +theme_paper() + 
  scale_fill_manual(values = escape_colors) + xlab("") + ylab("half-life (non-linear regression) [days]") + ggtitle("Escape half-life") + 
  ggrepel::geom_text_repel(aes(label = gene))

merged_df %>%
  dplyr::filter(!is.na(category)) %>%
  ggplot(aes(x = category, y = b_off, fill = category)) + geom_boxplot(width = 0.1, outlier.color = NA) + ggbeeswarm::geom_quasirandom() + theme_paper() + 
  scale_fill_manual(values = escape_colors) + xlab("") + ylab("offset parameter (non-linear regression)") + ggtitle("Residual escape") + 
  ggrepel::geom_text_repel(aes(label = gene))

