### We now work towards figure-level plots

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

source("./Scripts/auxiliary.R")

# Read known information about escapee-status. 
# This is based on a literature survey and classifies genes into constitutive, facultative and non-escaping
# apart from Xist (which isn't really an escapee bc it's only from the Xi)
# escapee_annotation <- readxl::read_excel("~/Desktop/Projects/XChromosome_Project/ProcessedData/ListOfEscapeeFromEdithLab.xlsx") %>%
#   dplyr::select(c("ENSEMBL_v102", "final status"))
# library("EnsDb.Mmusculus.v79") 
# symbols <- mapIds(EnsDb.Mmusculus.v79, keys = escapee_annotation$ENSEMBL_v102, keytype = "GENEID", column="SYMBOL")
# escapee_annotation$symbol <- symbols
# convert_status_names <- setNames(c("constitutive", "facultative", "variable"), c("E", "V", "S"))
# convert_status_names2 <- setNames(c("constitutive", "facultative", "silenced / variable"), c("E", "V", "S"))
# escapee_annotation$our_status <- unlist(convert_status_names[escapee_annotation$`final status`])
# escapee_annotation$our_status2 <- unlist(convert_status_names2[escapee_annotation$`final status`])

# Now read our pre-processed dataset
data <- readRDS("./ProcessedData/merged_dataset.rds")
data <- computeSumFactors(data)

# For the first analysis, exclude the Xist-knockdown experiments
data <- data[,!data$Guide %in% c("g13", "P3") & data$Experiment != "May2022"]
#data <- data[,data$Experiment != "September2022"]

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
hist(log10((rowSums(counts_inactive(data)) + rowSums(counts_active(data))) / ncol(data) + 1), breaks = 100)
data <- data[(rowSums(counts_inactive(data)) + rowSums(counts_active(data))) / ncol(data) > 10, ]

# subset on X-linked genes
data <- data[seqnames(data) == "X", ]

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Figure 1
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Xist overexpression silences escapees in post-silencing contexts
# First, we analyze which genes escape in the baseline in our 3 different clones
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

data_baseline_x <- data[,data$Condition == "Aux_0_Dox_0_WO_0_WOAuxNO"]

# Plot escape in individual clones
sample_nr = 1 # Clone 30

data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle("Clone 30") + 
  xlab("Total expression (allelic reads only)") + ylab("ASE (d-score)")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/cl30_escape.pdf", width = 12, height = 8)

data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>% 
  ggplot(aes(total, ase, col = escape_status)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle("Clone 30") + 
  xlab("Total expression (allelic reads only)") + ylab("ASE (d-score)")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/cl30_escape_annotation.pdf", width = 12, height = 8)

data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>%
  dplyr::filter(!is.na(escape_status)) %>%
  mutate(escape_status = factor(escape_status, levels = c("silenced / variable", "facultative", "constitutive"))) %>%
  ggplot(aes(x = escape_status, y = ase, fill = escape_status)) + geom_boxplot() + theme_paper() + 
  xlab("") + ylab("ASE (d-score)") + scale_fill_manual(values = rev(scales::hue_pal()(3))) + ggtitle("") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/cl30_escape_annotation_boxplot.pdf", width = 6, height = 8)

data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>%
  dplyr::filter(!is.na(escape_status)) %>%
  mutate(escape_status = factor(escape_status, levels = c("silenced / variable", "facultative", "constitutive"))) %>%
  ggplot(aes(x = escape_status, y = ase, col = escape_status)) + ggbeeswarm::geom_quasirandom() + theme_paper() + 
  xlab("") + ylab("ASE (d-score)") + scale_fill_manual(values = rev(scales::hue_pal()(3))) + ggtitle("") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/cl30_escape_annotation_distribution.pdf", width = 6, height = 8)

# Make the same plots for E6: (replicate 1)
sample_nr = 3

data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle("Clone E6 (replicate1)") + 
  xlab("Total expression (allelic reads only)") + ylab("ASE (d-score)")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/e6_escape.pdf", width = 12, height = 8)

# Prettier plot, with escapee-annotation
data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>% 
  ggplot(aes(total, ase, col = escape_status)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") +  ggtitle("Clone E6 (replicate1)") + 
  xlab("Total expression (allelic reads only)") + ylab("ASE (d-score)") + coord_fixed(ratio = 3)
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/e6_escape_annotation.pdf", width = 12, height = 8)

data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>%
  dplyr::filter(!is.na(escape_status)) %>%
  mutate(escape_status = factor(escape_status, levels = c("silenced / variable", "facultative", "constitutive"))) %>%
  ggplot(aes(x = escape_status, y = ase, fill = escape_status)) + geom_boxplot() + theme_paper() + 
  xlab("") + ylab("ASE (d-score)") + scale_fill_manual(values = rev(scales::hue_pal()(3))) + ggtitle("") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/cle6_escape_annotation_boxplot.pdf", width = 6, height = 8)

data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>%
  dplyr::filter(!is.na(escape_status)) %>%
  mutate(escape_status = factor(escape_status, levels = c("silenced / variable", "facultative", "constitutive"))) %>%
  ggplot(aes(x = escape_status, y = ase, col = escape_status)) + ggbeeswarm::geom_quasirandom() + theme_paper() + 
  xlab("") + ylab("ASE (d-score)") + scale_fill_manual(values = rev(scales::hue_pal()(3))) + ggtitle("") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/cle6_escape_annotation_distribution.pdf", width = 6, height = 8)


### Now compare all the clones
data_here <- data[,data$ConditionClean == "Control"]

# get d-scores across genes + samples
ratios <- counts_inactive(data_baseline_x) / (counts_inactive(data_baseline_x) + counts_active(data_baseline_x))

# exclude genes with 0 counts in individual samples
rownames( ratios[!apply(ratios, 1, function(x){!any(is.na(x))}), ] )
ratios <- ratios[apply(ratios, 1, function(x){!any(is.na(x))}), ]

colnames(ratios) <- c("CL30", "CL31", "E6 (rep1)", "E6 (rep1)", "E6 (rep3)", "CL31.16", "CL30.7")

genes_show <- c("Xist", "Mecp2", "Kdm6a", "Kdm5c")
genes_of_choice <- setNames(rep("", nrow(ratios)), rownames(ratios))
genes_of_choice[genes_show] <- genes_show

# Plot d-scores for genes after ordering the genes by chromosomale coordinate
ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  dplyr::filter(Condition == "Control") %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  ggplot(aes(Clone, y = reorder(name, position), fill = value)) + geom_tile(width = 0.9) +
    theme_paper() + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
    theme(axis.ticks.y = element_blank()) + xlab("") + ylab("Chromosome position") + 
    scale_fill_gradientn(colors = c("#808000", "pink", "#008080")) + 
    scale_y_discrete(labels = genes_of_choice, position = "right") + ylab("")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/heatmap_ordered.pdf", width = 10, height = 10)

# Plot d-scores for genes that escape in any sample, excluding genes with > 0.7 d-score
ratios <- ratios[rownames(ratios) == "Xist" | (apply(ratios, 1, function(x){any(x > 0.1)}) & !apply(ratios, 1, function(x){any(x > 0.7)})), ]
ratios <- ratios[order(rowMeans(ratios), decreasing = T), ]

colnames(ratios) <- c("CL30", "CL31", "E6 (rep1)", "E6 (rep1)", "E6 (rep3)", "CL31.16", "CL30.7")
pdf("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/heatmap_onlyescape.pdf", width = 10, height = 10)
pheatmap::pheatmap(ratios, cluster_rows = F, cluster_cols = F, fontsize = 20, show_rownames =  F, 
                   color=colorRampPalette(c("#808000", "pink", "#008080"))(50))
dev.off()
dev.off()

# Now we define the baseline escapees per sample as the genes with d-score > 0.1
escapees <- apply(ratios, 2, function(x){x > 0.1})
escapees_per_clone <- apply(escapees, 2, function(x){rownames(escapees)[x]})
escapees_per_clone <- lapply(escapees_per_clone, function(x){x[x != "Xist"]})
names(escapees_per_clone) <- c("CL30", "CL31", "E6 (rep1)", "E6 (rep1)", "E6 (rep3)", "CL31.16", "CL30.7")

# How many do we get per clone and do they overlap?
ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  dplyr::filter(Condition == "Control") %>%
  dplyr::filter(value > 0.1) %>%
  group_by(name) %>%
  mutate(escape_in_clones = n() > 3) %>%
  ggplot(aes(x = Clone, fill = escape_in_clones)) + geom_bar() + theme_paper() + 
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) + xlab("") + ylab("Number of genes with d-score > 0.1") + 
    scale_fill_manual(values = c("grey", "darkred"), name = "Escape in all clones")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new//Fig1/number_of_escapees_per_clone.pdf", width = 8, height = 8)

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Figure 2b
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Now, we look at the effect of Xist-treatment on escapee expression
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### Now we look at Xist overexpression and how that silences escapees
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO")]
data_here <- data_here[,!data_here$Experiment %in% c("September2022", "October2022")]

ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))
ratios <- ratios[apply(ratios, 1, function(x){!any(is.na(x))}), ]

# For Xist-expression, we need a scaling factor per sample to account for differences in seq depth
size_factors <- colSums(counts(data_here)) / colSums(counts(data_here))[[1]]

# check that Xist-overexpression works
data.frame(
  Dox = data_here$ndDox,
  Clone = paste0(data_here$Clone, " - ", data_here$Experiment), 
  Expression = as.numeric(counts(data_here["Xist", ]) / colSums(counts(data_here)) * 1e6)
) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", paste0("Dox (", .$Dox, " days)"))) %>%
  ggplot(aes(x = Condition, y = Expression)) + geom_boxplot() + geom_point() + coord_flip() + 
    geom_line(aes(group = Clone, linetype = Clone, col = Clone)) + theme_paper() + 
    xlab("") + ylab("Xist Expression (normalized CPM)")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/xist_overexpression.pdf", width = 8, height = 8)

# check that Xist-overexpression works
# data.frame(
#   Dox = data_here$ndDox,
#   Washout = data_here$washout, 
#   Clone = paste0(data_here$Clone, " - ", data_here$Experiment), 
#   Expression = as.numeric(counts(data_here["Xist", ]) / (colSums(counts(data_here)) * 1e6))
# ) %>%
#   add_column(TimeSeries = ifelse(.$Dox == 0, list(c("3", "7")), .$Dox)) %>%
#   unnest(TimeSeries) %>%
#   mutate(TimeSeries = paste0(unlist(TimeSeries), " days dox treatment")) %>%
#   add_column(Condition = ifelse(.$Dox == 0, "Control", .$Washout)) %>%
#   dplyr::filter(Condition %in% c("Control", "Dox only")) %>%
#   ggplot(aes(x = Condition, y = Expression)) + geom_boxplot() + geom_point() + coord_flip() + 
#     geom_line(aes(group = Clone, linetype = Clone, col = Clone)) + theme_paper() + 
#     xlab("Washout") + ylab("Xist Expression (normalized CPM)") + 
#     scale_x_discrete(labels = c("Control", "Dox only", "Dox + washout"))
# ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures/Fig2_prelim/xist_overexpression.pdf", width = 8, height = 8)

# Now look at the global effect of Xist-overexpression on escape
ratios %>%
  t() %>% data.frame() %>%
  add_column("Dox" = data_here$ndDox) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Dox, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", paste0("Dox (", .$Dox, " days)"))) %>%
  ggplot(aes(x = Clone, y = value, fill = Condition)) + geom_boxplot(col = "grey") + 
    theme_paper() + ylab("ASE (d-score)") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = c("black", "red", "blue", "darkred", "darkblue"), labels = c("Control", "Dox (3d)", "Dox (7d)"))
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/escapee_silencing.pdf", width = 8, height = 8)

# First, for CL30
ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  dplyr::filter(Condition %in% c("Control", "Dox (7d)", "Dox (7d) - washout") & Clone == "CL30 - December2021") %>%
  ggplot(aes(x = Condition, y = value, col = Condition, group = Condition)) + geom_point() + 
  geom_line(aes(group = name), col = "grey", alpha = 0.5) + theme_paper() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_manual(values = c("black", "red", "blue")) + xlab("") + ylab("ASE (d-score)") + 
  ggtitle("Reversibility after 7d (CL30)")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/gene_level_silencing_cl30.pdf", width = 8, height = 8)

ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  #dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  dplyr::filter(Condition %in% c("Control", "Dox (7d)", "Dox (7d) - washout") & Clone == "CL30 - December2021") %>%
  ggplot(aes(x = Condition, y = value, col = Condition, group = Condition)) + geom_point() + 
  geom_line(aes(group = name), col = "grey", alpha = 0.5) + theme_paper() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_manual(values = c("black", "red", "blue")) + xlab("") + ylab("ASE (d-score)") + 
  ggtitle("Reversibility after 7d (CL30)")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/gene_level_silencing_cl30_withallgenes.pdf", width = 8, height = 8)

# Second, for E6
ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  dplyr::filter(Condition %in% c("Control", "Dox (7d)", "Dox (7d) - washout") & Clone == "E6 - October2021") %>%
  ggplot(aes(x = Condition, y = value, col = Condition, group = Condition)) + geom_point() + 
  geom_line(aes(group = name), col = "grey", alpha = 0.5) + theme_paper() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_manual(values = c("black", "red", "blue")) + xlab("") + ylab("ASE (d-score)") + 
  ggtitle("Reversibility after 7d (E6)")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/gene_level_silencing_e6.pdf", width = 8, height = 8)

ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  #dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  dplyr::filter(Condition %in% c("Control", "Dox (7d)", "Dox (7d) - washout") & Clone == "E6 - October2021") %>%
  ggplot(aes(x = Condition, y = value, col = Condition, group = Condition)) + geom_point() + 
  geom_line(aes(group = name), col = "grey", alpha = 0.5) + theme_paper() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_manual(values = c("black", "red", "blue")) + xlab("") + ylab("ASE (d-score)") + 
  ggtitle("Reversibility after 7d (E6)")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/gene_level_silencing_e6_withallgenes.pdf", width = 8, height = 8)

# Now, check how silencing occurs across the entire chromosome

### Now compare all the clones
data_here <- data[,data$ConditionClean %in% c("Control", "Dox (3d)")]
data_here <- data_here[,!data_here$Experiment %in% c("September2022", "October2022")]

# get d-scores across genes + samples
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

# exclude genes with 0 counts in individual samples
rownames( ratios[!apply(ratios, 1, function(x){!any(is.na(x))}), ] )
ratios <- ratios[apply(ratios, 1, function(x){!any(is.na(x))}), ]

colnames(ratios) <- c("CL30", "CL31", "E6 (rep1)", "E6 (rep2)")

genes_show <- c("Xist", "Mecp2", "Kdm6a", "Kdm5c")
genes_of_choice <- setNames(rep("", nrow(ratios)), rownames(ratios))
genes_of_choice[genes_show] <- genes_show

# Plot d-scores for genes after ordering the genes by chromosomale coordinate
ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  #dplyr::filter(Condition == "Control") %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  ggplot(aes(Condition, y = reorder(name, position), fill = value)) + geom_tile(width = 0.9) +
    theme_paper() + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
    theme(axis.ticks.y = element_blank()) + xlab("") + ylab("Chromosome position") + 
    scale_fill_gradientn(colors = c("#808000", "pink", "#008080")) + 
    scale_y_discrete(labels = genes_of_choice, position = "right") + ylab("") + 
    facet_wrap(~Clone, nrow = 1)
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/heatmap_ordered_silencing.pdf", width = 10, height = 10)

# Plot d-scores for genes after ordering the genes by chromosomale coordinate
ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  #dplyr::filter(Condition == "Control") %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  dplyr::filter(Clone == "E6 - March2022") %>%
  ggplot(aes(Condition, y = reorder(name, position), fill = value)) + geom_tile(width = 0.9) +
    theme_paper() + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
    theme(axis.ticks.y = element_blank()) + xlab("") + ylab("Chromosome position") + 
    scale_fill_gradientn(colors = c("#808000", "pink", "#008080")) + 
    scale_y_discrete(labels = genes_of_choice, position = "right") + ylab("") + 
    facet_wrap(~Clone, nrow = 1)
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig1/heatmap_ordered_silencing_reduced_onlyE6march2022.pdf", width = 5, height = 10)

### In the actual paper, the SPEN part comes here... 


### This is for Figure 3, the reversibility aspect
### Now we look at Xist overexpression and how that silences escapees
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_4_WOAuxNO",
                                         "Aux_0_Dox_7_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_4_WOAuxNO")]
data_here <- data_here[,!data_here$Experiment %in% c("September2022", "October2022")]

ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))
ratios <- ratios[apply(ratios, 1, function(x){!any(is.na(x))}), ]

# For Xist-expression, we need a scaling factor per sample to account for differences in seq depth
size_factors <- colSums(counts(data_here)) / colSums(counts(data_here))[[1]]

# check that Xist-overexpression works
data.frame(
  Dox = data_here$ndDox,
  Washout = data_here$washout, 
  Clone = paste0(data_here$Clone, " - ", data_here$Experiment), 
  Expression = as.numeric(counts(data_here["Xist", ]) / colSums(counts(data_here)) * 1e6)
) %>%
  add_column(TimeSeries = ifelse(.$Dox == 0, list(c("3", "7")), .$Dox)) %>%
  unnest(TimeSeries) %>%
  mutate(TimeSeries = paste0(unlist(TimeSeries), " days dox treatment")) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", .$Washout)) %>%
  #dplyr::filter(Condition %in% c("Control", "Dox only")) %>% 
  ggplot(aes(x = Condition, y = Expression)) + geom_boxplot() + geom_point() + coord_flip() + 
    geom_line(aes(group = Clone, linetype = Clone, col = Clone)) + theme_paper() + 
    xlab("") + ylab("Xist Expression (normalized CPM)") + facet_wrap(~TimeSeries, nrow = 2) + 
    scale_x_discrete(labels = c("Control", "Dox", "Washout"))
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/xist_overexpression_reversibility.pdf", width = 8, height = 8)

# Now look at the global effect of Xist-overexpression on escape
ratios %>%
  t() %>% data.frame() %>%
  add_column("Dox" = data_here$ndDox) %>%
  add_column("Washout" = data_here$washout) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Dox, Washout, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  add_column(TimeSeries = ifelse(.$Dox == 0, list(c("3", "7")), .$Dox)) %>%
  unnest(TimeSeries) %>%
  mutate(TimeSeries = paste0(unlist(TimeSeries), " days dox treatment")) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", .$Washout)) %>%
  ggplot(aes(x = Clone, y = value, fill = Condition)) + geom_boxplot(col = "grey") + 
    theme_paper() + ylab("ASE (d-score)") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = c("black", "red", "blue", "darkred", "darkblue"), labels = c("Control", "Dox (3d)", "Dox (7d)")) + 
    facet_wrap(~TimeSeries)
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/escapee_silencing_reversible.pdf", width = 8, height = 8)

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Figure check washout after different days
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# Now, check how silencing occurs across the entire chromosome

### Now compare all the clones
data_here <- data[,data$Clone == "E6", data$ndAux == "0"]
#data_here <- data_here[,!data_here$Experiment %in% c("September2022", "October2022")]

# get d-scores across genes + samples
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

# exclude genes with 0 counts in individual samples
rownames( ratios[!apply(ratios, 1, function(x){!any(is.na(x))}), ] )
ratios <- ratios[apply(ratios, 1, function(x){!any(is.na(x))}), ]

ratios %>%
  t() %>% data.frame() %>%
  add_column("Dox" = data_here$ndDox) %>%
  add_column("Washout" = factor(data_here$ndWashout)) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Dox, Clone, Washout)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", paste0("Dox (", .$Dox, " days)"))) %>%
  mutate(Washout = factor(Washout, levels = c("0", "4", "7", "14"))) %>%
  mutate(Condition = factor(Condition, levels = c("Control", "Dox (3 days)", "Dox (7 days)", "Dox (14 days)", "Dox (21 days)"))) %>%
  ggplot(aes(x = Condition, y = value, fill = Washout)) + 
    geom_boxplot(col = "grey", position = position_dodge(preserve = "single")) + 
    theme_paper() + ylab("ASE (d-score)") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 1))

### In the actual paper, the SPEN part comes here... 

### This is for Figure 3, the reversibility aspect
### Now we look at Xist overexpression and how that silences escapees
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_4_WOAuxNO",
                                         "Aux_0_Dox_7_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_4_WOAuxNO")]
data_here <- data_here[,!data_here$Experiment %in% c("September2022", "October2022")]

ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))
ratios <- ratios[apply(ratios, 1, function(x){!any(is.na(x))}), ]

# For Xist-expression, we need a scaling factor per sample to account for differences in seq depth
size_factors <- colSums(counts(data_here)) / colSums(counts(data_here))[[1]]

# check that Xist-overexpression works
data.frame(
  Dox = data_here$ndDox,
  Washout = data_here$washout, 
  Clone = paste0(data_here$Clone, " - ", data_here$Experiment), 
  Expression = as.numeric(counts(data_here["Xist", ]) / colSums(counts(data_here)) * 1e6)
) %>%
  add_column(TimeSeries = ifelse(.$Dox == 0, list(c("3", "7")), .$Dox)) %>%
  unnest(TimeSeries) %>%
  mutate(TimeSeries = paste0(unlist(TimeSeries), " days dox treatment")) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", .$Washout)) %>%
  #dplyr::filter(Condition %in% c("Control", "Dox only")) %>% 
  ggplot(aes(x = Condition, y = Expression)) + geom_boxplot() + geom_point() + coord_flip() + 
  geom_line(aes(group = Clone, linetype = Clone, col = Clone)) + theme_paper() + 
  xlab("") + ylab("Xist Expression (normalized CPM)") + facet_wrap(~TimeSeries, nrow = 2) + 
  scale_x_discrete(labels = c("Control", "Dox", "Washout"))
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/xist_overexpression_reversibility.pdf", width = 8, height = 8)

# Now look at the global effect of Xist-overexpression on escape
ratios %>%
  t() %>% data.frame() %>%
  add_column("Dox" = data_here$ndDox) %>%
  add_column("Washout" = data_here$washout) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Dox, Washout, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  add_column(TimeSeries = ifelse(.$Dox == 0, list(c("3", "7")), .$Dox)) %>%
  unnest(TimeSeries) %>%
  mutate(TimeSeries = paste0(unlist(TimeSeries), " days dox treatment")) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", .$Washout)) %>%
  ggplot(aes(x = Clone, y = value, fill = Condition)) + geom_boxplot(col = "grey") + 
  theme_paper() + ylab("ASE (d-score)") + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("black", "red", "blue", "darkred", "darkblue"), labels = c("Control", "Dox (3d)", "Dox (7d)")) + 
  facet_wrap(~TimeSeries)
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/escapee_silencing_reversible.pdf", width = 8, height = 8)

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Figure 2c
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# We now use gLMs to test individual genes for Xist + washout effects
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

library(VGAM)

data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_4_WOAuxNO", 
                                         "Aux_0_Dox_7_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_4_WOAuxNO")]
data_here <- data_here[,!data_here$Experiment %in% c("September2022", "October2022")]

# Note that we define the escapees slightly differently here, harmonize at some point (0.1 / 0.05)
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))
escapees <- apply(ratios, 2, function(x){x > 0.05 & x < 0.7})
escapees[is.na(escapees)] <- FALSE

# Define genes that escape in all clones
general_escapees <- data_here[names(which(apply(escapees[,data_here$Condition == "Aux_0_Dox_0_WO_0_WOAuxNO"], 1, function(x){all(x)}))), ]

# Get count matrices etc
data_test_inactive <- counts_inactive(general_escapees)
data_test_active <- counts_active(general_escapees)

fitting_metadata_here <- colData(general_escapees)

# The model will be d ~ b0 + b1 * Dox * b2 * Washout
# We need to add a column with washout-day to fit independent coefficients to that
fitting_metadata_here <- fitting_metadata_here %>%
  data.frame() %>%
  add_column(washout_day = paste0(.$washout, "_", .$ndDox)) %>%
  mutate(washout_day = ifelse(grepl("0", washout_day), '0', .$washout_day))

all_coefs <- data.frame(do.call("rbind", lapply(1:nrow(data_test_active), function(i){
  tryCatch({
    y = as.numeric(data_test_inactive[i, ]) + 1
    N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ]) + 2
    fit1 <- glm(cbind(y, N-y) ~ Clone + ndDox + washout_day, family = "binomial", data = fitting_metadata_here)
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
all_coefs <- all_coefs[!apply(all_coefs, 1, function(x){all(is.na(x))}), !apply(all_coefs, 2, function(x){all(is.na(x))})]
all_coefs[is.na(all_coefs)] <- 0
#all_coefs[abs(all_coefs) > 5] <- sign(all_coefs[abs(all_coefs) > 5]) * 5

row_annotation <- data.frame(
  Gene = rownames(all_coefs), 
  escapee_annotation = escapee_annotation[match(rownames(all_coefs), escapee_annotation$symbol), ][,2]
) %>% column_to_rownames("Gene")

pheatmap::pheatmap(all_coefs, cluster_cols = F, color=colorRampPalette(c("navy", "white", "red"))(50), annotation_row = row_annotation, fontsize = 7)

all_coefs %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-c(Gene)) %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status)) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(x = name, y = value, col = escape_status)) + 
  geom_boxplot() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.2) +
  coord_flip() + theme_paper() + xlab("Covariate") +
  theme(aspect.ratio = 2) + theme(legend.position = "top") + geom_hline(linetype = "dashed", yintercept = 0)
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/All_model_coefficients.pdf")

# look at coefficients between washout and dox silencing -- looks cool! 
# To do: Find a model that does the LM and also considers the (known) uncertainties in x and y
all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status)) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(ndDox3, ndWashout, col = escape_status)) + geom_point() + 
  geom_smooth(method = "lm", col = "black", linetype = 'dashed', size = 0.5) + 
  ggrepel::geom_text_repel(aes(label = Gene)) + 
  geom_abline(linetype = "dashed", slope = -1) + 
  theme_paper() + xlab("Dox-effect (LM-coefficient)") + ylab("Washout-effect (LM-coefficient)") + 
  ggtitle("3d analysis")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/LM_coefs_3d_allClones.pdf", width = 8, height = 8)

all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status)) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(ndDox7, washout_dayYES_7, col = escape_status)) + geom_point() + 
  geom_smooth(method = "lm", col = "black", linetype = 'dashed', size = 0.5) + 
  ggrepel::geom_text_repel(aes(label = Gene)) + 
  geom_abline(linetype = "dashed", slope = -1) + 
  theme_paper() + xlab("Dox-effect (LM-coefficient)") + ylab("Washout-effect (LM-coefficient)") + 
  ggtitle("7d analysis")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/LM_coefs_7d_allClones.pdf", width = 8, height = 8)

# Now we plot example genes to illustrate the diversity of responses
plot_gene <- function(gene, save = F, path = NULL, suffix = "_escape_plot.pdf"){
  i = gene
  y = as.numeric(data_test_inactive[i, ])
  N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
  # fit1 <- glm(cbind(y, N-y) ~ Clone + ndDox + washout, family = "binomial", data = fitting_metadata_here)
  # fit0 <- glm(cbind(y, N-y) ~ 1, family = "binomial", data = fitting_metadata_here)
  
  pp <- data.frame(
    y = y,
    N = N,
    Dox = data_here$ndDox,
    Washout = data_here$washout, 
    Clone = paste0(data_here$Clone, " - ", data_here$Experiment)
  ) %>%
    add_column(TimeSeries = ifelse(.$Dox == 0, list(c("3", "7")), .$Dox)) %>%
    unnest(TimeSeries) %>%
    mutate(TimeSeries = paste0(unlist(TimeSeries), " days dox treatment")) %>%
    add_column(Condition = ifelse(.$Dox == 0, "Control", .$Washout)) %>%
    mutate(Condition = factor(Condition, levels = rev(c("Control", "NO", "YES")))) %>%
    ggplot(aes(x = Condition, y = y / N, col = Clone)) + geom_point(size = 4) + coord_flip() + ylim(0, 1) +
    geom_line(aes(group = Clone), linetype = 'dashed') + theme_paper() + facet_wrap(~TimeSeries, nrow = 2) + 
    xlab("") + ylab("ASE (d-score)") + 
    scale_x_discrete(labels = c("Dox + washout", "Dox only", "Control"))
  
  if (save){
    ggsave(paste0(path, "./", gene, suffix), pp)
  } 
  
  return(pp)
}

plot_gene("Fmr1", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")
plot_gene("Plxnb3", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")
plot_gene("Jpx", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")
plot_gene("Eif2s3x", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")
plot_gene("Kdm5c", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")
plot_gene("Ddx3x", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")

plot_gene("Slc6a8", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")
plot_gene("Bcap31", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")
plot_gene("Idh3g", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")
plot_gene("Ssr4", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")
plot_gene("Hcfc1", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")
plot_gene("Mecp2", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")
plot_gene("G6pdx", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")
plot_gene("Jpx", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/", "_escape_all_clones.pdf")

head(all_coefs[order(all_coefs$ndDox7), ], n = 20)


## Verify that change in allelic balance comes from increase in Xi, not decrease in Xa: 

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
    add_column(TimeSeries = ifelse(.$Dox == 0, list(c("3", "7")), .$Dox)) %>%
    unnest(TimeSeries) %>%
    mutate(TimeSeries = paste0(unlist(TimeSeries), " days dox treatment")) %>%
    add_column(Condition = ifelse(.$Dox == 0, "Control", .$Washout)) %>%
    mutate(Condition = factor(Condition, levels = rev(c("Control", "NO", "YES")))) %>%
    pivot_longer(cols = c("inactive", "active")) %>%
    mutate(group = paste0(Clone, "_", name)) %>%
    ggplot(aes(x = Condition, y = value, col = Clone, shape = name)) + 
      geom_point(size = 4) + 
      geom_line(aes(group = group)) + 
      coord_flip() + 
      theme_paper() + facet_wrap(~TimeSeries, nrow = 2) + 
      xlab("") + ylab("ASE (d-score)") + 
      scale_x_discrete(labels = c("Dox + washout", "Dox only", "Control")) + 
      scale_y_log10()
  
  if (save){
    ggsave(paste0(path, "./", gene, suffix), pp)
  }
  
  return(pp)
}

plot_gene_total("Jpx")

## systematic analysis: DESeq on both chromosomes seperately, block clone and test e.g. over time

## 

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Figure 3 something
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Now, we do the LM fits, but only for E6, as it gives a larger space of escapees to explore
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
##################################################################################################
general_escapees <- data_here[names(which(apply(
  escapees[,data_here$Condition == "Aux_0_Dox_0_WO_NO_WOAuxNO" & data_here$Clone == "E6"], 
  1, function(x){all(x)}))), ]

data_test_inactive <- counts_inactive(general_escapees[,data_here$Clone == "E6"])
data_test_active <- counts_active(general_escapees[,data_here$Clone == "E6"])

fitting_metadata_here <- colData(general_escapees)[data_here$Clone == "E6", ]

#genes_out <- c("Car5b", "Apoo", "2210013O21Rik)
genes_out <- c("")
data_test_inactive <- data_test_inactive[!rownames(data_test_inactive) %in% genes_out, ]
data_test_active <- data_test_active[!rownames(data_test_active) %in% genes_out, ]

# The model will be d ~ b0 + b1 * Dox * b2 * Washout
# We need to add a column with washout-day to fit independent coefficients to that
fitting_metadata_here <- fitting_metadata_here %>%
  data.frame() %>%
  add_column(washout_day = paste0(.$washout, "_", .$ndDox)) %>%
  mutate(washout_day = ifelse(grepl("NO", washout_day), 'NO', .$washout_day))

all_coefs <- data.frame(do.call("rbind", lapply(1:nrow(data_test_active), function(i){
  tryCatch({
    y = as.numeric(data_test_inactive[i, ]) + 1
    N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ]) + 2
    fit1 <- glm(cbind(y, N-y) ~ ndDox + washout_day, family = "binomial", data = fitting_metadata_here)
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
all_coefs <- all_coefs[!apply(all_coefs, 1, function(x){all(is.na(x))}), !apply(all_coefs, 2, function(x){all(is.na(x))})]
all_coefs[is.na(all_coefs)] <- 0
#all_coefs[abs(all_coefs) > 5] <- sign(all_coefs[abs(all_coefs) > 5]) * 5
all_coefs_silencing <- all_coefs

# save coefficients for downstream analysis
saveRDS(all_coefs, "./ProcessedData/LM_coefficients_E6.rds")

row_annotation <- data.frame(
  Gene = rownames(all_coefs), 
  escapee_annotation = escapee_annotation[match(rownames(all_coefs), escapee_annotation$symbol), ][,2]
) %>% column_to_rownames("Gene")

pheatmap::pheatmap(all_coefs, cluster_cols = F, color=colorRampPalette(c("navy", "white", "red"))(50), annotation_row = row_annotation, fontsize = 7)

all_coefs %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-c(Gene)) %>%
  ggplot(aes(x = name, y = value)) + geom_boxplot() +
  coord_flip() + theme_paper() + xlab("Covariate") +
  theme(aspect.ratio = 2)

all_coefs %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-c(Gene)) %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status)) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(x = name, y = value, col = escape_status)) + geom_boxplot() + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.2) +
    coord_flip() + theme_paper() + xlab("Covariate") +
    theme(aspect.ratio = 2) + theme(legend.position = "top") + geom_hline(yintercept = 0, linetype = "dashed") + 
    ylab("Model coefficient") + xlab("") + scale_x_discrete(labels = rev(c("Intercept", "Washout (3d)", "Washout (7d)", "7d dox", "3d dox")))
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/E6_only_model_coefficients.pdf")

## Compare Intercept vs Dox effects
# If we wanna check consistency between average baseline-escape and estimated intercepts
# plot(logit(rowMeans(data_test_inactive / (data_test_inactive + data_test_active))), all_coefs$X.Intercept.)
all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status)) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(X.Intercept., washout_dayYES_3, col = escape_status)) + geom_point() + 
  geom_smooth(method = "lm", col = "black", linetype = 'dashed', size = 0.5) + 
  ggrepel::geom_text_repel(aes(label = Gene)) + 
  theme_paper() + xlab("Basal effect (X-intercept)") + ylab("Washout-effect (LM-coefficient)") + 
  ggtitle("3d analysis, intercept vs dox")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/LM_coefs_3d_E6_intercept_vs_dox.pdf")

all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status)) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(X.Intercept., washout_dayYES_7, col = escape_status)) + geom_point() + 
  geom_smooth(method = "lm", col = "black", linetype = 'dashed', size = 0.5) + 
  ggrepel::geom_text_repel(aes(label = Gene)) + 
  theme_paper() + xlab("Basal effect (X-intercept)") + ylab("Washout-effect (LM-coefficient)") + 
  ggtitle("7d analysis, intercept vs dox")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/LM_coefs_7d_E6_intercept_vs_dox.pdf")

## Compare Dox vs Washout effects
all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status)) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(ndDox3, washout_dayYES_3, col = escape_status)) + geom_point() + 
  geom_smooth(method = "lm", col = "black", linetype = 'dashed', size = 0.5) + 
  ggrepel::geom_text_repel(aes(label = Gene)) + 
  geom_abline(linetype = "dashed", slope = -1) + 
  theme_paper() + xlab("Dox-effect (LM-coefficient)") + ylab("Washout-effect (LM-coefficient)") + 
  ggtitle("3d analysis")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig3/LM_coefs_3d_E6.pdf")

all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status)) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(ndDox7, washout_dayYES_7, col = escape_status)) + geom_point() + 
  geom_smooth(method = "lm", col = "black", linetype = 'dashed', size = 0.5) + 
  ggrepel::geom_text_repel(aes(label = Gene)) + 
  geom_abline(linetype = "dashed", slope = -1) + 
  theme_paper() + xlab("Dox-effect (LM-coefficient)") + ylab("Washout-effect (LM-coefficient)") + 
  ggtitle("7d analysis")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures/Fig2_prelim/LM_coefs_7d_E6.pdf")

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### Auxin Treatment
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

### Figure 3: 
# Xist silencing in post-silencing contexts is spen-dependent
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_NO_WOAuxNO", "Aux_0_Dox_3_WO_NO_WOAuxNO", "Aux_0_Dox_7_WO_NO_WOAuxNO", 
                                         "Aux_2_Dox_0_WO_NO_WOAuxNO", "Aux_5_Dox_0_WO_NO_WOAuxNO", "Aux_9_Dox_0_WO_NO_WOAuxNO", 
                                         "Aux_5_Dox_3_WO_NO_WOAuxNO", "Aux_9_Dox_7_WO_NO_WOAuxNO")]

data_here <- data_here[,!data_here$Clone == "E6"]

# check that xist-overexpression works

rename_conditions <- c(
  "NA_NA" = "Control", 
  "Aux_NA" = "- Dox / + Aux", 
  "NA_Dox" = "+ Dox / - Aux", 
  "Aux_Dox" = "+ Dox / + Aux"
)

data.frame(
  Condition = factor(data_here$ConditionClean, levels = c("Control", "Aux (2d)", "Aux (5d)", "Aux (9d)", "Dox (3d)", 
                                                          "Dox (3d), Aux (5d)", "Dox (7d)", "Dox (7d), Aux (9d)")), 
  Clone = paste0(data_here$Clone, " - ", data_here$Experiment), 
  Expression = as.numeric(counts(data_here["Xist", ]) / colSums(counts(data_here)) * 1e6), 
  Dox = data_here$ndDox, 
  Aux = pmax(as.numeric(data_here$ndAux) - 2, 0)
) %>%
  add_column(Timepoint = pmax(.$Dox, .$Aux)) %>%
  add_column(Condition2 = paste0(str_extract(.$Condition, "Aux"), "_", str_extract(.$Condition, "Dox"))) %>%
  add_column(TimeSeries = ifelse(.$Dox == 0, list(c("3", "7")), .$Dox)) %>%
  unnest(TimeSeries) %>%
  mutate(Timepoint = paste0(Timepoint, " days")) %>%
  mutate(Condition2 = factor(rename_conditions[.$Condition2], levels = rev(rename_conditions))) %>%
  #dplyr::filter(Timepoint == "7 days" & Condition == "NA_NA") %>%
  ggplot(aes(x = Condition2, y = Expression, col = Clone)) + geom_point() + 
  geom_line(aes(group = Clone), linetype = "dashed") + 
  coord_flip() + theme_paper() + 
    xlab("") + ggtitle("Xist Expression (normalized CPM)") + 
  facet_wrap(~Timepoint, ncol = 1, scales = "free_y") + 
  theme_paper()
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig2/Xist_try2.pdf")

# check changes in escapee expression
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

#escapees_per_clone <- apply(escapees, 2, function(x){rownames(escapees)[x]})
escapees_per_clone <- lapply(escapees_per_clone, function(x){x[x != "Xist"]})
names(escapees_per_clone) <- c("CL30 - December2021", "CL31 - December2021", "E6 - October2021", "E6 - March2022")

# new plot into timepoints
ratios %>%
  t() %>% data.frame() %>%
  add_column(Dox = data_here$ndDox) %>%
  add_column(Aux = pmax(as.numeric(data_here$ndAux) - 2, 0)) %>%
  add_column(Clone = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  add_column(Condition = data_here$ConditionClean) %>%
  pivot_longer(-c(Dox, Aux, Clone, Condition)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  ungroup() %>%
  add_column(Timepoint = pmax(.$Dox, .$Aux)) %>%
  add_column(Condition2 = paste0(str_extract(.$Condition, "Aux"), "_", str_extract(.$Condition, "Dox"))) %>%
  add_column(TimeSeries = ifelse(.$Dox == 0, list(c("3", "7")), .$Dox)) %>%
  unnest(TimeSeries) %>%
  mutate(Timepoint = paste0(Timepoint, " days")) %>%
  mutate(Condition2 = factor(rename_conditions[.$Condition2], levels = rename_conditions)) %>%
  ggplot(aes(x = Condition2, y = value, fill = Clone)) + geom_boxplot(col = "grey") + 
    theme_paper() + ylab("ASE (d-score)") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
    facet_wrap(~Timepoint, scale = "free_x") + 
  xlab("") + theme_paper()
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig2/escape_total_new.pdf")

#### LM stuff for SPEN-experiment
general_escapees <- data_here[intersect(escapees_per_clone$`CL30 - December2021`, escapees_per_clone$`CL31 - December2021`), ]

data_test_inactive <- counts_inactive(general_escapees)
data_test_active <- counts_active(general_escapees)

fitting_metadata_here <- colData(general_escapees)

genes_out <- names(which(rowMeans(data_test_inactive / (data_test_active + data_test_inactive)) > 0.6))
genes_out <- c(genes_out, "Rpl3-ps1")
data_test_inactive <- data_test_inactive[!rownames(data_test_inactive) %in% genes_out, ]
data_test_active <- data_test_active[!rownames(data_test_active) %in% genes_out, ]

all_coefs <- data.frame(do.call("rbind", lapply(1:nrow(data_test_active), function(i){
  tryCatch({
    y = as.numeric(data_test_inactive[i, ]) + 1
    N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ]) + 2
    fit1 <- glm(cbind(y, N-y) ~ ndDox * ndAux, family = "binomial", data = fitting_metadata_here)
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
all_coefs <- all_coefs[!apply(all_coefs, 1, function(x){all(is.na(x))}), !apply(all_coefs, 2, function(x){all(is.na(x))})]
all_coefs[is.na(all_coefs)] <- 0
#all_coefs[abs(all_coefs) > 5] <- sign(all_coefs[abs(all_coefs) > 5]) * 5

# save coefficients for downstream analysis
#saveRDS(all_coefs, "./ProcessedData/LM_coefficients_E6.rds")

row_annotation <- data.frame(
  Gene = rownames(all_coefs), 
  escapee_annotation = escapee_annotation[match(rownames(all_coefs), escapee_annotation$symbol), ][,2]
) %>% column_to_rownames("Gene")

pheatmap::pheatmap(all_coefs, cluster_cols = F, color=colorRampPalette(c("navy", "white", "red"))(50), annotation_row = row_annotation, fontsize = 7)

all_coefs %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-c(Gene)) %>%
  ggplot(aes(x = name, y = value)) + geom_boxplot() +
  coord_flip() + theme_paper() + xlab("Covariate") +
  theme(aspect.ratio = 2)

all_coefs %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-c(Gene)) %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status)) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(x = name, y = value, col = escape_status)) + geom_boxplot() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.2) +
  coord_flip() + theme_paper() + xlab("Covariate") +
  theme(aspect.ratio = 2) + theme(legend.position = "top") + geom_hline(yintercept = 0, linetype = "dashed") + 
  ylab("Model coefficient") + xlab("") #+ scale_x_discrete(labels = rev(c("Intercept", "Washout", "7d dox", "3d dox")))
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig2/spen_model_coefs.pdf")

all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status)) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(ndDox3, ndDox3.ndAux5, col = escape_status)) + geom_point() + 
  geom_smooth(method = "lm", col = "black", linetype = 'dashed', size = 0.5) + 
  ggrepel::geom_text_repel(aes(label = Gene)) + 
  geom_abline(linetype = "dashed", slope = -1) + 
  theme_paper() + xlab("Dox-effect (LM-coefficient)") + ylab("SPEN-effect (LM-coefficient)") + 
  ggtitle("3d analysis")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig2/spen_vs_dox_model_3d.pdf")

all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status)) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(ndDox7, ndDox7.ndAux9, col = escape_status)) + geom_point() + 
  geom_smooth(method = "lm", col = "black", linetype = 'dashed', size = 0.5) + 
  ggrepel::geom_text_repel(aes(label = Gene)) + 
  geom_abline(linetype = "dashed", slope = -1) + 
  theme_paper() + xlab("Dox-effect (LM-coefficient)") + ylab("SPEN-effect (LM-coefficient)") + 
  ggtitle("7d analysis")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig2/spen_vs_dox_model_7d.pdf")

plot_gene <- function(gene, save = F, path = NULL, suffix = "_escape_plot.pdf"){
  i = gene
  y = as.numeric(data_test_inactive[i, ])
  N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
  fit1 <- glm(cbind(y, N-y) ~ Clone + ndDox + ndAux, family = "binomial", data = fitting_metadata_here)
  fit0 <- glm(cbind(y, N-y) ~ 1, family = "binomial", data = fitting_metadata_here)
  
  pp <- data.frame(
    y = y,
    N = N,
    Clone = paste0(data_here$Clone, " - ", data_here$Experiment),
    Cov = data_here$ConditionClean
  ) %>% ggplot(aes(Cov, y / N, col = Clone)) + geom_point(size = 4) + coord_flip() + ylim(0, 1) +
    ggtitle(i) +
    geom_hline(yintercept = rev_logit(coef(fit1)[[1]]), linetype = "dashed") + 
    geom_line(aes(group = Clone), linetype = "dashed") + 
    theme_paper() + xlab("") + ylab("ASE (d-score)")

  if (save){
    ggsave(paste0(path, "./", gene, suffix), pp)
  } 
  
  return(pp)
}

plot_gene("Kdm6a", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig2/", "_escape_all_clonesAux.pdf")
plot_gene("Fmr1", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig2/", "_escape_all_clonesAux.pdf")
plot_gene("Plxnb3", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig2/", "_escape_all_clonesAux.pdf")
plot_gene("Ddx3x", save = T, "~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig2/", "_escape_all_clonesAux.pdf")

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### Knockdown
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

### Figure 4: 

data <- readRDS("./ProcessedData/merged_dataset.rds")
data_here <- data[,data$Guide %in% c("g13", "P3")]
data_here <- data_here[as.character(seqnames(data_here)) == "X", ]
data_here$ConditionClean <- ifelse(data_here$Guide == "P3", "non-targeting control", "sgXist")

#data_here <- data_here[rowSums(counts_inactive(data_here) + counts_active(data_here)) > 200, ]

# check that xist-knockdown works
data.frame(
  Condition = data_here$ConditionClean, 
  Clone = paste0(data_here$Clone, " - ", data_here$Experiment), 
  Xist_Expression = as.numeric(counts(data_here["Xist", ]) / colSums(counts(data_here)) * 1e6)
) %>%
  ggplot(aes(x = factor(Condition, levels = c("sgXist", "non-targeting control")),  y = Xist_Expression)) + geom_boxplot() + geom_point() + 
  geom_line(aes(group = Clone, linetype = Clone, col = Clone)) + theme_paper() + 
  coord_flip() + xlab("") + ylab("Xist Expression (normalized CPM)")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig5/xist.pdf")

data_baseline_x <- data_here

data_baseline_x <- data_here[,data_here$Guide == "P3"]

ratios <- counts_inactive(data_baseline_x) / (counts_inactive(data_baseline_x) + counts_active(data_baseline_x))
escapees_per_clone <- apply(escapees, 2, function(x){rownames(escapees)[x]})
escapees_per_clone <- lapply(escapees_per_clone, function(x){x[x != "Xist" & !is.na(x)]})
names(escapees_per_clone) <- c("C3March2022", "D11_March2022")

ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  #dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  ggplot(aes(x = Clone, y = value, fill = Condition)) + geom_boxplot(col = "grey") + 
  theme_paper() + ylab("ASE (d-score)") + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("black", "red")) + xlab("")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig5/escape_all_genes.pdf")

ratios <- counts_inactive(data_baseline_x) / (counts_inactive(data_baseline_x) + counts_active(data_baseline_x))
escapees <- apply(ratios, 2, function(x){x > 0})
escapees_per_clone <- apply(escapees, 2, function(x){rownames(escapees)[x]})
escapees_per_clone <- lapply(escapees_per_clone, function(x){x[x != "Xist" & !is.na(x)]})
names(escapees_per_clone) <- c("C3March2022", "D11_March2022")

ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  ggplot(aes(x = Clone, y = value, fill = Condition)) + geom_boxplot(col = "grey") + 
  theme_paper() + ylab("ASE (d-score)") + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("black", "red")) + xlab("")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new//Fig5/escape_only_escapees.pdf")

ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(Clone == "C3 - March2022") %>%
  ggplot(aes(x = Condition, y = value, col = Condition, group = Condition)) + geom_point() + 
  geom_line(aes(group = name), col = "grey", alpha = 0.5) + ggpubr::stat_compare_means(paired = T) + theme_paper() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_manual(values = c("black", "red", "blue")) + xlab("") + ylab("ASE (d-score)") + 
  ggtitle("Clone3")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig5/clone3.pdf")

ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(Clone == "D11 - March2022") %>%
  ggplot(aes(x = Condition, y = value, col = Condition, group = Condition)) + geom_point() + 
  geom_line(aes(group = name), col = "grey", alpha = 0.5) + ggpubr::stat_compare_means(paired = T) + theme_paper() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_manual(values = c("black", "red", "blue")) + xlab("") + ylab("ASE (d-score)") + 
  ggtitle("Clone11")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig5/clone11.pdf")

# Clone C3
ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(Clone == "C3 - March2022") %>%
  pivot_wider(values_from = value, names_from = Condition) %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$name, escapee_annotation$symbol), ]$our_status2)) %>%
  ggplot(aes(x = `non-targeting control`, y = sgXist, col = escape_status)) + geom_point() + geom_abline(linetype = "dashed") + 
    scale_x_sqrt() + scale_y_sqrt() + theme_paper() + ggrepel::geom_text_repel(aes(label = name))
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig5/clone3_maplot.pdf")

# Clone D11
ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(Clone == "D11 - March2022") %>%
  pivot_wider(values_from = value, names_from = Condition) %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$name, escapee_annotation$symbol), ]$our_status2)) %>%
  ggplot(aes(x = `non-targeting control`, y = sgXist, col = escape_status)) + geom_point() + geom_abline(linetype = "dashed") + 
    scale_x_sqrt() + scale_y_sqrt() + theme_paper() + ggrepel::geom_text_repel(aes(label = name))
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig5/clone11_maplot.pdf")

### subset on other genes here? not sure
ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  pivot_wider(values_from = value, names_from = Condition) %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$name, escapee_annotation$symbol), ]$our_status2)) %>%
  add_column(difference = .$sgXist - .$`non-targeting control`) %>%
  dplyr::filter(Clone = "D11 - March2022") %>%
  dplyr::filter(`non-targeting control` < 0.6) %>%
  dplyr::filter(!is.na(escape_status)) %>%
  ggplot(aes(x = cut(`non-targeting control`, breaks = 10), y = difference)) + 
    geom_boxplot() + coord_flip()

ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  pivot_wider(values_from = value, names_from = Condition) %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$name, escapee_annotation$symbol), ]$our_status2)) %>%
  add_column(difference = .$sgXist - .$`non-targeting control`) %>%
  dplyr::filter(Clone == "D11 - March2022") %>%
  dplyr::filter(`non-targeting control` < 0.6) %>%
  dplyr::filter(!is.na(escape_status)) %>%
  ggplot(aes(x = cut(`non-targeting control`, breaks = 10), y = difference, fill = escape_status)) + 
    geom_boxplot() + coord_flip()

###

library(VGAM)

# define escapees here
coefs_dox_experiment <- readRDS("./ProcessedData/LM_coefficients_E6.rds")

ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

data_test_inactive <- counts_inactive(data_here[rownames(coefs_dox_experiment), ])
data_test_active <- counts_active(data_here[rownames(coefs_dox_experiment), ])

fitting_metadata_here <- colData(data_here)

all_coefs <- data.frame(do.call("rbind", lapply(1:nrow(data_test_active), function(i){
  tryCatch({
    y = as.numeric(data_test_inactive[i, ]) + 1
    N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ]) + 2
    fit1 <- glm(cbind(y, N-y) ~ Clone + Guide, family = "binomial", data = fitting_metadata_here)
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
all_coefs <- all_coefs[!apply(all_coefs, 1, function(x){all(is.na(x))}), !apply(all_coefs, 2, function(x){all(is.na(x))})]
all_coefs[is.na(all_coefs)] <- 0
#all_coefs[abs(all_coefs) > 5] <- sign(all_coefs[abs(all_coefs) > 5]) * 5

pheatmap::pheatmap(all_coefs)

data.frame(
  gene = rownames(all_coefs), 
  xist_kd = all_coefs$GuideP3, 
  xist_oe = coefs_dox_experiment[rownames(all_coefs), ]$ndDox7
) %>% 
  dplyr::filter(xist_kd > -4) %>% #### beware of the filtering
  ggplot(aes(xist_kd, xist_oe)) + geom_point() + 
    geom_smooth(method = "lm", linetype = "dashed", col = "grey") + 
    ggrepel::geom_text_repel(aes(label = gene)) + 
    theme_paper() + xlab("Xist knockdown") + ylab("Xist overexpression") + 
    ggpubr::stat_cor(label.x = -0.1, label.y = -9)
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig5/xist_oe_xist_ko.pdf")

plot_gene <- function(gene, save = F, path = NULL, suffix = "_escape_plot.pdf"){
  i = gene
  y = as.numeric(data_test_inactive[i, ])
  N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
  fit1 <- glm(cbind(y, N-y) ~ Clone + Guide, family = "binomial", data = fitting_metadata_here)

  pp <- data.frame(
    y = y,
    N = N,
    Clone = paste0(data_here$Clone, " - ", data_here$Experiment),
    Cov = data_here$ConditionClean
  ) %>% ggplot(aes(Cov, y / N, col = Clone)) + geom_point(size = 4) + coord_flip() + ylim(0, 1) +
    ggtitle(i) +
    geom_hline(yintercept = rev_logit(coef(fit1)[[1]]), linetype = "dashed") + 
    geom_line(aes(group = Clone), linetype = "dashed") + 
    theme_paper() + xlab("") + ylab("ASE (d-score)")
  
  if (save){
    ggsave(paste0(path, "./", gene, suffix), pp)
  } 
  
  return(pp)
}

plot_gene("Mid1ip1", save = F, path = NULL, suffix = "_escape_plot.pdf")

#### 
genes_show <- c("Xist", "Mecp2", "Kdm6a", "Kdm5c")
genes_of_choice <- setNames(rep("", nrow(ratios)), rownames(ratios))
genes_of_choice[genes_show] <- genes_show

ratios[apply(ratios, 1, function(x){sum(is.na(x))}) < 1, ] %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  #dplyr::filter(Condition == "Control") %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  ggplot(aes(Condition, y = reorder(name, position), fill = value)) + geom_tile(width = 0.9) +
    facet_wrap(~Clone) + 
    theme_paper() + theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
    theme(axis.ticks.y = element_blank()) + xlab("") + ylab("Chromosome position") + 
    scale_fill_gradientn(colors = c("#808000", "pink", "#008080")) + 
    scale_y_discrete(labels = genes_of_choice, position = "right") + ylab("")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig5/spatial_maps_xist_kd.pdf")

difference_c1 <- ratios[apply(ratios, 1, function(x){sum(is.na(x))}) < 1, ][,1] - 
  ratios[apply(ratios, 1, function(x){sum(is.na(x))}) < 1, ][,2]

difference_c2 <- ratios[apply(ratios, 1, function(x){sum(is.na(x))}) < 1, ][,3] - 
  ratios[apply(ratios, 1, function(x){sum(is.na(x))}) < 1, ][,4]

data.frame(
  x = difference_c1, 
  y = difference_c2, 
  gene = rownames(ratios[apply(ratios, 1, function(x){sum(is.na(x))}) < 1, ])
) %>%
  add_column(position = start(rowRanges(data_here[as.character(.$gene), ]))) %>%
  pivot_longer(-c(gene, position)) %>%
  ggplot(aes(x = position, y = value, col = name)) + geom_point() + 
  geom_smooth(span = 0.1)
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures_new/Fig5/spatial_maps_xist_kd_fc.pdf")

all_coefs %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-c(Gene)) %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status)) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(x = name, y = value, col = escape_status)) + geom_boxplot() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.2) +
  coord_flip() + theme_paper() + xlab("Covariate") +
  theme(aspect.ratio = 2) + theme(legend.position = "top") + geom_hline(yintercept = 0, linetype = "dashed") + 
  ylab("Model coefficient") + xlab("") #+ scale_x_discrete(labels = rev(c("Intercept", "Washout", "7d dox", "3d dox")))
#ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures/Fig3_prelim/spen_model_coefs.pdf")






### ### ### ### ### ### ### ### ### ### ### ### 
### Long term experiments
### First do analysis like this, then decide where this stuff comes in
### ### ### ### ### ### ### ### ### ### ### ### 

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


data <- readRDS("./ProcessedData/merged_dataset.rds")
data <- data[,!data$Guide %in% c("g13", "P3") & data$Experiment != "May2022"]

data[,data$Experiment == "September2022"]$Condition

# clean up names
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
  "Aux_0_Dox_14_WO_0_WOAuxNO" = "Dox (14d)",
  
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

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Figure 1
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Xist overexpression silences escapees in post-silencing contexts
# First, we analyze which genes escape in the baseline in our 3 different clones
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

data_baseline_x <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_NO_WOAuxNO","Aux_0_Dox_0_WO_NO_WOAuxNA" )]

# Check that escape in new E6 sample is comparable to the other two
# Make the same plots for E6: (replicate 1)
sample_nr = 3

# Prettier plot, with escapee-annotation
data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>% 
  ggplot(aes(total, ase, col = escape_status)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") +  ggtitle("Clone E6 (replicate1)") + 
  xlab("Total expression (allelic reads only)") + ylab("ASE (d-score)") + coord_fixed(ratio = 3)

data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>%
  dplyr::filter(!is.na(escape_status)) %>%
  mutate(escape_status = factor(escape_status, levels = c("silenced / variable", "facultative", "constitutive"))) %>%
  ggplot(aes(x = escape_status, y = ase, fill = escape_status)) + geom_boxplot() + theme_paper() + 
  xlab("") + ylab("ASE (d-score)") + scale_fill_manual(values = rev(scales::hue_pal()(3))) + ggtitle("") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## Rep2
sample_nr = 4

# Prettier plot, with escapee-annotation
data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>% 
  ggplot(aes(total, ase, col = escape_status)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") +  ggtitle("Clone E6 (replicate2)") + 
  xlab("Total expression (allelic reads only)") + ylab("ASE (d-score)") + coord_fixed(ratio = 3)

data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>%
  dplyr::filter(!is.na(escape_status)) %>%
  mutate(escape_status = factor(escape_status, levels = c("silenced / variable", "facultative", "constitutive"))) %>%
  ggplot(aes(x = escape_status, y = ase, fill = escape_status)) + geom_boxplot() + theme_paper() + 
  xlab("") + ylab("ASE (d-score)") + scale_fill_manual(values = rev(scales::hue_pal()(3))) + ggtitle("") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# Make the same plots for E6: (replicate 3, the new one)
sample_nr = 5

# Prettier plot, with escapee-annotation
data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>% 
  ggplot(aes(total, ase, col = escape_status)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") +  ggtitle("Clone E6 (replicate3)") + 
  xlab("Total expression (allelic reads only)") + ylab("ASE (d-score)") + coord_fixed(ratio = 3)

data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>%
  dplyr::filter(!is.na(escape_status)) %>%
  mutate(escape_status = factor(escape_status, levels = c("silenced / variable", "facultative", "constitutive"))) %>%
  ggplot(aes(x = escape_status, y = ase, fill = escape_status)) + geom_boxplot() + theme_paper() + 
  xlab("") + ylab("ASE (d-score)") + scale_fill_manual(values = rev(scales::hue_pal()(3))) + ggtitle("") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

### Now compare all the clones
data_here <- data[,data$ConditionClean == "Control"]

# get d-scores across genes + samples
ratios <- counts_inactive(data_baseline_x) / (counts_inactive(data_baseline_x) + counts_active(data_baseline_x))

# exclude genes with 0 counts in individual samples
rownames( ratios[!apply(ratios, 1, function(x){!any(is.na(x))}), ] )
ratios <- ratios[apply(ratios, 1, function(x){!any(is.na(x))}), ]

colnames(ratios) <- c("CL30", "CL31", "E6 (rep1)", "E6 (rep2)", "E6 (rep3)")

genes_show <- c("Xist", "Mecp2", "Kdm6a", "Kdm5c")
genes_of_choice <- setNames(rep("", nrow(ratios)), rownames(ratios))
genes_of_choice[genes_show] <- genes_show

# Plot d-scores for genes after ordering the genes by chromosomale coordinate
ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  dplyr::filter(Condition == "Control") %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  ggplot(aes(Clone, y = reorder(name, position), fill = value)) + geom_tile(width = 0.9) +
    theme_paper() + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
    theme(axis.ticks.y = element_blank()) + xlab("") + ylab("Chromosome position") + 
    scale_fill_gradientn(colors = c("#808000", "pink", "#008080")) + 
    scale_y_discrete(labels = genes_of_choice, position = "right") + ylab("")

# Plot d-scores for genes that escape in any sample, excluding genes with > 0.7 d-score
ratios <- ratios[rownames(ratios) == "Xist" | (apply(ratios, 1, function(x){any(x > 0.1)}) & !apply(ratios, 1, function(x){any(x > 0.7)})), ]
ratios <- ratios[order(rowMeans(ratios), decreasing = T), ]

colnames(ratios) <- c("CL30", "CL31", "E6 (rep1)", "E6 (rep2)", "E6 (rep3)")
pheatmap::pheatmap(ratios, cluster_rows = F, cluster_cols = F, fontsize = 20, show_rownames =  F, 
                   color=colorRampPalette(c("#808000", "pink", "#008080"))(50))

# Now we define the baseline escapees per sample as the genes with d-score > 0.1
escapees <- apply(ratios, 2, function(x){x > 0.1})
escapees_per_clone <- apply(escapees, 2, function(x){rownames(escapees)[x]})
escapees_per_clone <- lapply(escapees_per_clone, function(x){x[x != "Xist"]})
names(escapees_per_clone) <- c("CL30 - December2021", "CL31 - December2021", "E6 - October2021", "E6 - March2022", "E6 - September2022")

# How many do we get per clone and do they overlap?
ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  dplyr::filter(Condition == "Control") %>%
  dplyr::filter(value > 0.1) %>%
  group_by(name) %>%
  mutate(escape_in_clones = n() > 3) %>%
  ggplot(aes(x = Clone, fill = escape_in_clones)) + geom_bar() + theme_paper() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) + xlab("") + ylab("Number of genes with d-score > 0.1") + 
  scale_fill_manual(values = c("grey", "darkred"), name = "Escape in all clones")

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Figure 2b
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Now, we look at the effect of Xist-treatment on escapee expression
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### Now we look at Xist overexpression and how that silences escapees
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_NO_WOAuxNO", "Aux_0_Dox_3_WO_NO_WOAuxNO", "Aux_0_Dox_7_WO_NO_WOAuxNO", 
                                         "Aux_0_Dox_14_WO_NO_WOAuxNO", "Aux_0_Dox_21_WO_NO_WOAuxNO")]

ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))
ratios <- ratios[apply(ratios, 1, function(x){!any(is.na(x))}), ]

# For Xist-expression, we need a scaling factor per sample to account for differences in seq depth
size_factors <- colSums(counts(data_here)) / colSums(counts(data_here))[[1]]

# check that Xist-overexpression works
data.frame(
  Dox = data_here$ndDox,
  Clone = paste0(data_here$Clone, " - ", data_here$Experiment), 
  Expression = as.numeric(counts(data_here["Xist", ]) / colSums(counts(data_here)) * 1e6)
) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", paste0("Dox (", .$Dox, " days)"))) %>%
  mutate(Condition = factor(Condition, levels = c("Control", paste0("Dox (", c(3, 7, 14, 21), " days)")))) %>%
  ggplot(aes(x = Condition, y = Expression)) + geom_boxplot() + geom_point() + coord_flip() + 
  geom_line(aes(group = Clone, linetype = Clone, col = Clone)) + theme_paper() + 
  xlab("") + ylab("Xist Expression (normalized CPM)")

# Now look at the global effect of Xist-overexpression on escape
ratios %>%
  t() %>% data.frame() %>%
  add_column("Dox" = data_here$ndDox) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Dox, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", paste0("Dox (", .$Dox, " days)"))) %>%
  mutate(Condition = factor(Condition, levels = c("Control", paste0("Dox (", c(3, 7, 14, 21), " days)")))) %>%
  ggplot(aes(x = Clone, y = value, fill = Condition)) + geom_boxplot(col = "grey") + 
  theme_paper() + ylab("ASE (d-score)") + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("black", "red", "blue", "darkred", "darkblue"), labels = c("Control", "Dox (3d)", "Dox (7d)", "Dox (14d)", "Dox (21d)"))

# Look at gene-wise silencing, after n days
ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  dplyr::filter(Condition %in% c("Control", "Dox (7d)", "Dox (7d) - washout") & Clone == "E6 - September2022") %>%
  ggplot(aes(x = Condition, y = value, col = Condition, group = Condition)) + geom_point() + 
  geom_line(aes(group = name), col = "grey", alpha = 0.5) + theme_paper() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_manual(values = c("black", "red", "blue")) + xlab("") + ylab("ASE (d-score)") + 
  ggtitle("Silencing after 7d (E6)")

ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  dplyr::filter(Condition %in% c("Control", "Dox (14d)", "Dox (14d) - washout") & Clone == "E6 - September2022") %>%
  ggplot(aes(x = Condition, y = value, col = Condition, group = Condition)) + geom_point() + 
  geom_line(aes(group = name), col = "grey", alpha = 0.5) + theme_paper() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_manual(values = c("black", "red", "blue")) + xlab("") + ylab("ASE (d-score)") + 
  ggtitle("Silencing after 14d (E6)")

ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  dplyr::filter(Condition %in% c("Control", "Dox (21d)", "Dox (21d) - washout") & Clone == "E6 - September2022") %>%
  ggplot(aes(x = Condition, y = value, col = Condition, group = Condition)) + geom_point() + 
  geom_line(aes(group = name), col = "grey", alpha = 0.5) + theme_paper() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_manual(values = c("black", "red", "blue")) + xlab("") + ylab("ASE (d-score)") + 
  ggtitle("Silencing after 21d (E6)")

# Highlight genes that still escape: 
sample_nr = data$ConditionClean == "Dox (21d)"

data.frame(
  total = rowSums(counts_inactive(data[,sample_nr])) + rowSums(counts_active(data[,sample_nr])),
  ase = rowSums(counts_inactive(data[,sample_nr])) / (rowSums(counts_inactive(data[,sample_nr])) +
                                                                   rowSums(counts_active(data[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>% 
  ggplot(aes(total, ase, col = escape_status)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") +  ggtitle("E6, rep3, 21d Xist") + 
  xlab("Total expression (allelic reads only)") + ylab("ASE (d-score)") + coord_fixed(ratio = 3)

# Now, check how silencing occurs across the entire chromosome
### Now compare all the clones
data_here <- data[,data$ConditionClean %in% c("Control", paste0("Dox (", c(3, 7, 14, 21), "d)"))]

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
  dplyr::filter(Condition == "Control") %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  mutate(Condition = factor(Condition, levels = c("Control", paste0("Dox (", c(3, 7, 14, 21), "d)")))) %>%
  ggplot(aes(Condition, y = reorder(name, position), fill = value)) + geom_tile(width = 0.9) +
    theme_paper() + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
    theme(axis.ticks.y = element_blank()) + xlab("") + ylab("Chromosome position") + 
    scale_fill_gradientn(colors = c("#808000", "pink", "#008080")) + 
    scale_y_discrete(labels = genes_of_choice, position = "right") + ylab("") + 
    facet_wrap(~Clone, nrow = 1)

# ask specifically, what is the difference in escape between 14d and 21d

data_here <- data[,data$Experiment == "September2022"]
data_here <- data_here[,data_here$ConditionClean %in% c("Dox (14d)", "Dox (21d)")]

plot_scatter <- function(data_subset){
  ratios <- counts_inactive(data_subset) / (counts_inactive(data_subset) + counts_active(data_subset))
  # exclude genes with 0 counts in individual samples
  #rownames( ratios[!apply(ratios, 1, function(x){!any(is.na(x))}), ] )
  #ratios <- ratios[apply(ratios, 1, function(x){!any(is.na(x))}), ]
  
  df_plot <- data.frame(Gene = rownames(ratios), sample1 = ratios[,1], sample2 = ratios[,2] ) %>% 
    add_column(highlight = .$sample1 > 0.1 | .$sample2 > 0.1) %>%
    add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2))
  
  df_plot %>% 
    ggplot(aes(x = sample1, y = sample2, col = escape_status)) + geom_point() + theme_paper() + geom_abline(linetype = "dashed") +
      ggrepel::geom_text_repel(aes(label = Gene), data = df_plot %>% dplyr::filter(highlight)) + 
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey") + 
      geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey") + labs("")
}

# Make clear pairwise comparisons

# 0d vs 21d 
data_here <- data[,data$Experiment == "September2022"]
data_here <- data_here[,data_here$ConditionClean %in% c("Control", "Dox (21d)")]
plot_scatter(data_here) + xlab("E6 rep3, Control (d-score)") + ylab("E6 rep3, 21d dox (d-score)")
ggsave("./Figures_new/Fig3/LongXist_Control_vs_21d.pdf")

# 7d vs 21d 
data_here <- data[,data$Experiment == "September2022"]
data_here <- data_here[,data_here$ConditionClean %in% c("Dox (7d)", "Dox (21d)")]
plot_scatter(data_here) + xlab("E6 rep3, 7d dox (d-score)") + ylab("E6 rep3, 21d dox (d-score)")
ggsave("./Figures_new/Fig3//LongXist_7d_vs_21d.pdf")

# 14d vs 21d 
data_here <- data[,data$Experiment == "September2022"]
data_here <- data_here[,data_here$ConditionClean %in% c("Dox (14d)", "Dox (21d)")]
plot_scatter(data_here) + xlab("E6 rep3, 14d dox (d-score)") + ylab("E6 rep3, 21d dox (d-score)")
ggsave("./Figures_new/Fig3//LongXist_14d_vs_21d.pdf")

# washout comparison 14d
data_here <- data[,data$Experiment == "September2022"]
data_here <- data_here[,data_here$ConditionClean %in% c("Dox (14d) - washout", "Dox (14d) - washout_long")]
plot_scatter(data_here[,c(2,1)]) + xlab("E6 rep3, 14d dox + washout (d-score)") + ylab("E6 rep3, 14d dox + washout long (d-score)")
ggsave("./Figures_new/Fig3//LongXist_14d_washout_vs_washout_long.pdf")

# washout comparison 21d
data_here <- data[,data$Experiment == "September2022"]
data_here <- data_here[,data_here$ConditionClean %in% c("Dox (21d) - washout", "Dox (21d) - washout_long")]
plot_scatter(data_here) + xlab("E6 rep3, 21d dox + washout (d-score)") + ylab("E6 rep3, 21d dox + washout long (d-score)")
ggsave("./Figures_new/Fig3//LongXist_21d_washout_vs_washout_long.pdf")


# export dataframe in excel to be able to check ratios offline
data_here <- data[,data$Experiment == "September2022"]
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))
colnames(ratios) <- data_here$ConditionClean
ratios$escape_status <- unlist(escapee_annotation[match(rownames(ratios), escapee_annotation$symbol), ]$our_status2)

ratios$dASE_Cont_7d <- ratios$`Dox (7d)` - ratios$Control
ratios$dASE_Cont_14d <- ratios$`Dox (14d)` - ratios$Control
ratios$dASE_Cont_21d <- ratios$`Dox (21d)` - ratios$Control
ratios$dASE_7d_14d <- ratios$`Dox (14d)` - ratios$`Dox (7d)`
ratios$dASE_7d_21d <- ratios$`Dox (21d)` - ratios$`Dox (7d)`
ratios$dASE_14d_21d <- ratios$`Dox (21d)` - ratios$`Dox (14d)`
ratios$dASE_14d_washout_washoutlong <- ratios$`Dox (14d) - washout` - ratios$`Dox (14d) - washout_long`
ratios$dASE_21d_washout_washoutlong <- ratios$`Dox (21d) - washout` - ratios$`Dox (21d) - washout_long`

write.csv(ratios, "./ProcessedData/LongXist_experiment.csv")

### This is for Figure 3, the reversibility aspect
### Now we look at Xist overexpression and how that silences escapees
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_NO_WOAuxNO", "Aux_0_Dox_3_WO_NO_WOAuxNO", "Aux_0_Dox_3_WO_YES_WOAuxNO", 
                                         "Aux_0_Dox_7_WO_NO_WOAuxNO", "Aux_0_Dox_7_WO_YES_WOAuxNO", "Aux_0_Dox_14_WO_NO_WOAuxNO", 
                                         "Aux_0_Dox_21_WO_YES_WOAuxNO", "Aux_0_Dox_7_WO_NO_WOAuxNO", "Aux_0_Dox_7_WO_YES_WOAuxNO", 
                                         "Aux_0_Dox_14_WO_YES_WOAuxNO", "Aux_0_Dox_21_WO_YES_long_WOAuxNO", "Aux_0_Dox_0_WO_NO_WOAuxNO", 
                                         "Aux_0_Dox_21_WO_NO_WOAuxNO", "Aux_0_Dox_14_WO_YES_WOAuxNO", "Aux_0_Dox_14_WO_YES_long_WOAuxNO")]

ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))
ratios <- ratios[apply(ratios, 1, function(x){!any(is.na(x))}), ]

# For Xist-expression, we need a scaling factor per sample to account for differences in seq depth
size_factors <- colSums(counts(data_here)) / colSums(counts(data_here))[[1]]

# check that Xist-overexpression works
data.frame(
  Dox = data_here$ndDox,
  Washout = data_here$washout, 
  Clone = paste0(data_here$Clone, " - ", data_here$Experiment), 
  Expression = as.numeric(counts(data_here["Xist", ]) / colSums(counts(data_here)) * 1e6)
) %>%
  add_column(TimeSeries = ifelse(.$Dox == 0, list(c("3", "7", "14", "21")), .$Dox)) %>%
  unnest(TimeSeries) %>%
  mutate(TimeSeries = paste0(unlist(TimeSeries), " days dox treatment")) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", .$Washout)) %>%
  mutate(TimeSeries = factor(TimeSeries, levels = paste0(c(3, 7, 14, 21), " days dox treatment"))) %>%
  ggplot(aes(x = Condition, y = Expression)) + geom_boxplot() + geom_point() + coord_flip() + 
    geom_line(aes(group = Clone, linetype = Clone, col = Clone)) + theme_paper() + 
    xlab("") + ylab("Xist Expression (normalized CPM)") + facet_wrap(~TimeSeries, nrow = 2) + 
    scale_x_discrete(labels = c("Control", "Dox", "Washout", "Washout - long"))

# Now look at the global effect of Xist-overexpression on escape
ratios %>%
  t() %>% data.frame() %>%
  add_column("Dox" = data_here$ndDox) %>%
  add_column("Washout" = data_here$washout) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Dox, Washout, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  add_column(TimeSeries = ifelse(.$Dox == 0, list(c("3", "7", "14", "21")), .$Dox)) %>%
  unnest(TimeSeries) %>%
  mutate(TimeSeries = paste0(unlist(TimeSeries), " days dox treatment")) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", .$Washout)) %>%
  mutate(TimeSeries = factor(TimeSeries, levels = paste0(c(3, 7, 14, 21), " days dox treatment"))) %>%
  ggplot(aes(x = Clone, y = value, fill = Condition)) + geom_boxplot(col = "grey") + 
    theme_paper() + ylab("ASE (d-score)") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
    #scale_fill_manual(values = c("black", "red", "blue", "darkred", "darkblue"), labels = c("Control", "Dox (3d)", "Dox (7d)")) + 
    facet_wrap(~TimeSeries)


sample_nr = data$ConditionClean == "Control"
data.frame(
  total = rowSums(counts_inactive(data[,sample_nr])) + rowSums(counts_active(data[,sample_nr])),
  ase = rowSums(counts_inactive(data[,sample_nr])) / (rowSums(counts_inactive(data[,sample_nr])) +
                                                        rowSums(counts_active(data[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>% 
  ggplot(aes(total, ase, col = escape_status)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") +  ggtitle("Clone E6 (replicate3)") + 
  xlab("Total expression (allelic reads only)") + ylab("ASE (d-score)") + coord_fixed(ratio = 3)

sample_nr = data$ConditionClean == "Dox (21d)"
data.frame(
  total = rowSums(counts_inactive(data[,sample_nr])) + rowSums(counts_active(data[,sample_nr])),
  ase = rowSums(counts_inactive(data[,sample_nr])) / (rowSums(counts_inactive(data[,sample_nr])) +
                                                        rowSums(counts_active(data[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>% 
  ggplot(aes(total, ase, col = escape_status)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") +  ggtitle("Clone E6 - silencing 21d") + 
  xlab("Total expression (allelic reads only)") + ylab("ASE (d-score)") + coord_fixed(ratio = 3)

sample_nr = data$ConditionClean == "Dox (21d) - washout_long"
data.frame(
  total = rowSums(counts_inactive(data[,sample_nr])) + rowSums(counts_active(data[,sample_nr])),
  ase = rowSums(counts_inactive(data[,sample_nr])) / (rowSums(counts_inactive(data[,sample_nr])) +
                                                        rowSums(counts_active(data[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ]$our_status2)) %>% 
  ggplot(aes(total, ase, col = escape_status)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") +  ggtitle("Clone E6, silencing + 21d washout") + 
  xlab("Total expression (allelic reads only)") + ylab("ASE (d-score)") + coord_fixed(ratio = 3)

## Now analyze reveersion after washout
data_here <- data[,data$ConditionClean %in% c("Control", paste0("Dox (", c(3, 7, 14, 21), "d)"), 
                                              paste0("Dox (", c(3, 7, 14, 21), "d) - washout"), 
                                              paste0("Dox (", c(3, 7, 14, 21), "d) - washout_long"))]

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
factor_order <-  c("Control", "Dox (7d)", "Dox (7d) - washout", "Dox (14d)", "Dox (14d) - washout", 
                   "Dox (21d)", "Dox (21d) - washout", "Dox (21d) - washout_long")

ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  dplyr::filter(Clone == "E6 - September2022") %>%
  #mutate(Condition = factor(Condition, levels =)) %>%
  mutate(Condition = factor(Condition, levels = factor_order)) %>%
  ggplot(aes(Condition, y = reorder(name, position), fill = value)) + geom_tile(width = 0.9) +
  theme_paper() + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
  theme(axis.ticks.y = element_blank()) + xlab("") + ylab("Chromosome position") + 
  scale_fill_gradientn(colors = c("#808000", "pink", "#008080")) + 
  scale_y_discrete(labels = genes_of_choice, position = "right") + ylab("") + 
  facet_wrap(~Clone, nrow = 1)

### bottom line: 
# experiment looks good as a replicate (3d is missing unfortunately)
# silencing mostly saturates at 14d
# washout doesnt come back anymore
# need to quantify relative reversibility of different regions



