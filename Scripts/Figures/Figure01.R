#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### Plots for figure  1
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 

setwd("~/Desktop/Projects/XChromosome_Antonia/")

library(tidyverse)
library(ggplot2)
library(scran)
library(scater)


theme_paper <- function(){
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        text = element_text(size = 20))
}


source("./Scripts/auxiliary.R")

theme_paper <- function(){
  theme(panel.background = element_blank(),
        panel.grid = element_line(colour= "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5)),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        text = element_text(size = 20))
}

# Now read our pre-processed dataset
data <- readRDS("./ProcessedData/merged_dataset.rds")
#data <- computeSumFactors(data)

# Which samples do we have in the dataset?
table(data$Clone)
table(data$ndTreatment)
table(data$ndAux)
table(data$ndDox)
table(data$washout)
table(data$WOwithAux)
table(data$Experiment)

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

# save data before filtering
data_before_filtering <-  data[seqnames(data) == "X", ]

# remove genes with < 10 allelic reads per sample 
data <- data[(rowSums(counts_inactive(data)) + rowSums(counts_active(data))) / ncol(data) > 10, ]

# 
data_with_autosomes <- data
data <- data[seqnames(data) == "X", ]

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Xist overexpression silences escapees in post-silencing contexts
# First, we analyze which genes escape in the baseline in our 3 different clones
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

data_baseline_x <- data[,data$Condition == "Aux_0_Dox_0_WO_0_WOAuxNO"]

# In the basal state, show allelic coverage: 
# which genes are expressed and for which do we have allele-specific signal?
# distribution of allelic ratios on the x and in autosomes

# Show for one C30 example: 

sample_show = 1
data.frame(
  counts_total = counts(data_before_filtering[,sample_show])[[1]], 
  counts_allelic = (counts_active(data_before_filtering[,sample_show]) +  counts_inactive(data_before_filtering[,sample_show]))[[1]]
) %>%
  ggplot(aes(x = counts_total + 1, y = counts_allelic + 1)) + geom_point() + scale_x_log10() + scale_y_log10() + 
  theme_paper() + geom_vline(xintercept = 10, linetype = "dashed") + xlab("Total reads mapped") + ylab("Allele-specific reads mapped") + 
  geom_abline(linetype = "dashed")
ggsave("./FiguresIllustrator/FigS1/mapped_reads_c30.pdf")

data.frame(
  chromosome = as.character(seqnames(rowRanges(data_with_autosomes))) == "X", 
  total = (counts_active(data_with_autosomes[,sample_show]) + counts_inactive(data_with_autosomes[,sample_show]))[[1]], 
  ase =  (counts_active(data_with_autosomes[,sample_show]) / (counts_active(data_with_autosomes[,sample_show]) + counts_inactive(data_with_autosomes[,sample_show])))[[1]]
) %>%
  ggplot(aes(x = ase, fill = chromosome)) + geom_density() + theme_paper() + ylab("Density") + xlab("B6 / (B6 + CAST)")
ggsave("./FiguresIllustrator/FigS1/density_ase_c30.pdf")

# Show for one E6 example: 

sample_show = 31
data.frame(
  counts_total = counts(data_before_filtering[,sample_show])[[1]], 
  counts_allelic = (counts_active(data_before_filtering[,sample_show]) +  counts_inactive(data_before_filtering[,sample_show]))[[1]]
) %>%
  ggplot(aes(x = counts_total + 1, y = counts_allelic + 1)) + geom_point() + scale_x_log10() + scale_y_log10() + 
  theme_paper() + geom_vline(xintercept = 10, linetype = "dashed") + xlab("Total reads mapped") + ylab("Allele-specific reads mapped") + 
  geom_abline(linetype = "dashed")
ggsave("./FiguresIllustrator/FigS1/mapped_reads_e6.pdf")

data.frame(
  chromosome = as.character(seqnames(rowRanges(data_with_autosomes))) == "X", 
  total = (counts_active(data_with_autosomes[,sample_show]) + counts_inactive(data_with_autosomes[,sample_show]))[[1]], 
  ase =  (counts_active(data_with_autosomes[,sample_show]) / (counts_active(data_with_autosomes[,sample_show]) + counts_inactive(data_with_autosomes[,sample_show])))[[1]]
) %>%
  ggplot(aes(x = ase, fill = chromosome)) + geom_density() + theme_paper() + ylab("Density") + xlab("B6 / (B6 + CAST)")
ggsave("./FiguresIllustrator/FigS1/density_ase_e6.pdf")

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
  xlab("Total expression (allelic reads)") + ylab("ASE (d-score)")
ggsave("./FiguresIllustrator/FigS1/cl30_escape.pdf", width = 12, height = 8)

sample_nr = 3 # E6

data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle("Clone 30") +
  xlab("Total expression (allelic reads)") + ylab("ASE (d-score)")
ggsave("./FiguresIllustrator/FigS1/e6_escape.pdf", width = 12, height = 8)


### Now compare all the clones
data_here <- data[,data$ConditionClean == "Control"]

# get d-scores across genes + samples
ratios <- counts_inactive(data_baseline_x) / (counts_inactive(data_baseline_x) + counts_active(data_baseline_x))

# Genes with average ASE > 0.7 are likely mapping artefacts: 
genes_out <- names(rowMeans(ratios, na.rm = T)[rowMeans(ratios, na.rm = T) > 0.8 & rownames(ratios) != "Xist"])
genes_out_na <- rownames( ratios[!apply(ratios, 1, function(x){!any(is.na(x))}), ] )

genes_out <- c(genes_out, genes_out_na)

# exclude genes with 0 counts in individual samples
ratios <- ratios[!rownames(ratios) %in% genes_out, ]

colnames(ratios) <- c("CL30", "CL31", "E6 (rep1)", "E6 (rep2)", "E6 (rep3)", "CL31.16", "CL30.7")

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
  scale_fill_gradient2(midpoint = 0.5, mid = "#117733", low = "#BEBEBE", high = "black", limits = c(0, 1), na.value = "white") + 
  scale_y_discrete(labels = genes_of_choice, position = "right") + ylab("")
ggsave("./FiguresIllustrator/FigS1/heatmap_ordered.pdf", width = 10, height = 10)

# Plot d-scores for genes that escape in any sample, excluding genes with > 0.7 d-score
ratios <- ratios[rownames(ratios) == "Xist" | (apply(ratios, 1, function(x){any(x > 0.1)}) & !apply(ratios, 1, function(x){any(x > 0.7)})), ]
ratios <- ratios[order(rowMeans(ratios), decreasing = T), ]

# Now we define the baseline escapees per sample as the genes with d-score > 0.1
escapees <- apply(ratios, 2, function(x){x > 0.1})
escapees_per_clone <- apply(escapees, 2, function(x){rownames(escapees)[x]})
escapees_per_clone <- lapply(escapees_per_clone, function(x){x[x != "Xist"]})
names(escapees_per_clone) <- c("CL30 - December2021", "CL31 - December2021", "E6 - October2021", "E6 - March2022", "E6 - September2022", 
                               "CL30 - October2022", "CL31 - October2022")

# How many do we get per clone and do they overlap?
ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  dplyr::filter(Condition == "Control") %>%
  dplyr::filter(value > 0.1) %>%
  group_by(name) %>%
  mutate(escape_in_clones = n() >= 7) %>%
  ggplot(aes(x = Clone, fill = escape_in_clones)) + geom_bar() + theme_paper() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) + xlab("") + ylab("Number of genes with d-score > 0.1") + 
  scale_fill_manual(values = c("grey", "darkred"), name = "Escape in all clones")
ggsave("./FiguresIllustrator/FigS1/number_of_escapees_per_clone.pdf", width = 8, height = 8)

ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  dplyr::filter(Condition == "Control") %>%
  dplyr::filter(value > 0.1) %>%
  group_by(name) %>%
  mutate(escape_in_clones = n() >= 7) %>%
  add_column(escape_group = rowData(data[.$name, ])$EscapeAnnotation) %>%
  ggplot(aes(x = Clone, fill = escape_in_clones)) + geom_bar() + theme_paper() + 
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) + xlab("") + ylab("Number of genes with d-score > 0.1") + 
    scale_fill_manual(values = c("grey", "darkred"), name = "Escape in all clones") + facet_wrap(~escape_group)
ggsave("./FiguresIllustrator/FigS1/number_of_escapees_per_clone_stratified.pdf", width = 8, height = 8)

# We are using a meta-analysis of escapees across studies to classify our escapees.
# We plot this in generate_sce.R

## For the main figure: 
# a) experimental outline
# b) Xist RNA fish + dox
# c) heatmap of silencing across x chromosome

### Now we look at Xist overexpression and how that silences escapees
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO", 
                                         "Aux_0_Dox_14_WO_0_WOAuxNO", "Aux_0_Dox_21_WO_0_WOAuxNO")]
data_here <- data_here[,data_here$Clone == "E6"]

ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))
ratios <- ratios[!rownames(ratios) %in% genes_out, ]

# For Xist-expression, we need a scaling factor per sample to account for differences in seq depth
size_factors <- colSums(counts(data_here)) / colSums(counts(data_here))[[1]]

# check that Xist-overexpression works
data.frame(
  Dox = data_here$ndDox,
  Clone = paste0(data_here$Clone, " - ", data_here$Experiment), 
  Expression = as.numeric(counts(data_here["Xist", ]) / colSums(counts(data_here)) * 1e6)
) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", paste0("Dox (", .$Dox, " days)"))) %>%
  mutate(Condition = factor(Condition, levels = rev(c("Control", paste0("Dox (", c(3, 7, 14, 21), " days)"))))) %>%
  ggplot(aes(x = Condition, y = Expression)) + 
  stat_summary(stat = "mean", geom = "bar", fill = "grey", col = "black") + 
    coord_flip() + theme_paper() + 
    xlab("") + ylab("Expression (CPM)") + geom_jitter(size = 3, width = 0.1) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 320000)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("./FiguresIllustrator/Fig1/xist_overexpression.pdf", width = 4, height = 6)

# ratios %>%
#   t() %>% data.frame() %>%
#   add_column("Dox" = data_here$ndDox) %>%
#   add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
#   pivot_longer(-c(Dox, Clone)) %>%
#   group_by(Clone) %>%
#   dplyr::filter(name %in% intersect(intersect(escapees_per_clone[[3]], escapees_per_clone[[4]]), escapees_per_clone[[5]])) %>%
#   add_column(Condition = ifelse(.$Dox == 0, "Control", paste0("Dox (", .$Dox, " days)"))) %>%
#   mutate(Condition = factor(Condition, levels = rev(c("Control", paste0("Dox (", c(3, 7, 14, 21), " days)"))))) %>%
#   ggplot(aes(x = reorder(name, -value), y = Condition, fill = value)) + geom_tile() + 
#     scale_fill_gradient2(midpoint = 0.5, mid = "#117733", low = "#BEBEBE", high = "black", limits = c(0, 1), na.value = "white") + 
#     theme_paper() + theme(axis.line.x =element_blank(), axis.line.y =element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) + 
#     theme(legend.position = "left") + ylab("") + theme(panel.border=element_blank()) 
# ggsave("./FiguresIllustrator/Fig1/heatmap_escape.pdf", width = 8, height = 5.5)

# make the same heatmap with all genes, not just escapees: 

genes_show <- c("Xist", "Mecp2", "Kdm6a", "Kdm5c")
genes_of_choice <- setNames(rep("", nrow(ratios)), rownames(ratios))
genes_of_choice[genes_show] <- genes_show

ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  mutate(Condition = factor(Condition, levels = rev(c("Control", "Dox (3d)", "Dox (7d)", "Dox (14d)", "Dox (21d)")))) %>%
  ggplot(aes(reorder(name, position), y = Condition, fill = value, col = value)) + geom_tile(width = 0.9) +
    scale_fill_gradient2(midpoint = 0.5, mid = "#117733", low = "#BEBEBE", high = "black", limits = c(0, 1), na.value = "white") + 
    scale_colour_gradient2(midpoint = 0.5, mid = "#117733", low = "#BEBEBE", high = "black", limits = c(0, 1), na.value = "white") + 
    scale_y_discrete(labels = genes_of_choice, position = "right") + 
    theme_paper() + 
    theme(axis.ticks.x = element_blank(), panel.grid.major = element_blank()) + 
    xlab("") + 
    ylab("") + theme(panel.border = element_blank(), axis.line=element_blank()) + 
    scale_x_discrete(labels = genes_of_choice, position = "bottom")
ggsave("./FiguresIllustrator/Fig1/heatmap_ordered.pdf", width = 10, height = 10)

ratios %>%
  t() %>% data.frame() %>%
  add_column("Dox" = data_here$ndDox) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Dox, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", paste0("Dox (", .$Dox, " days)"))) %>%
  mutate(Condition = factor(Condition, levels = (c("Control", paste0("Dox (", c(3, 7, 14, 21), " days)"))))) %>%
  add_column(escape_status = unlist(rowData(data_here)[match(.$name, rownames(data_here)), ]$EscapeAnnotation)) %>%
  mutate(escape_status = factor(escape_status, c("silenced / variable", "facultative", "constitutive"))) %>%
  dplyr::filter(!is.na(escape_status)) %>%
  ggplot(aes(x = Condition, y = value, fill = escape_status)) + 
    #geom_jitter(position = position_jitterdodge(jitter.width = 0.2)) + 
    geom_boxplot(col = "black", outlier.color = "grey") + 
    theme_paper() + ylab("ASE (d-score)") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
    scale_fill_manual(values = c("grey", "darkgreen", "orange"), labels = c("Silenced / Variable", "Facultative", "Constitutive")) + 
    labs(fill = "Escape Category") + xlab("")
ggsave("./FiguresIllustrator/Fig1/escapee_silencing_boxplots.pdf", width = 9, height = 6)

# For the supplement, also plot stratified by replicates: 

ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  mutate(Condition = factor(Condition, levels = rev(c("Control", "Dox (3d)", "Dox (7d)", "Dox (14d)", "Dox (21d)")))) %>%
  ggplot(aes(reorder(name, position), y = interaction(Clone, Condition, sep = " -- "), group = Clone, fill = value, col = value)) + geom_tile(width = 0.9) +
    scale_fill_gradient2(midpoint = 0.5, mid = "#117733", low = "#BEBEBE", high = "black", limits = c(0, 1), na.value = "white") + 
    scale_colour_gradient2(midpoint = 0.5, mid = "#117733", low = "#BEBEBE", high = "black", limits = c(0, 1), na.value = "white") + 
    scale_y_discrete(labels = genes_of_choice, position = "right") + 
    theme_paper() + 
    theme(axis.ticks.x = element_blank(), panel.grid.major = element_blank()) + 
    xlab("") + 
    ylab("") + theme(panel.border = element_blank(), axis.line=element_blank()) + 
    scale_x_discrete(labels = genes_of_choice, position = "bottom")
ggsave("./FiguresIllustrator/FigS1/heatmap_ordered_by_clones.pdf", width = 10, height = 10)

ratios %>%
  t() %>% data.frame() %>%
  add_column("Dox" = data_here$ndDox) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Dox, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", paste0("Dox (", .$Dox, " days)"))) %>%
  mutate(Condition = factor(Condition, levels = (c("Control", paste0("Dox (", c(3, 7, 14, 21), " days)"))))) %>%
  add_column(escape_status = unlist(rowData(data_here)[match(.$name, rownames(data_here)), ]$EscapeAnnotation)) %>%
  mutate(escape_status = factor(escape_status, c("silenced / variable", "facultative", "constitutive"))) %>%
  dplyr::filter(!is.na(escape_status)) %>%
  ggplot(aes(x = Condition, y = value, fill = Clone)) + 
    #geom_jitter(position = position_jitterdodge(jitter.width = 0.2)) + 
    geom_boxplot(col = "black", outlier.color = "grey") + 
    theme_paper() + ylab("ASE (d-score)") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
    #scale_fill_manual(values = c("grey", "darkgreen", "orange"), labels = c("Silenced / Variable", "Facultative", "Constitutive")) + 
    labs(fill = "Escape Category") + xlab("")
ggsave("./FiguresIllustrator/FigS1/escapee_silencing_boxplots_per_clone.pdf", width = 9, height = 6)

### For supplement, plot the heatmap etc in CL30/CL31
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO", 
                                         "Aux_0_Dox_14_WO_0_WOAuxNO", "Aux_0_Dox_21_WO_0_WOAuxNO")]
data_here <- data_here[,data_here$Clone != "E6"]

ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))
ratios <- ratios[!rownames(ratios) %in% genes_out, ]

# For Xist-expression, we need a scaling factor per sample to account for differences in seq depth
size_factors <- colSums(counts(data_here)) / colSums(counts(data_here))[[1]]

# check that Xist-overexpression works
data.frame(
  Dox = data_here$ndDox,
  Clone = paste0(data_here$Clone, " - ", data_here$Experiment), 
  Expression = as.numeric(counts(data_here["Xist", ]) / colSums(counts(data_here)) * 1e6)
) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", paste0("Dox (", .$Dox, " days)"))) %>%
  mutate(Condition = factor(Condition, levels = rev(c("Control", paste0("Dox (", c(3, 7, 14, 21), " days)"))))) %>%
  ggplot(aes(x = Condition, y = Expression)) + 
  stat_summary(stat = "mean", geom = "bar", fill = "grey", col = "black") + 
  coord_flip() + theme_paper() + 
  xlab("") + ylab("Xist Expression (normalized CPM)") + geom_jitter(size = 3, width = 0.1) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 260000))
ggsave("./FiguresIllustrator/FigS1/xist_overexpression_cl30_31.pdf", width = 6, height = 8)

genes_show <- c("Xist", "Mecp2", "Kdm6a", "Kdm5c")
genes_of_choice <- setNames(rep("", nrow(ratios)), rownames(ratios))
genes_of_choice[genes_show] <- genes_show

ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  #dplyr::filter(value < 0.7) %>%
  #dplyr::filter(Clone == "E6") %>%
  mutate(Condition = factor(Condition, levels = rev(c("Control", "Dox (3d)", "Dox (7d)", "Dox (14d)", "Dox (21d)")))) %>%
  ggplot(aes(reorder(name, position), y = Condition, fill = value, col = value)) + geom_tile(width = 0.9) +
  scale_fill_gradient2(midpoint = 0.5, mid = "#117733", low = "#BEBEBE", high = "black", limits = c(0, 1), na.value = "white") + 
  scale_colour_gradient2(midpoint = 0.5, mid = "#117733", low = "#BEBEBE", high = "black", limits = c(0, 1), na.value = "white") + 
  scale_y_discrete(labels = genes_of_choice, position = "right") + 
  theme_paper() + 
  theme(axis.ticks.x = element_blank(), panel.grid.major = element_blank()) + 
  xlab("") + 
  ylab("") + theme(panel.border = element_blank(), axis.line=element_blank()) + 
  scale_x_discrete(labels = genes_of_choice, position = "bottom")
ggsave("./FiguresIllustrator/FigS1/heatmap_ordered_cl30_cl31.pdf", width = 10, height = 10)

ratios %>%
  t() %>% data.frame() %>%
  add_column("Dox" = data_here$ndDox) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Dox, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", paste0("Dox (", .$Dox, " days)"))) %>%
  mutate(Condition = factor(Condition, levels = (c("Control", paste0("Dox (", c(3, 7, 14, 21), " days)"))))) %>%
  add_column(escape_status = unlist(rowData(data_here)[match(.$name, rownames(data_here)), ]$EscapeAnnotation)) %>%
  mutate(escape_status = factor(escape_status, c("silenced / variable", "facultative", "constitutive"))) %>%
  dplyr::filter(!is.na(escape_status)) %>%
  ggplot(aes(x = Condition, y = value, fill = escape_status)) + 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2)) + 
  geom_boxplot(col = "black", outlier.color = "grey") + 
  theme_paper() + ylab("ASE (d-score)") + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
  scale_fill_manual(values = c("grey", "darkgreen", "orange"), labels = c("Silenced / Variable", "Facultative", "Constitutive")) + 
  labs(fill = "Escape Category") + xlab("")
ggsave("./FiguresIllustrator/Fig1/escapee_silencing_boxplots_cl30_cl31.pdf", width = 9, height = 6)

ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  #dplyr::filter(value < 0.7) %>%
  #dplyr::filter(Clone == "E6") %>%
  mutate(Condition = factor(Condition, levels = rev(c("Control", "Dox (3d)", "Dox (7d)", "Dox (14d)", "Dox (21d)")))) %>%
  ggplot(aes(reorder(name, position), y = interaction(Clone, Condition, sep = " -- "), group = Clone, fill = value, col = value)) + geom_tile(width = 0.9) +
    scale_fill_gradient2(midpoint = 0.5, mid = "#117733", low = "#BEBEBE", high = "black", limits = c(0, 1), na.value = "white") + 
    scale_colour_gradient2(midpoint = 0.5, mid = "#117733", low = "#BEBEBE", high = "black", limits = c(0, 1), na.value = "white") + 
    scale_y_discrete(labels = genes_of_choice, position = "right") + 
    theme_paper() + 
    theme(axis.ticks.x = element_blank(), panel.grid.major = element_blank()) + 
    xlab("") + 
    ylab("") + theme(panel.border = element_blank(), axis.line=element_blank()) + 
    scale_x_discrete(labels = genes_of_choice, position = "bottom")
ggsave("./FiguresIllustrator/FigS1/heatmap_ordered_by_clones_cl30_cl31.pdf", width = 10, height = 10)

ratios %>%
  t() %>% data.frame() %>%
  add_column("Dox" = data_here$ndDox) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Dox, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  add_column(Condition = ifelse(.$Dox == 0, "Control", paste0("Dox (", .$Dox, " days)"))) %>%
  mutate(Condition = factor(Condition, levels = (c("Control", paste0("Dox (", c(3, 7, 14, 21), " days)"))))) %>%
  add_column(escape_status = unlist(rowData(data_here)[match(.$name, rownames(data_here)), ]$EscapeAnnotation)) %>%
  mutate(escape_status = factor(escape_status, c("silenced / variable", "facultative", "constitutive"))) %>%
  dplyr::filter(!is.na(escape_status)) %>%
  ggplot(aes(x = Condition, y = value, fill = Clone)) + 
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2)) + 
  geom_boxplot(col = "black", outlier.color = "grey") + 
  theme_paper() + ylab("ASE (d-score)") + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
  #scale_fill_manual(values = c("grey", "darkgreen", "orange"), labels = c("Silenced / Variable", "Facultative", "Constitutive")) + 
  labs(fill = "Escape Category") + xlab("")
ggsave("./FiguresIllustrator/FigS1/escapee_silencing_boxplots_per_clone_cl30_cl31.pdf", width = 9, height = 6)


### check whether silencing escape depends on initial escape level?

data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO", 
                                         "Aux_0_Dox_14_WO_0_WOAuxNO", "Aux_0_Dox_21_WO_0_WOAuxNO")]
data_here <- data_here[,data_here$Clone == "E6"]

ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))
ratios <- ratios[!rownames(ratios) %in% genes_out, ]

p1 <- data.frame(
  TotalExpression = rowSums(counts(data_here[rownames(ratios), ])), 
  ControlEscape = rowMeans(ratios[,data_here$ConditionClean == "Control"]), 
  d7_Escape = rowMeans(ratios[,data_here$ConditionClean == "Dox (7d)"]), 
  escape_status = rowData(data_here[rownames(ratios), ])$EscapeAnnotation
) %>%
  ggplot(aes(x = ControlEscape, d7_Escape)) + geom_point() + 
    geom_abline(linetype = 'dashed') + theme_paper() + 
    xlab("Control (escape, in / (in + act)") + ylab("7d dox (escape, in / (in + act)") + coord_fixed()
p1
ggMarginal(p1, type = "histogram")

escape_pos = start(rowRanges(data_here[rownames(ratios), ]))

data.frame(
  TotalExpression = rowSums(counts(data_here[rownames(ratios), ])), 
  ControlEscape = rowMeans(ratios[,data_here$ConditionClean == "Control"]), 
  d7_Escape = rowMeans(ratios[,data_here$ConditionClean == "Dox (7d)"]), 
  escape_status = rowData(data_here[rownames(ratios), ])$EscapeAnnotation
) %>%
  add_column(distance_to_next = unlist(lapply(escape_pos, function(x){xx = abs(escape_pos - x); return(min(xx[xx != 0]))}))) %>%
  ggplot(aes(x = cut_number(log10(distance_to_next), n = 5), ControlEscape, col = log10(distance_to_next))) + geom_boxplot(oputlier.color = NA) + 
    ggbeeswarm::geom_quasirandom() + theme_paper() + 
    xlab("Distance to nearest expressed gene") + ylab("d-score") + scale_color_viridis()

# data.frame(
#   TotalExpression = rowSums(counts(data_here[rownames(ratios), ])), 
#   ControlEscape = rowMeans(ratios[,data_here$ConditionClean == "Control"]), 
#   d7_Escape = rowMeans(ratios[,data_here$ConditionClean == "Dox (7d)"]), 
#   escape_status = rowData(data_here[rownames(ratios), ])$EscapeAnnotation
# ) %>%
#   add_column(distance_to_next = unlist(lapply(escape_pos, function(x){xx = abs(escape_pos - x); return(min(xx[xx != 0]))}))) %>%
#   ggplot(aes(x = cut_number(log10(distance_to_next), n = 5), ControlEscape, col = log10(distance_to_next))) + geom_half_point(side = "r") +
#     geom_half_violin(side = "l") + theme_paper() + 
#     xlab("Distance to nearest expressed gene") + ylab("d-score") + scale_color_viridis()

data.frame(
  TotalExpression = rowSums(counts(data_here[rownames(ratios), ])), 
  ControlEscape = rowMeans(ratios[,data_here$ConditionClean == "Control"]), 
  d7_Escape = rowMeans(ratios[,data_here$ConditionClean == "Dox (7d)"]), 
  escape_status = rowData(data_here[rownames(ratios), ])$EscapeAnnotation
) %>%
  add_column(distance_to_next = unlist(lapply(escape_pos, function(x){xx = abs(escape_pos - x); return(min(xx[xx != 0]))}))) %>%
  ggplot(aes(x = cut_number(log10(distance_to_next), n = 5), d7_Escape - ControlEscape, col = log10(distance_to_next))) + geom_boxplot(oputlier.color = NA) + 
    ggbeeswarm::geom_quasirandom() + theme_paper() + 
    xlab("Distance to nearest expressed gene") + ylab("d-score difference (7d - Control") + scale_color_viridis()

data.frame(
  TotalExpression = rowSums(counts(data_here[rownames(ratios), ])), 
  ControlEscape = rowMeans(ratios[,data_here$ConditionClean == "Control"]), 
  d7_Escape = rowMeans(ratios[,data_here$ConditionClean == "Dox (7d)"]), 
  escape_status = rowData(data_here[rownames(ratios), ])$EscapeAnnotation
) %>%
  ggplot(aes(x = cut(ControlEscape, breaks = 10), d7_Escape)) + geom_boxplot() + 
  geom_abline(linetype = 'dashed') + theme_paper() + 
  xlab("Control (escape, in / (in + act)") + ylab("7d dox (escape, in / (in + act)") + coord_fixed()

data.frame(
  TotalExpression = rowSums(counts(data_here[rownames(ratios), ])), 
  ControlEscape = rowMeans(ratios[,data_here$ConditionClean == "Control"]), 
  d7_Escape = rowMeans(ratios[,data_here$ConditionClean == "Dox (7d)"]), 
  escape_status = rowData(data_here[rownames(ratios), ])$EscapeAnnotation
) %>%
  ggplot(aes(x = ControlEscape, d7_Escape - ControlEscape)) + geom_point() + 
  theme_paper() + xlab("Control (escape, in / (in + act)") + ylab("7d dox (escape, in / (in + act)")

data.frame(
  TotalExpression = rowSums(counts(data_here[rownames(ratios), ])), 
  ControlEscape = rowMeans(ratios[,data_here$ConditionClean == "Control"]), 
  d7_Escape = rowMeans(ratios[,data_here$ConditionClean == "Dox (7d)"]), 
  escape_status = rowData(data_here[rownames(ratios), ])$EscapeAnnotation
) %>%
  ggplot(aes(x = TotalExpression, ControlEscape - d7_Escape, col = escape_status)) + geom_point() + 
  geom_abline(linetype = 'dashed') + theme_paper() + scale_x_log10()













