### We now work towards figure-level plots

### Read data
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

escapee_annotation <- readxl::read_excel("~/Desktop/Projects/XChromosome_Project/ProcessedData/ListOfEscapeeFromEdithLab.xlsx") %>%
  dplyr::select(c("ENSEMBL_v102", "final status"))
library("EnsDb.Mmusculus.v79") 
symbols <- mapIds(EnsDb.Mmusculus.v79, keys = escapee_annotation$ENSEMBL_v102, keytype = "GENEID", column="SYMBOL")
escapee_annotation$symbol <- symbols

data <- readRDS("./ProcessedData/merged_dataset.rds")
data <- data[,!data$Guide %in% c("g13", "P3") & data$Experiment != "May2022"]

# clean up names
name_conversion_condition <- list(
  "Aux_0_Dox_0_WO_NO_WOAuxNO" = "Control", 
  "Aux_0_Dox_3_WO_NO_WOAuxNO" = "Dox (3d)", 
  "Aux_0_Dox_3_WO_YES_WOAuxNO" = "Dox (3d) - washout", 
  "Aux_0_Dox_7_WO_NO_WOAuxNO" = "Dox (7d)", 
  "Aux_0_Dox_7_WO_YES_WOAuxNO" = "Dox (7d) - washout", 
  "Aux_0_Dox_7_WO_YES_WOAuxYES" = "Dox (7d) - washout (with Aux)", 
  "Aux_2_Dox_0_WO_NO_WOAuxNO" = "Aux (2d)", 
  "Aux_2_Dox_0_WO_YES_WOAuxNO" = "Aux (2d) - washout", 
  "Aux_5_Dox_0_WO_NO_WOAuxNO" = "Aux (5d)", 
  "Aux_5_Dox_3_WO_NO_WOAuxNO" = "Dox (3d), Aux (5d)", 
  "Aux_5_Dox_3_WO_YES_WOAuxNO" = "Dox (3d) - washout, Aux (5d)", 
  "Aux_9_Dox_0_WO_NO_WOAuxNO" = "Aux (9d)", 
  "Aux_9_Dox_7_WO_NO_WOAuxNO" = "Dox (7d), Aux (9d)", 
  "Aux_9_Dox_7_WO_YES_WOAuxNO" = "Dox (7d) - washout, Aux (9d)", 
  "Aux_9_Dox_7_WO_YES_WOAuxYES" = "Dox (7d) - washout (with Aux), Aux (9d)"
)

data$ConditionClean <- unlist(name_conversion_condition[data$Condition])

# remove genes with < 10 allelic reads per sample 
hist(log10((rowSums(counts_inactive(data)) + rowSums(counts_active(data))) / ncol(data) + 1), breaks = 100)
data <- data[(rowSums(counts_inactive(data)) + rowSums(counts_active(data))) / ncol(data) > 10, ]

# subset on x-linked genes

data <- data[seqnames(data) == "X", ]

### Figure 1: 
# NPC clones with higher Xist levels show less facultative escape and vice versa
# This we can't do here, bc not nearly enough N
# Possibly get the data to do it

### Figure 2: 
# Xist overexpression silences escapees in post-silencing contexts
# First, we analyze which genes escape in the baseline in our 3 different clones

data_baseline_x <- data[,data$Condition == "Aux_0_Dox_0_WO_NO_WOAuxNO"]

# Plot escape in individual clones
sample_nr = 1
data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle(colnames(data_baseline_x)[[sample_nr]])

# Prettier plot
data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle("Clone 30") + 
  xlab("Total expression (allelic reads only)") + ylab("ASE (d-score)")

# Prettier plot, with escapee-annotation
data.frame(
  total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
  ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
                                                                   rowSums(counts_active(data_baseline_x[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>% 
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ][,2])) %>% 
  mutate(escape_status = ifelse(escape_status == "NA", "S", escape_status)) %>%
  ggplot(aes(total, ase, col = escape_status)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
    theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle("Clone 30") + 
    xlab("Total expression (allelic reads only)") + ylab("ASE (d-score)") + coord_fixed(ratio = 3)
ggsave("~/Desktop/testytest.pdf", width = 20, height = 8)

# sample_nr = 2
# data.frame(
#   total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
#   ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
#                                                                    rowSums(counts_active(data_baseline_x[,sample_nr])))
# ) %>%
#   rownames_to_column("Gene") %>%
#   ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
#   theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle(colnames(data_baseline_x)[[sample_nr]])
# 
# sample_nr = 3
# data.frame(
#   total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
#   ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
#                                                                    rowSums(counts_active(data_baseline_x[,sample_nr])))
# ) %>%
#   rownames_to_column("Gene") %>%
#   ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
#   theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle(colnames(data_baseline_x)[[sample_nr]])
# 
# sample_nr = 4
# data.frame(
#   total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
#   ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
#                                                                    rowSums(counts_active(data_baseline_x[,sample_nr])))
# ) %>%
#   rownames_to_column("Gene") %>%
#   ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
#   theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle(colnames(data_baseline_x)[[sample_nr]])

# look at escapees across samples
ratios <- counts_inactive(data_baseline_x) / (counts_inactive(data_baseline_x) + counts_active(data_baseline_x))
ratios <- ratios[apply(ratios, 1, function(x){!any(is.na(x))}), ]
# remove all genes without detectable escape, 
ratios <- ratios[rownames(ratios) == "Xist" | (apply(ratios, 1, function(x){any(x > 0.1)}) & !apply(ratios, 1, function(x){any(x > 0.7)})), ]
ratios <- ratios[order(rowMeans(ratios), decreasing = T), ]

colnames(ratios) <- c("CL30", "CL31", "E6 (rep1)", "E6 (rep2)")
pheatmap::pheatmap(ratios, cluster_rows = F, cluster_cols = F, fontsize = 20, show_rownames =  F)

ratios <- counts_inactive(data_baseline_x) / (counts_inactive(data_baseline_x) + counts_active(data_baseline_x))
escapees <- apply(ratios, 2, function(x){x > 0.1})

# now we can look at the xist-oe + washout experiment
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_NO_WOAuxNO", "Aux_0_Dox_3_WO_NO_WOAuxNO", "Aux_0_Dox_3_WO_YES_WOAuxNO", 
                                         "Aux_0_Dox_7_WO_NO_WOAuxNO", "Aux_0_Dox_7_WO_YES_WOAuxNO")]

# check that xist-overexpression works
data.frame(
  Condition = data_here$ConditionClean, 
  Clone = paste0(data_here$Clone, " - ", data_here$Experiment), 
  Expression = as.numeric(counts(data_here["Xist", ]) / colSums(counts(data_here)) * 1e6)
) %>%
  ggplot(aes(x = Condition, y = Expression)) + geom_boxplot() + geom_point() + coord_flip() + 
  geom_line(aes(group = Clone, linetype = Clone, col = Clone)) + theme_paper() + 
  xlab("") + ggtitle("Xist expression (total reads)")

# check changes in escapee expression
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

escapees_per_clone <- apply(escapees, 2, function(x){rownames(escapees)[x]})
escapees_per_clone <- lapply(escapees_per_clone, function(x){x[x != "Xist"]})
names(escapees_per_clone) <- c("CL30 - December2021", "CL31 - December2021", "E6 - October2021", "E6 - March2022")

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
  scale_fill_manual(values = c("black", "red", "blue", "darkred", "darkblue"))

ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  ggplot(aes(Condition, y = reorder(name, position), fill = value)) + facet_wrap(~Clone, ncol = 4) + geom_tile() +
  theme_paper() + theme(axis.text.x=element_text(angle = 45, hjust = 1)) + scale_fill_gradientn(colors = c("blue", "grey", "red")) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + xlab("") + ylab("Chromosome position")

# ggsave("~/Desktop/test.pdf", width = 5, height = 50, limitsize = F)
  
# now quantify in a paired manner, focussing on the differences
ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  dplyr::filter(Condition %in% c("Control", "Dox (3d)", "Dox (3d) - washout") & Clone == "CL30 - December2021") %>%
  ggplot(aes(x = Condition, y = value, col = Condition, group = Condition)) + geom_point() + 
  geom_line(aes(group = name), col = "grey", alpha = 0.5) + ggpubr::stat_compare_means(paired = T) + theme_paper() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_manual(values = c("black", "red", "blue")) + xlab("") + ylab("ASE (d-score)") + 
  ggtitle("Reversibility after 3d (CL30)")

# For 7d
ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  dplyr::filter(Condition %in% c("Control", "Dox (7d)", "Dox (7d) - washout") & Clone == "CL30 - December2021") %>%
  ggplot(aes(x = Condition, y = value, col = Condition, group = Condition)) + geom_point() + 
  geom_line(aes(group = name), col = "grey", alpha = 0.5) + ggpubr::stat_compare_means(paired = T) + theme_paper() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_manual(values = c("black", "red", "blue")) + xlab("") + ylab("ASE (d-score)") + 
  ggtitle("Reversibility after 7d (CL30)")

# For 7d
ratios %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  #dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  dplyr::filter(Condition %in% c("Control", "Dox (7d)", "Dox (7d) - washout") & Clone == "CL30 - December2021") %>%
  ggplot(aes(x = Condition, y = value, col = Condition, group = Condition)) + geom_point() + 
  geom_line(aes(group = name), col = "grey", alpha = 0.5) + ggpubr::stat_compare_means(paired = T) + theme_paper() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_color_manual(values = c("black", "red", "blue")) + xlab("") + ylab("ASE (d-score)") + 
  ggtitle("Reversibility after 7d (CL30)")

# For all samples normalize to untreated, and collapse to mean / median
ratios_clone1 <- ratios[,paste0(colData(data)[colnames(ratios), ]$Clone, "_", colData(data)[colnames(ratios), ]$Experiment) == "CL30_December2021"]
ratios_clone1 <- ratios_clone1 - ratios_clone1[,"CL30.NoDoxNoAux" ]
ratios_clone1 <- ratios_clone1 %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here[,rownames(.)]$ConditionClean) %>%
  add_column("Clone" = paste0(data_here[,rownames(.)]$Clone, " - ", data_here[,rownames(.)]$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  group_by(Condition) %>%
  mutate(mean_escape = median(value, na.rm = T))

ratios_clone2 <- ratios[,paste0(colData(data)[colnames(ratios), ]$Clone, "_", colData(data)[colnames(ratios), ]$Experiment) == "CL31_December2021"]
ratios_clone2 <- ratios_clone2 - ratios_clone2[,"CL31.NoDoxNoAux"]
ratios_clone2 <- ratios_clone2  %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here[,rownames(.)]$ConditionClean) %>%
  add_column("Clone" = paste0(data_here[,rownames(.)]$Clone, " - ", data_here[,rownames(.)]$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  group_by(Condition) %>%
  mutate(mean_escape = median(value, na.rm = T))

ratios_clone3 <- ratios[,paste0(colData(data)[colnames(ratios), ]$Clone, "_", colData(data)[colnames(ratios), ]$Experiment) == "E6_October2021"]
ratios_clone3 <- ratios_clone3 - ratios_clone3[,"C.E6.NoDox"]
ratios_clone3 <- ratios_clone3 %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here[,rownames(.)]$ConditionClean) %>%
  add_column("Clone" = paste0(data_here[,rownames(.)]$Clone, " - ", data_here[,rownames(.)]$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  group_by(Condition) %>%
  mutate(mean_escape = median(value, na.rm = T)) %>%
  ungroup()

ratios_clone4 <- ratios[,paste0(colData(data)[colnames(ratios), ]$Clone, "_", colData(data)[colnames(ratios), ]$Experiment) == "E6_March2022"]
ratios_clone4 <- ratios_clone4 - ratios_clone4[,"C.E6_NoDox_r2"]
ratios_clone4 <- ratios_clone4 %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here[,rownames(.)]$ConditionClean) %>%
  add_column("Clone" = paste0(data_here[,rownames(.)]$Clone, " - ", data_here[,rownames(.)]$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  group_by(Condition) %>%
  mutate(mean_escape = median(value, na.rm = T)) %>%
  ungroup()

do.call("rbind", list(ratios_clone1, ratios_clone2, ratios_clone3, ratios_clone4)) %>%
  ggplot(aes(x = Condition, y = mean_escape, col = Clone)) + geom_point() + geom_line(aes(group = Clone), linetype= 'dashed') + 
  theme_paper() + xlab("") + ylab("Difference in ASE (d-score)") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# From the previous plots, it looks like genes react differently to Xist-OE and that silencing is reversible to different extents
# we therefore assess gene-specific effects of dox + washout using linear modelling
library(VGAM)

data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_NO_WOAuxNO", "Aux_0_Dox_3_WO_NO_WOAuxNO", "Aux_0_Dox_3_WO_YES_WOAuxNO", 
                                         "Aux_0_Dox_7_WO_NO_WOAuxNO", "Aux_0_Dox_7_WO_YES_WOAuxNO")]

# define interesting genes: genes, where dox explains something
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

data_test_inactive <- counts_inactive(data_here)
data_test_active <- counts_active(data_here)

fitting_metadata_here <- colData(data_here)

all_coefs_anova <- data.frame(do.call("rbind", lapply(1:nrow(data_test_active), function(i){
  tryCatch({
    y = as.numeric(data_test_inactive[i, ])
    N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
    fit1 <- glm(cbind(y, N-y) ~ Clone + ndDox + washout, family = "binomial", data = fitting_metadata_here)
    anova_dev = data.frame(anova(fit1))[,2][-1]
    resid_dev = anova(fit1)[,4][1] - sum(anova_dev)
    return(c(anova_dev, resid_dev) / sum(c(anova_dev, resid_dev)))
  }, error = function(cond) {
    return(c(NA))
  })
})))

rownames(all_coefs_anova) <- rownames(data_test_inactive)
colnames(all_coefs_anova) <- c("Clone", "Dox", "Washout", "Residual")

all_coefs_anova %>%
  pivot_longer(-c()) %>%
  ggplot(aes(x = factor(name, levels = c("Clone", "Dox", "Washout", "Residual")), y = value)) + geom_violin(fill = "grey") + 
  geom_boxplot(width = 0.05, outlier.color = NA) + 
  xlab("Component") + ylab("Proportion of Deviance") + theme_paper()

head(all_coefs_anova[order(all_coefs_anova$Residual, decreasing = T), ], n = 20)

# define escapees here
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))
escapees <- apply(ratios, 2, function(x){x > 0.05 & x < 0.7})
escapees[is.na(escapees)] <- FALSE

general_escapees <- data_here[names(which(apply(escapees[,data_here$Condition == "Aux_0_Dox_0_WO_NO_WOAuxNO"], 1, function(x){all(x)}))), ]

data_test_inactive <- counts_inactive(general_escapees)
data_test_active <- counts_active(general_escapees)

fitting_metadata_here <- colData(general_escapees)

# 
all_coefs <- data.frame(do.call("rbind", lapply(1:nrow(data_test_active), function(i){
  tryCatch({
    y = as.numeric(data_test_inactive[i, ])
    N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
    # fit1 <- glm(cbind(y, N-y) ~ Clone + ndDox + washout, family = "binomial", data = fitting_metadata_here)
    fit1 <- glm(cbind(y, N-y) ~ Clone * ndDox + washout, family = "binomial", data = fitting_metadata_here)
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
  coord_flip() + theme_paper() + xlab("Covariate") +
  theme(aspect.ratio = 2)

all_coefs %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-c(Gene)) %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ][,2])) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(x = name, y = value, col = escape_status)) + 
    geom_boxplot() + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.2) +
    coord_flip() + theme_paper() + xlab("Covariate") +
    theme(aspect.ratio = 2) + theme(legend.position = "top")

# look at coefficients between washout and dox silencing -- looks cool! 
all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ][,2])) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(ndDox3, washoutYES, col = escape_status)) + geom_point() + 
    geom_smooth(method = "lm", col = "black", linetype = 'dashed', size = 0.5) + 
    ggrepel::geom_text_repel(aes(label = Gene)) + 
    geom_abline(linetype = "dashed", slope = -1) + 
    theme_paper()

all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ][,2])) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(ndDox7, washoutYES, col = escape_status)) + geom_point() + 
  geom_smooth(method = "lm", col = "black", linetype = 'dashed', size = 0.5) + 
  ggrepel::geom_text_repel(aes(label = Gene)) + 
  geom_abline(linetype = "dashed", slope = -1) + 
  theme_paper()

# example genes

plot_gene <- function(gene){
  i = gene
  y = as.numeric(data_test_inactive[i, ])
  N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
  fit1 <- glm(cbind(y, N-y) ~ Clone + ndDox + washout, family = "binomial", data = fitting_metadata_here)
  fit0 <- glm(cbind(y, N-y) ~ 1, family = "binomial", data = fitting_metadata_here)
  
  data.frame(
    y = y,
    N = N,
    Clone = paste0(data_here$Clone, " - ", data_here$Experiment),
    Cov = data_here$ConditionClean
  ) %>% ggplot(aes(Cov, y / N, col = Clone)) + geom_point(size = 4) + coord_flip() + ylim(0, 1) +
    ggtitle(i) +
    geom_hline(yintercept = rev_logit(coef(fit1)[[1]]), linetype = "dashed") + 
    geom_line(aes(group = Clone), linetype = "dashed") + 
    theme_paper() + xlab("") + ylab("ASE (d-score)")
}

plot_gene("Fmr1")
plot_gene("Plxnb3")
plot_gene("Jpx")
plot_gene("Eif2s3x")

antonia_genes <- c("Plxnb3", "Slc6a8", "Bcap31", "Idh3g", "Ssr4", "Hcfc1", "Mecp2", "G6pdx")

plot_gene("Plxnb3")
plot_gene("Slc6a8")
plot_gene("Bcap31")
plot_gene("Idh3g")
plot_gene("Ssr4")
plot_gene("Hcfc1")
plot_gene("Mecp2")
plot_gene("G6pdx")

## 


##################################################################################################

general_escapees <- data_here[names(which(apply(
  escapees[,data_here$Condition == "Aux_0_Dox_0_WO_NO_WOAuxNO" & data_here$Clone == "E6"], 
  1, function(x){all(x)}))), ]

data_test_inactive <- counts_inactive(general_escapees[,data_here$Clone == "E6"])
data_test_active <- counts_active(general_escapees[,data_here$Clone == "E6"])

fitting_metadata_here <- colData(general_escapees)[data_here$Clone == "E6", ]

genes_out <- c("Car5b", "Apoo", "2210013O21Rik")
data_test_inactive <- data_test_inactive[!rownames(data_test_inactive) %in% genes_out, ]
data_test_active <- data_test_active[!rownames(data_test_active) %in% genes_out, ]

all_coefs <- data.frame(do.call("rbind", lapply(1:nrow(data_test_active), function(i){
  tryCatch({
    y = as.numeric(data_test_inactive[i, ])
    N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
    fit1 <- glm(cbind(y, N-y) ~ ndDox + washout, family = "binomial", data = fitting_metadata_here)
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
all_coefs_silencing <- all_coefs

# all_coefs[abs(all_coefs) > 10] <- sign(all_coefs[abs(all_coefs) > 10]) * 10

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
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ][,2])) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(x = name, y = value, col = escape_status)) + 
  geom_boxplot() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.2) +
  coord_flip() + theme_paper() + xlab("Covariate") +
  theme(aspect.ratio = 2) + theme(legend.position = "top")

# look at coefficients between washout and dox silencing -- looks cool! 
all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ][,2])) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(ndDox3, washoutYES, col = escape_status)) + geom_point() + 
  geom_smooth(method = "lm", col = "black", linetype = 'dashed', size = 0.5) + 
  ggrepel::geom_text_repel(aes(label = Gene)) + 
  geom_abline(linetype = "dashed", slope = -1) + 
  theme_paper()

all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_status = unlist(escapee_annotation[match(.$Gene, escapee_annotation$symbol), ][,2])) %>% 
  dplyr::filter(escape_status != "NA") %>%
  ggplot(aes(ndDox7, washoutYES, col = escape_status)) + geom_point() + 
  geom_smooth(method = "lm", col = "black", linetype = 'dashed', size = 0.5) + 
  ggrepel::geom_text_repel(aes(label = Gene)) + 
  geom_abline(linetype = "dashed", slope = -1) + 
  theme_paper()

# example genes

plot_gene <- function(gene){
  i = gene
  y = as.numeric(data_test_inactive[i, ])
  N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
  fit1 <- glm(cbind(y, N-y) ~ ndDox + washout, family = "binomial", data = fitting_metadata_here)

  data.frame(
    y = y,
    N = N,
    Clone = paste0(data_here$Clone, " - ", data_here$Experiment),
    Cov = data_here$ConditionClean
  ) %>% dplyr::filter(Clone %in% c("E6 - March2022", "E6 - October2021")) %>%
    ggplot(aes(Cov, y / N, col = Clone)) + geom_point(size = 4) + coord_flip() + ylim(0, 1) +
    ggtitle(i) +
    geom_hline(yintercept = rev_logit(coef(fit1)[[1]]), linetype = "dashed") + 
    geom_line(aes(group = Clone), linetype = "dashed") + 
    theme_paper() + xlab("") + ylab("ASE (d-score)")
}

plot_gene("Prdx4")
plot_gene("Jpx")
plot_gene("Rbbp7")
plot_gene("Kdm5c")

# look at "variably" escaping genes



## look for clustered / non-clustered responses

all_coefs_mod <- all_coefs %>%
  add_column(chromosomal_position = start(rowRanges(general_escapees)[rownames(.), ]))

# somewhat ad hoc: median distance to the closest k escapees
all_coefs_mod$closest_escapee <- unlist(lapply(rownames(all_coefs_mod), function(gene){
  distances = abs(all_coefs_mod$chromosomal_position - all_coefs_mod[gene, ]$chromosomal_position)
  distances = sort(distances, decreasing = F)[1:6] # k = 5
  min(distances[distances > 0])
}))
  
all_coefs_mod %>% 
  add_column(reversibility = abs(all_coefs_mod$ndDox7) - all_coefs_mod$washoutYES) %>%
  pivot_longer(c("ndDox3", "ndDox7", "washoutYES", "reversibility")) %>%
  ggplot(aes(x = closest_escapee, y = value)) + 
    geom_point() + scale_x_log10() + geom_smooth(method = "lm") + 
    facet_wrap(~name, scales = "free") + 
    ggpubr::stat_cor(label.x = 3, label.y = 6, method = 'pearson')


#############################################################################################################################
# now we check the genomic distributions of coefficients across the linear genome

ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("Condition" = data_here$ConditionClean) %>%
  add_column("Clone" = paste0(data_here$Clone, " - ", data_here$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  add_column(position = start(rowRanges(data_here[as.character(.$name), ]))) %>%
  ggplot(aes(Condition, y = reorder(name, position), fill = value)) + facet_wrap(~Clone, ncol = 4) + geom_tile() +
  theme_paper() + theme(axis.text.x=element_text(angle = 45, hjust = 1)) + scale_fill_gradientn(colors = c("blue", "grey", "red")) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + xlab("") + ylab("Chromosome position")

all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_ratio = abs(.$ndDox7 / .$washoutYES)) %>%
  dplyr::filter(!escape_ratio < -2) %>%
  add_column(position = start(rowRanges(data_here[as.character(.$Gene), ]))) %>%
  pivot_longer(-c(Gene, position)) %>%
  #dplyr::filter(name == "escape_ratio") %>%
  ggplot(aes(x = reorder(Gene, position), y = name, fill = value)) + geom_tile() + 
  theme_paper() + scale_fill_gradient2(low = "blue", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5))

all_coefs %>%
  rownames_to_column("Gene") %>%
  add_column(escape_ratio = .$washoutYES / .$ndDox7) %>%
  dplyr::filter(!escape_ratio < -2) %>%
  add_column(position = start(rowRanges(data_here[as.character(.$Gene), ]))) %>%
  pivot_longer(-c(Gene, position)) %>%
  ggplot(aes(x = position, y = name, col = value)) + geom_point(size = 5) + 
  theme_paper() + scale_color_gradient2(low = "blue", high = "red") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5))

# look for associated features

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
  
  data_fit <- data_set_fit_scaled
  
  correlations <- apply(data_fit, 2, function(x){cor(coef_vector, as.numeric(x), method = "spearman")})
  
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

coefs_fit_here$ratio <- abs(coefs_fit_here$ndDox7 / coefs_fit_here$washoutYES)

rf_results <- apply(coefs_fit_here, 2, run_rf_per_vector)

rsqs <- unlist(lapply(rf_results, function(x){mean(x[[2]])}))
mses <- unlist(lapply(rf_results, function(x){mean(x[[3]])}))

i = 1
colnames(coefs_fit_here)[i]
testy <- rf_results[[i]][[1]]

# combine correlations
all_cors_heatmap <- do.call("rbind", lapply(rf_results, function(x){x[[1]][,"cor"]}))
colnames(all_cors_heatmap) <- colnames(data_set_fit)

pheatmap::pheatmap(all_cors_heatmap, cluster_rows = F, color=colorRampPalette(c("navy", "white", "red"))(50))

all_cors_heatmap_transposed <- t(all_cors_heatmap)
head(all_cors_heatmap_transposed[order(abs(all_cors_heatmap_transposed[,5]), decreasing = T), ], n = 20)

data.frame(
  gene = rownames(coefs_fit_here), 
  ratio = coefs_fit_here$ratio, 
  feature = data_set_fit$gene_density
) %>%
  ggplot(aes(x = feature, y = ratio)) + geom_point() + geom_smooth(method = "lm") + ggrepel::geom_text_repel(aes(label = gene))

## look at ratio predictivity

data_fit <- data_set_fit_scaled
names_save <- colnames(data_fit)
colnames(data_fit) <- paste0("Feature", 1:ncol(data_fit))

data_fit <- cbind(data_fit, "objective" = coefs_fit_here$ratio)

rf <- randomForest(
  objective ~ .,
  data=data_fit,
  importance = T, 
)

importance_df <- data.frame(rf$importance)
rownames(importance_df) <- names_save

######################################################################################################################################################

### Figure 3: 
# Xist silencing in post-silencing contexts is spen-dependent
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_NO_WOAuxNO", "Aux_0_Dox_3_WO_NO_WOAuxNO", "Aux_0_Dox_7_WO_NO_WOAuxNO", 
                                         "Aux_2_Dox_0_WO_NO_WOAuxNO", "Aux_5_Dox_0_WO_NO_WOAuxNO", "Aux_9_Dox_0_WO_NO_WOAuxNO", 
                                         "Aux_5_Dox_3_WO_NO_WOAuxNO", "Aux_9_Dox_7_WO_NO_WOAuxNO")]

data_here <- data_here[,!data_here$Clone == "E6"]

# check that xist-overexpression works
data.frame(
  Condition = factor(data_here$ConditionClean, levels = c("Control", "Aux (2d)", "Aux (5d)", "Aux (9d)", "Dox (3d)", 
                                                          "Dox (3d), Aux (5d)", "Dox (7d)", "Dox (7d), Aux (9d)")), 
  Clone = paste0(data_here$Clone, " - ", data_here$Experiment), 
  Expression = as.numeric(counts(data_here["Xist", ]) / colSums(counts(data_here)) * 1e6)
) %>%
  ggplot(aes(x = Condition, y = Expression)) + geom_boxplot() + geom_point() + coord_flip() + 
  geom_line(aes(group = Clone, linetype = Clone, col = Clone)) + theme_paper() + 
  xlab("") + ggtitle("Xist expression (total reads)")

# check changes in escapee expression
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

escapees_per_clone <- apply(escapees, 2, function(x){rownames(escapees)[x]})
escapees_per_clone <- lapply(escapees_per_clone, function(x){x[x != "Xist"]})
names(escapees_per_clone) <- c("CL30 - December2021", "CL31 - December2021", "E6 - October2021", "E6 - March2022")

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
  scale_fill_manual(values = c("black", "lightgrey", "grey", "darkgrey", "red", "blue", "darkred", "darkblue"))

ratios_clone1 <- ratios[,paste0(colData(data)[colnames(ratios), ]$Clone, "_", colData(data)[colnames(ratios), ]$Experiment) == "CL30_December2021"]
ratios_clone1 <- ratios_clone1 - ratios_clone1[,"CL30.NoDoxNoAux" ]
ratios_clone1 <- ratios_clone1 %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here[,rownames(.)]$ConditionClean) %>%
  add_column("Clone" = paste0(data_here[,rownames(.)]$Clone, " - ", data_here[,rownames(.)]$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  group_by(Condition) %>%
  mutate(mean_escape = median(value))

ratios_clone2 <- ratios[,paste0(colData(data)[colnames(ratios), ]$Clone, "_", colData(data)[colnames(ratios), ]$Experiment) == "CL31_December2021"]
ratios_clone2 <- ratios_clone2 - ratios_clone2[,"CL31.NoDoxNoAux"]
ratios_clone2 <- ratios_clone2  %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here[,rownames(.)]$ConditionClean) %>%
  add_column("Clone" = paste0(data_here[,rownames(.)]$Clone, " - ", data_here[,rownames(.)]$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  group_by(Condition) %>%
  mutate(mean_escape = median(value))

levels_here <- c("Control", "Aux (2d)", "Aux (5d)", "Aux (9d)", "Dox (3d)", "Dox (3d), Aux (5d)", "Dox (7d)", "Dox (7d), Aux (9d)")
do.call("rbind", list(ratios_clone1, ratios_clone2)) %>%
  ggplot(aes(x = factor(Condition, levels = levels_here), y = mean_escape, col = Clone)) + geom_point() + coord_flip() + geom_line(aes(group = Clone), linetype= 'dashed') + 
  theme_paper() + xlab("Condition") + ylab("Average escape difference (to control)")

# If we look at all genes (not just escapees), we see that spen desilences genes to an extent: 
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

ratios_clone1 <- ratios[,paste0(colData(data)[colnames(ratios), ]$Clone, "_", colData(data)[colnames(ratios), ]$Experiment) == "CL30_December2021"]
ratios_clone1 <- ratios_clone1 - ratios_clone1[,"CL30.NoDoxNoAux" ]
ratios_clone1 <- ratios_clone1 %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here[,rownames(.)]$ConditionClean) %>%
  add_column("Clone" = paste0(data_here[,rownames(.)]$Clone, " - ", data_here[,rownames(.)]$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Condition) %>%
  mutate(mean_escape = mean(value, na.rm = T))

ratios_clone2 <- ratios[,paste0(colData(data)[colnames(ratios), ]$Clone, "_", colData(data)[colnames(ratios), ]$Experiment) == "CL31_December2021"]
ratios_clone2 <- ratios_clone2 - ratios_clone2[,"CL31.NoDoxNoAux"]
ratios_clone2 <- ratios_clone2  %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here[,rownames(.)]$ConditionClean) %>%
  add_column("Clone" = paste0(data_here[,rownames(.)]$Clone, " - ", data_here[,rownames(.)]$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Condition) %>%
  mutate(mean_escape = mean(value, na.rm = T))

levels_here <- c("Control", "Aux (2d)", "Aux (5d)", "Aux (9d)", "Dox (3d)", "Dox (3d), Aux (5d)", "Dox (7d)", "Dox (7d), Aux (9d)")
do.call("rbind", list(ratios_clone1, ratios_clone2)) %>%
  ggplot(aes(x = factor(Condition, levels = levels_here), y = mean_escape, col = Clone)) + geom_point() + coord_flip() + geom_line(aes(group = Clone), linetype= 'dashed') + 
  theme_paper() + xlab("Condition") + ylab("Average escape difference (to control)")


# Is the re-establishment of expression spen-dependent / does depletion of spen accelerate re-activation? 
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_NO_WOAuxNO", "Aux_0_Dox_7_WO_NO_WOAuxNO", 
                                         "Aux_0_Dox_7_WO_YES_WOAuxYES", "Aux_0_Dox_7_WO_YES_WOAuxNO")]

# check changes in escapee expression
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

ratios_clone1 <- ratios[,paste0(colData(data)[colnames(ratios), ]$Clone, "_", colData(data)[colnames(ratios), ]$Experiment) == "CL30_December2021"]
ratios_clone1 <- ratios_clone1 - ratios_clone1[,"CL30.NoDoxNoAux" ]
ratios_clone1 <- ratios_clone1 %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here[,rownames(.)]$ConditionClean) %>%
  add_column("Clone" = paste0(data_here[,rownames(.)]$Clone, " - ", data_here[,rownames(.)]$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  group_by(Condition) %>%
  mutate(mean_escape = mean(value, na.rm = T))

ratios_clone2 <- ratios[,paste0(colData(data)[colnames(ratios), ]$Clone, "_", colData(data)[colnames(ratios), ]$Experiment) == "CL31_December2021"]
ratios_clone2 <- ratios_clone2 - ratios_clone2[,"CL31.NoDoxNoAux"]
ratios_clone2 <- ratios_clone2  %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here[,rownames(.)]$ConditionClean) %>%
  add_column("Clone" = paste0(data_here[,rownames(.)]$Clone, " - ", data_here[,rownames(.)]$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  group_by(Condition) %>%
  mutate(mean_escape = mean(value, na.rm = T))

levels_here <- c("Control", "Dox (7d)", "Dox (7d) - washout", "Dox (7d) - washout (with Aux)")
do.call("rbind", list(ratios_clone1, ratios_clone2)) %>%
  ggplot(aes(x = factor(Condition, levels = levels_here), y = mean_escape, col = Clone)) + geom_point() + coord_flip() + geom_line(aes(group = Clone), linetype= 'dashed') + 
  theme_paper() + xlab("Condition") + ylab("Average escape difference (to control)")

# maybe there is some effect, but quite minimal

### Figure 4: 
# Xist-depletion leads to upregulation of escapees / X-linked genes in post-silencing contexts

data <- readRDS("./ProcessedData/merged_dataset.rds")
data_here <- data[,data$Guide %in% c("g13", "P3")]
data_here <- data_here[as.character(seqnames(data_here)) == "X", ]
data_here$ConditionClean <- ifelse(data_here$Guide == "P3", "non-targeting control", "sgXist")

data_here <- data_here[rowSums(counts_inactive(data_here) + counts_active(data_here)) > 200, ]

# check that xist-knockdown works
data.frame(
  Condition = data_here$ConditionClean, 
  Clone = paste0(data_here$Clone, " - ", data_here$Experiment), 
  Xist_Expression = as.numeric(counts(data_here["Xist", ]) / colSums(counts(data_here)) * 1e6)
) %>%
  ggplot(aes(x = factor(Condition, levels = c("sgXist", "non-targeting control")),  y = Xist_Expression)) + geom_boxplot() + geom_point() + 
  geom_line(aes(group = Clone, linetype = Clone, col = Clone)) + theme_paper() + 
  coord_flip() + xlab("") + ylab("Xist (normalized counts)")

data_baseline_x <- data_here

# sample_nr = 1
# data.frame(
#   total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
#   ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
#                                                                    rowSums(counts_active(data_baseline_x[,sample_nr])))
# ) %>%
#   rownames_to_column("Gene") %>%
#   ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
#   theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle(colnames(data_baseline_x)[[sample_nr]])
# 
# sample_nr = 2
# data.frame(
#   total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
#   ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
#                                                                    rowSums(counts_active(data_baseline_x[,sample_nr])))
# ) %>%
#   rownames_to_column("Gene") %>%
#   ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
#   theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle(colnames(data_baseline_x)[[sample_nr]])
# 
# sample_nr = 3
# data.frame(
#   total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
#   ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
#                                                                    rowSums(counts_active(data_baseline_x[,sample_nr])))
# ) %>%
#   rownames_to_column("Gene") %>%
#   ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
#   theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle(colnames(data_baseline_x)[[sample_nr]])
# 
# sample_nr = 4
# data.frame(
#   total = rowSums(counts_inactive(data_baseline_x[,sample_nr])) + rowSums(counts_active(data_baseline_x[,sample_nr])),
#   ase = rowSums(counts_inactive(data_baseline_x[,sample_nr])) / (rowSums(counts_inactive(data_baseline_x[,sample_nr])) +
#                                                                    rowSums(counts_active(data_baseline_x[,sample_nr])))
# ) %>%
#   rownames_to_column("Gene") %>%
#   ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
#   theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle(colnames(data_baseline_x)[[sample_nr]])

# 
data_baseline_x <- data_here[,data_here$Guide == "P3"]

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
  scale_fill_manual(values = c("black", "red"))

# For all samples normalize to untreated, and collapse to mean / median
ratios_clone1 <- ratios[,paste0(colData(data)[colnames(ratios), ]$Clone, "_", colData(data)[colnames(ratios), ]$Experiment) == "C3_March2022"]
ratios_clone1 <- ratios_clone1 - ratios_clone1$C.C12C3.KRAB_P3_r1
ratios_clone1 <- ratios_clone1 %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here[,rownames(.)]$ConditionClean) %>%
  add_column("Clone" = paste0(data_here[,rownames(.)]$Clone, " - ", data_here[,rownames(.)]$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  # dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  group_by(Condition) %>%
  mutate(mean_escape = mean(value, na.rm = T))

ratios_clone2 <- ratios[,paste0(colData(data)[colnames(ratios), ]$Clone, "_", colData(data)[colnames(ratios), ]$Experiment) == "D11_March2022"]
ratios_clone2 <- ratios_clone2 - ratios_clone2$C.D11.KRAB_P3_r1
ratios_clone2 <- ratios_clone2  %>%
  t() %>% data.frame() %>%
  add_column("Condition" = data_here[,rownames(.)]$ConditionClean) %>%
  add_column("Clone" = paste0(data_here[,rownames(.)]$Clone, " - ", data_here[,rownames(.)]$Experiment)) %>%
  pivot_longer(-c(Condition, Clone)) %>%
  group_by(Clone) %>%
  # dplyr::filter(name %in% escapees_per_clone[[cur_group_id()]]) %>%
  group_by(Condition) %>%
  mutate(mean_escape = mean(value, na.rm = T))

do.call("rbind", list(ratios_clone1, ratios_clone2)) %>%
  ggplot(aes(x = factor(Condition, levels = c("sgXist", "non-targeting control")), y = mean_escape, col = Clone)) + geom_point() + coord_flip() + geom_line(aes(group = Clone), linetype= 'dashed') + 
  theme_paper() + xlab("") + ylab("Average escape difference (to control)")

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
  ggtitle("Clone22")

# Now we use LMs to check activation effects
library(VGAM)

# define escapees here
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))

data_test_inactive <- counts_inactive(data_here)
data_test_active <- counts_active(data_here)

fitting_metadata_here <- colData(data_here)

all_coefs <- data.frame(do.call("rbind", lapply(1:nrow(data_test_active), function(i){
  tryCatch({
    y = as.numeric(data_test_inactive[i, ])
    N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
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

row_annotation <- data.frame(
  Gene = rownames(all_coefs), 
  escapee_annotation = escapee_annotation[match(rownames(all_coefs), escapee_annotation$symbol), ][,2]
) %>% column_to_rownames("Gene")

pheatmap::pheatmap(all_coefs, cluster_cols = F, color=colorRampPalette(c("navy", "white", "red"))(50), annotation_row = row_annotation, fontsize = 7)

head(all_coefs[order(all_coefs$GuideP3, decreasing = F), ], n = 20)

all_coefs %>%
  rownames_to_column("Gene") %>% 
  arrange(GuideP3) %>%
  ggplot(aes(x = 1:nrow(all_coefs), y = GuideP3)) + geom_point() + 
  ggrepel::geom_text_repel(aes(label = Gene), max.overlaps = 30)

# example genes
i = "Mid1ip1"
y = as.numeric(data_test_inactive[i, ])
N = as.numeric(data_test_active[i, ] + data_test_inactive[i, ])
fit1 <- glm(cbind(y, N-y) ~ Clone + Guide, family = "binomial", data = fitting_metadata_here)

summary(fit1)
coef(fit1)

data.frame(coef(fit1)) %>%
  rownames_to_column("Variable") %>%
  ggplot(aes(Variable, coef.fit1.)) + geom_bar(stat = "identity") +
  coord_flip()

data.frame(
  y = y,
  N = N,
  Cov = data_here$Guide,
  Clone = data_here$Clone
) %>% ggplot(aes(Cov, y / N, col = Clone)) + geom_point() + coord_flip() + ylim(0, 1) +
  ggtitle(i) +
  geom_hline(yintercept = rev_logit(coef(fit1)[[1]]), linetype = "dashed")

# Is there a relationship between baseline escape and increase?
escapee_annotation <- escapee_annotation %>%
  filter(!duplicated(symbol)) %>%
  filter(!is.na(symbol)) %>%
  column_to_rownames("symbol")

ratios[,1:2] %>%
  rownames_to_column("Gene") %>%
  add_column(anno = escapee_annotation[.$Gene, ]$`final status`) %>%
  ggplot(aes(C.C12C3.KRAB_P3_r1, C.C12C3.KRAB_g13_r1, col = anno)) + geom_point() + 
    geom_abline(linetype = "dashed") + 
  theme_paper() + ggrepel::geom_text_repel(aes(label = Gene)) + 
  ggtitle("Clone 3") + xlab("NTC") + ylab("sgXist") + coord_fixed()

ratios[,1:2] %>%
  rownames_to_column("Gene") %>%
  add_column(anno = escapee_annotation[.$Gene, ]$`final status`) %>%
  ggplot(aes(x = anno, C.C12C3.KRAB_g13_r1 - C.C12C3.KRAB_P3_r1, col = anno)) +
  geom_boxplot(width = 0.2, outlier.color = NA) + ggbeeswarm::geom_quasirandom() + theme_paper() + coord_flip() + 
  ylab("Change in escape (d_NTC - d_sgXist)") + xlab("")

ratios[,3:4] %>%
  rownames_to_column("Gene") %>%
  add_column(anno = escapee_annotation[.$Gene, ]$`final status`) %>%
  ggplot(aes(C.D11.KRAB_P3_r1, C.D11.KRAB_g13_r1, col = anno)) + geom_point() + 
  geom_abline(linetype = "dashed") + 
  theme_paper() + ggrepel::geom_text_repel(aes(label = Gene)) + 
  ggtitle("Clone 22") + xlab("NTC") + ylab("sgXist") + coord_fixed()

ratios[,3:4] %>%
  rownames_to_column("Gene") %>%
  add_column(anno = escapee_annotation[.$Gene, ]$`final status`) %>%
  ggplot(aes(x = anno, C.D11.KRAB_g13_r1 - C.D11.KRAB_P3_r1, col = anno)) +
  geom_boxplot(width = 0.2, outlier.color = NA) + ggbeeswarm::geom_quasirandom() + theme_paper() + coord_flip() + 
  ylab("Change in escape (d_NTC - d_sgXist)") + xlab("")

ratios %>%
  ggplot(aes(x = C.D11.KRAB_g13_r1 - C.D11.KRAB_P3_r1,
             y = C.C12C3.KRAB_g13_r1 - C.C12C3.KRAB_P3_r1)) + geom_point() + theme_paper() + 
  geom_smooth(method = "lm")

# what's the relationship between responsiveness to Xist-OE and knockdown?

xist_oe_response <- rowMeans(ratios[,c(1, 3)]) - rowMeans(ratios[,c(2, 4)])
xist_oe_response <- xist_oe_response[rowMeans(ratios[,c(2, 4)]) > 0.05]

both.genes <- intersect(names(xist_oe_response), rownames(all_coefs_silencing))

data.frame(
  gene = both.genes, 
  xist_kd = xist_oe_response[both.genes], 
  xist_oe = all_coefs_silencing[both.genes, ]$ndDox7
) %>% 
  dplyr::filter(xist_kd > -20) %>%
  ggplot(aes(xist_kd, xist_oe)) + geom_point() + 
    ggrepel::geom_text_repel(aes(label = gene)) + 
    theme_paper() + xlab("Xist knockdown") + ylab("Xist overexpression") + 
  geom_smooth(method = "lm") + ggpubr::stat_cor(label.x = -0.1, label.y = -7.5)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Look for per clone / Xist effects
##### ##### ##### ##### ##### ##### ##### ##### ##### #####

data <- readRDS("./ProcessedData/merged_dataset.rds")
data <- logNormCounts(data)
data <- data[, data$Experiment == "May2022" | data$Condition == "Aux_0_Dox_0_WO_NO_WOAuxNO" | grepl("P3", data$CorrectedSampleName)]
data <- data[seqnames(data) == "X", ]
ratios <- counts_inactive(data) / (counts_inactive(data) + counts_active(data))
ratios <- ratios[ rowMeans(counts_inactive(data) + counts_active(data)) > 10, ]

escapees_per_sample <- apply(ratios, 2, function(x){x = rownames(ratios)[x < 0.7 & x > 0.1]; x[!is.na(x)]})
all_escapees <- Reduce("union", escapees_per_sample)
ratios_escapee <- ratios[c(all_escapees, "Xist"), ]

ComplexHeatmap::Heatmap(ratios_escapee[order(rowMeans(ratios_escapee), decreasing = T), ], 
                        cluster_rows = F, cluster_columns = F, row_names_gp = grid::gpar(fontsize = 4), 
                        circlize::colorRamp2(c(0, 0.5, 1), c("blue", "lightyellow", "red")))

antonia_genes <- c("Plxnb3", "Slc6a8", "Bcap31", "Idh3g", "Ssr4", "Hcfc1", "Mecp2", "G6pdx")

ratios %>%
  t() %>% data.frame(check.names = F) %>%
  add_column("SampleName" = data$SampleName) %>%
  pivot_longer(-c(SampleName)) %>%
  add_column(position = start(rowRanges(data[as.character(.$name), ]))) %>%
  ggplot(aes(SampleName, y = reorder(name, position), fill = value)) + geom_tile() +
    theme_paper() + theme(axis.text.x=element_text(angle = 45, hjust = 1)) + scale_fill_gradientn(colors = c("blue", "grey", "red")) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + xlab("") + ylab("Chromosome position")

diffd_with_xist <- list(
  "B1" = "no_xist", 
  "C3" = "unclear", 
  "C5" = "no_xist", 
  "CL30" = "with_xist", 
  "CL31" = "with_xist", 
  "D11" = "unclear", 
  "E6" = "no_xist", 
  "H4" = "no_xist", 
  "JTG" = "unclear"
)

data.frame(
  Clone = data$Clone, 
  xist_levels = logcounts(data)["Xist", ], 
  average_escape = colMeans(ratios, na.rm = T), 
  with_xist = unlist(diffd_with_xist[data$Clone])
) %>%
  ggplot(aes(x = xist_levels, y = average_escape, col = with_xist)) + geom_point(size = 5) + 
  theme_classic() + geom_smooth(method = "lm") + ggrepel::geom_text_repel(aes(label = Clone))

data.frame(
  Clone = data$Clone, 
  xist_levels = logcounts(data)["Xist", ], 
  average_escape = colMeans(ratios, na.rm = T),
  with_xist = unlist(diffd_with_xist[data$Clone])
) %>%
  ggplot(aes(x = xist_levels, y = average_escape, col = with_xist)) + geom_point(size = 5) + 
  geom_smooth(method = "lm", se = F) + ggrepel::geom_text_repel(aes(label = Clone)) + 
  theme_paper() + ylim(0, 0.25)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### check karyotypes!!!
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

data <- readRDS("./ProcessedData/merged_dataset.rds")

data <- data[rowSums(counts_active(data)) + rowSums(counts_inactive(data)) > 100, ]

data.frame(counts_active(data) / (counts_active(data) + counts_inactive(data))) %>%
  add_column("chromosome" = as.character(seqnames(rowRanges(data)))) %>%
  rownames_to_column("gene") %>% pivot_longer(-c(chromosome, gene)) %>%
  add_column("clone" = colData(data)[.$name, ]$Clone) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::filter(chromosome != "Y") %>%
  ggplot(aes(x = factor(chromosome, levels = c(as.character(1:19), c("X", "Y", "MT"))), 
             y = value, col = clone)) + stat_summary(stat = "mean") + facet_wrap(~name) +
  geom_hline(yintercept = 0.5) + theme_paper()

# ranges_data <- rowRanges(data)
# ranges_data@elementMetadata <- DataFrame(data.frame(counts_active(data) / (counts_active(data) + counts_inactive(data))))
# ranges(ranges_data) <- IRanges(start = start(ranges_data), end = start(ranges_data), names = names(ranges_data))
# 
# tiled_genome <- keepStandardChromosomes(seqinfo(EnsDb.Mmusculus.v79))
# tiled_genome <- tiled_genome[c(as.character(1:19), "X"), ]
# tiled_genome_here <- tileGenome(tiled_genome, tilewidth = 1e7, cut.last.tile.in.chrom = T)
# 
# hits <- findOverlaps(ranges_data, tiled_genome_here)
# 
# overlaps <- findOverlaps(ranges_data, tiled_genome_here)
# ranges_data <- ranges_data[from(overlaps)]
# averagedSignal <- aggregate(ranges_data@elementMetadata, to(overlaps), function(x){print(x)})

# 
# ggsave("~/Desktop/karyotypes.pdf")


#### 

# df <- data.frame(
#   x = rnorm(100), 
#   y = rnorm(100)
# ) 
# 
# df %>%
#   ggplot(aes(x = x, y = y)) + geom_point() + 
#   theme( 
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank()
#   ) + 
#   theme(
#     axis.title.x=element_text(hjust = 0.02),
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank(), 
#     axis.title.y=element_text(hjust = 0.02),
#     axis.text.y=element_blank(),
#     axis.ticks.y=element_blank()
#   ) + 
#   theme( 
#     axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm")), ), 
#     axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"))), 
#   ) + 
#   xlab("UMAP1") + 
#   ylab("UMAP2") + 
#   guides(x = ggh4x::guide_axis_truncated(trunc_upper = unit(0.1, "npc")), 
#          y = ggh4x::guide_axis_truncated(trunc_upper = unit(0.1, "npc")), 
#   )
  





