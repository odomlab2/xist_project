# fit rates to 3d - 21d timecourse: 

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

# Now read our pre-processed dataset
data <- readRDS("./ProcessedData/merged_dataset.rds")
data <- computeSumFactors(data)

# For the first analysis, exclude the Xist-knockdown experiments
data <- data[,!data$Guide %in% c("g13", "P3")]
data <- data[,data$Experiment != "May2022"]

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
# Now, we look at the effect of Xist-treatment on escapee expression
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### Now we look at Xist overexpression and how that silences escapees
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO", 
                                         "Aux_0_Dox_14_WO_0_WOAuxNO", "Aux_0_Dox_21_WO_0_WOAuxNO")]

#data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO")]

total_counts <- rowSums(counts_active(data_here) + counts_inactive(data_here))
# average_escape <- rowMeans(counts_active(data_here) + counts_inactive(data_here))

data_here <- data_here[,data_here$Clone == "E6"]

ratios <- counts_inactive(data_here) / (counts_active(data_here) + counts_inactive(data_here))
escapees <- apply(ratios, 2, function(x){x > 0.05 & x < 0.7})
escapees[is.na(escapees)] <- FALSE
general_escapees <- data_here[names(which(apply(escapees[,data_here$Condition == "Aux_0_Dox_0_WO_0_WOAuxNO"], 1, function(x){all(x)}))), ]

escapees <- apply(ratios, 2, function(x){x > 0.1})
escapees_per_clone <- apply(escapees, 2, function(x){rownames(escapees)[x]})
escapees_per_clone <- lapply(escapees_per_clone, function(x){x[x != "Xist"]})
names(escapees_per_clone) <- c("CL30 - December2021", "CL31 - December2021", "E6 - October2021", "E6 - March2022", "E6 - September2022", 
                               "CL30 - October2022", "CL31 - October2022")
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
    labs(fill = "Escape Category") + xlab("")

plot_gene <- function(gene){
  data.frame(
    active = as.numeric(counts_active(data_here)[gene, ]), 
    inactive = as.numeric(counts_inactive(data_here)[gene, ]), 
    timepoint = as.numeric(data_here$ndDox), 
    experiment = data_here$Experiment
  ) %>% mutate(rates = inactive / (inactive + active)) %>%
    ggplot(aes(x = timepoint, y = rates, col = experiment)) + geom_point() + theme_paper()
}

genes_out <- c()
#genes_out <- c("RP23-135F22.6", "Gpkow")
#genes_out <- c(genes_out)

general_escapees <- general_escapees[!rownames(general_escapees) %in% genes_out, ]

fit_curve_for_gene <- function(gene){
  print(gene)
  df_here <- data.frame(
    active = as.numeric(counts_active(data_here)[gene, ]) + 1, 
    inactive = as.numeric(counts_inactive(data_here)[gene, ]) + 1, 
    timepoint = as.numeric(data_here$ndDox), 
    experiment = data_here$Experiment
  ) %>% mutate(total = active + inactive) %>%
    mutate(rate = inactive / (active + inactive))
  #bin.glm <- glm(cbind(inactive, total) ~ timepoint, family = binomial(link = log), data = df_here)
  rates.nls <- lm(log(rate) ~ 1 + timepoint, data = df_here)
  return(c(coef(rates.nls), summary(rates.nls)$adj.r.squared))
}

all_coefs <- data.frame(do.call("rbind", lapply(rownames(general_escapees), fit_curve_for_gene)))
all_coefs$gene <- rownames(general_escapees)
all_coefs$halflifes <- log(1/2) / ( all_coefs$timepoint )
all_coefs$total <- total_counts[all_coefs$gene]

all_coefs %>% ggplot() + geom_histogram(aes(x = halflifes), bins = 100) 
all_coefs %>% ggplot(aes(x = reorder(gene, halflifes), y = halflifes)) + geom_point() + ylim(c(-1, 20)) + coord_flip()

all_coefs$category <- rowData(data_here[all_coefs$gene, ])$EscapeAnnotation

all_coefs %>%
  ggplot(aes(x = halflifes, y = V3, col = category)) + geom_point() + xlab("Halflifes") + ylab("Regression R^2")

# all_coefs %>% ggplot(aes(x = category, y = halflifes)) + geom_boxplot() + geom_jitter() + 
#   ylim(c(-1, 20))
# 
# all_coefs %>% 
#   dplyr::arrange(V3) %>%
#   head(n = 50)
# 
# all_coefs %>% 
#   dplyr::arrange(halflifes) %>%
#   head(n = 50)
# 
# all_coefs %>% ggplot(aes(x = total, y = halflifes)) + geom_point() + scale_x_log10() + ylim(c(-1, 20))

all_coefs_filt <- all_coefs %>% dplyr::filter(V3 > 0.3)

all_coefs %>% 
  dplyr::arrange(halflifes) %>%
  head(n = 50)

all_coefs_filt %>%
  ggplot(aes(x = category, y = halflifes, fill = category)) + geom_boxplot(width = 0.1, outlier.color = NA) + ggbeeswarm::geom_quasirandom()

all_coefs_filt %>% ggplot(aes(x = category, y = halflifes)) + geom_boxplot() + geom_jitter()

all_coefs_filt %>% 
  dplyr::arrange(halflifes) %>%
  head(n = 30)

gene = "Rps4x"

df_here <- data.frame(
  active = as.numeric(counts_active(data_here)[gene, ]) + 1, 
  inactive = as.numeric(counts_inactive(data_here)[gene, ]) + 1, 
  timepoint = as.numeric(data_here$ndDox), 
  experiment = data_here$Experiment
) %>% mutate(total = active + inactive) %>% mutate(rate = inactive / (active + inactive))
#bin.glm <- glm(cbind(inactive, total) ~ timepoint, family = binomial(link = log), data = df_here)
rates.nls <- lm(log(rate) ~ timepoint, data = df_here)
#bin.glm <- lm(inactive / total ~ timepoint, data = df_here)

neg_exp = function(x){
  exp(coef(rates.nls)[[1]]) * exp( (coef(rates.nls)[[2]]) * x)
}

data.frame(
  active = as.numeric(counts_active(data_here)[gene, ]) + 1, 
  inactive = as.numeric(counts_inactive(data_here)[gene, ]) + 1, 
  timepoint = as.numeric(data_here$ndDox), 
  experiment = data_here$Experiment
) %>% mutate(rates = inactive / (inactive + active)) %>%
  ggplot(aes(x = timepoint, y = rates, col = experiment)) + geom_point() + theme_paper() + 
  geom_function(fun = neg_exp, color = "red") + ylab("allelic ratio") + 
  ylim(c(0, 1)) + ggtitle(gene) + 
  annotate("text", label = paste0("R^2: ", round(all_coefs_filt[all_coefs_filt$gene == gene, ]$V3, digits = 4)), x = 2, y = 1, size = 8) + 
  annotate("text", label = paste0("t1/2: ", round(all_coefs_filt[all_coefs_filt$gene == gene, ]$halflifes, digits = 4)), x = 2, y = 0.9, size = 8)

# data.frame(
#   active = as.numeric(counts_active(data_here)[gene, ]) + 1,
#   inactive = as.numeric(counts_inactive(data_here)[gene, ]) + 1,
#   timepoint = as.numeric(data_here$ndDox),
#   experiment = data_here$Experiment
# ) %>% mutate(rates = inactive / (inactive + active)) %>%
#   ggplot(aes(x = timepoint, y = rates, col = experiment)) + geom_point() + theme_paper() +
#   geom_function(fun = neg_exp, color = "red") + scale_y_log10() + ylab("allelic ratio")




### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Now, we look at the effect of Xist-treatment on escapee expression
### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# some things look a bit weird - try also with only 0-7d timepoints

### Now we look at Xist overexpression and how that silences escapees
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO", 
                                         "Aux_0_Dox_14_WO_0_WOAuxNO", "Aux_0_Dox_21_WO_0_WOAuxNO")]

data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO")]

total_counts <- rowSums(counts_active(data_here) + counts_inactive(data_here))
# average_escape <- rowMeans(counts_active(data_here) + counts_inactive(data_here))

data_here <- data_here[,data_here$Clone == "E6"]

ratios <- counts_inactive(data_here) / (counts_active(data_here) + counts_inactive(data_here))
escapees <- apply(ratios, 2, function(x){x > 0.05 & x < 0.7})
escapees[is.na(escapees)] <- FALSE
general_escapees <- data_here[names(which(apply(escapees[,data_here$Condition == "Aux_0_Dox_0_WO_0_WOAuxNO"], 1, function(x){all(x)}))), ]

escapees <- apply(ratios, 2, function(x){x > 0.1})
escapees_per_clone <- apply(escapees, 2, function(x){rownames(escapees)[x]})
escapees_per_clone <- lapply(escapees_per_clone, function(x){x[x != "Xist"]})
names(escapees_per_clone) <- c("CL30 - December2021", "CL31 - December2021", "E6 - October2021", "E6 - March2022", "E6 - September2022", 
                               "CL30 - October2022", "CL31 - October2022")
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
  labs(fill = "Escape Category") + xlab("")

plot_gene <- function(gene){
  data.frame(
    active = as.numeric(counts_active(data_here)[gene, ]), 
    inactive = as.numeric(counts_inactive(data_here)[gene, ]), 
    timepoint = as.numeric(data_here$ndDox), 
    experiment = data_here$Experiment
  ) %>% mutate(rates = inactive / (inactive + active)) %>%
    ggplot(aes(x = timepoint, y = rates, col = experiment)) + geom_point() + theme_paper()
}

genes_out <- c()
#genes_out <- c("RP23-135F22.6", "Gpkow")
#genes_out <- c(genes_out)

general_escapees <- general_escapees[!rownames(general_escapees) %in% genes_out, ]

fit_curve_for_gene <- function(gene){
  print(gene)
  df_here <- data.frame(
    active = as.numeric(counts_active(data_here)[gene, ]) + 1, 
    inactive = as.numeric(counts_inactive(data_here)[gene, ]) + 1, 
    timepoint = as.numeric(data_here$ndDox), 
    experiment = data_here$Experiment
  ) %>% mutate(total = active + inactive) %>%
    mutate(rate = inactive / (active + inactive))
  #bin.glm <- glm(cbind(inactive, total) ~ timepoint, family = binomial(link = log), data = df_here)
  rates.nls <- lm(log(rate) ~ 1 + timepoint, data = df_here)
  return(c(coef(rates.nls), summary(rates.nls)$adj.r.squared))
}

all_coefs <- data.frame(do.call("rbind", lapply(rownames(general_escapees), fit_curve_for_gene)))
all_coefs$gene <- rownames(general_escapees)
all_coefs$halflifes <- log(1/2) / ( all_coefs$timepoint )
all_coefs$total <- total_counts[all_coefs$gene]

all_coefs %>% ggplot() + geom_histogram(aes(x = halflifes), bins = 100) 
all_coefs %>% ggplot(aes(x = reorder(gene, halflifes), y = halflifes)) + geom_point() + ylim(c(-1, 20)) + coord_flip()

all_coefs$category <- rowData(data_here[all_coefs$gene, ])$EscapeAnnotation

all_coefs %>%
  dplyr::filter(!is.na(category)) %>%
  ggplot(aes(x = halflifes, y = V3, col = category)) + geom_point() + xlab("Halflifes") + ylab("Regression R^2")

# all_coefs %>% ggplot(aes(x = category, y = halflifes)) + geom_boxplot() + geom_jitter() + 
#   ylim(c(-1, 20))
# 
# all_coefs %>% 
#   dplyr::arrange(V3) %>%
#   head(n = 50)
# 
# all_coefs %>% 
#   dplyr::arrange(halflifes) %>%
#   head(n = 50)
# 
# all_coefs %>% ggplot(aes(x = total, y = halflifes)) + geom_point() + scale_x_log10() + ylim(c(-1, 20))

all_coefs_filt <- all_coefs %>% dplyr::filter(V3 > 0.3)

all_coefs %>% 
  dplyr::arrange(halflifes) %>%
  head(n = 50)

all_coefs_filt %>%
  ggplot(aes(x = category, y = halflifes, fill = category)) + geom_boxplot(width = 0.1, outlier.color = NA) + ggbeeswarm::geom_quasirandom()

all_coefs_filt %>% ggplot(aes(x = category, y = halflifes)) + geom_boxplot() + geom_jitter()

all_coefs_filt %>% 
  dplyr::arrange(halflifes) %>%
  head(n = 30)

gene = "Eif2s3x"

df_here <- data.frame(
  active = as.numeric(counts_active(data_here)[gene, ]) + 1, 
  inactive = as.numeric(counts_inactive(data_here)[gene, ]) + 1, 
  timepoint = as.numeric(data_here$ndDox), 
  experiment = data_here$Experiment
) %>% mutate(total = active + inactive) %>% mutate(rate = inactive / (active + inactive))
#bin.glm <- glm(cbind(inactive, total) ~ timepoint, family = binomial(link = log), data = df_here)
rates.nls <- lm(log(rate) ~ timepoint, data = df_here)
#bin.glm <- lm(inactive / total ~ timepoint, data = df_here)

neg_exp = function(x){
  exp(coef(rates.nls)[[1]]) * exp( (coef(rates.nls)[[2]]) * x)
}

data.frame(
  active = as.numeric(counts_active(data_here)[gene, ]) + 1, 
  inactive = as.numeric(counts_inactive(data_here)[gene, ]) + 1, 
  timepoint = as.numeric(data_here$ndDox), 
  experiment = data_here$Experiment
) %>% mutate(rates = inactive / (inactive + active)) %>%
  ggplot(aes(x = timepoint, y = rates, col = experiment)) + geom_point() + theme_paper() + 
  geom_function(fun = neg_exp, color = "red") + ylab("allelic ratio") + 
  ylim(c(0, 1)) + ggtitle(gene) + 
  annotate("text", label = paste0("R^2: ", round(all_coefs_filt[all_coefs_filt$gene == gene, ]$V3, digits = 4)), x = 2, y = 1, size = 8) + 
  annotate("text", label = paste0("t1/2: ", round(all_coefs_filt[all_coefs_filt$gene == gene, ]$halflifes, digits = 4)), x = 2, y = 0.9, size = 8)


### ### ### ### ### ### ### ### ### ### ### 
### Try fitting nls model with and without offset
### ### ### ### ### ### ### ### ### ### ### 


### Now we look at Xist overexpression and how that silences escapees
data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO", 
                                         "Aux_0_Dox_14_WO_0_WOAuxNO", "Aux_0_Dox_21_WO_0_WOAuxNO")]

#data_here <- data[,data$Condition %in% c("Aux_0_Dox_0_WO_0_WOAuxNO", "Aux_0_Dox_3_WO_0_WOAuxNO", "Aux_0_Dox_7_WO_0_WOAuxNO")]

total_counts <- rowSums(counts_active(data_here) + counts_inactive(data_here))
# average_escape <- rowMeans(counts_active(data_here) + counts_inactive(data_here))

data_here <- data_here[,data_here$Clone == "E6"]

ratios <- counts_inactive(data_here) / (counts_active(data_here) + counts_inactive(data_here))
escapees <- apply(ratios, 2, function(x){x > 0.05 & x < 0.7})
escapees[is.na(escapees)] <- FALSE
general_escapees <- data_here[names(which(apply(escapees[,data_here$Condition == "Aux_0_Dox_0_WO_0_WOAuxNO"], 1, function(x){all(x)}))), ]

escapees <- apply(ratios, 2, function(x){x > 0.1})
escapees_per_clone <- apply(escapees, 2, function(x){rownames(escapees)[x]})
escapees_per_clone <- lapply(escapees_per_clone, function(x){x[x != "Xist"]})
names(escapees_per_clone) <- c("CL30 - December2021", "CL31 - December2021", "E6 - October2021", "E6 - March2022", "E6 - September2022", 
                               "CL30 - October2022", "CL31 - October2022")

genes_out <- c()
#genes_out <- c("RP23-135F22.6", "Gpkow")
#genes_out <- c(genes_out)

general_escapees <- general_escapees[!rownames(general_escapees) %in% genes_out, ]

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
  #rates.nls <- nls(rate ~ a * exp(- k * timepoint), data = df_here)
  tryCatch({
    rates.nls <- nls(rate ~ a * exp(- k * timepoint), data = df_here, algorithm = "port")
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
  #rates.nls <- nls(rate ~ a * exp(- k * timepoint), data = df_here)
  tryCatch({
    rates.nls.offset <- nls(rate ~ a * exp(- k * timepoint) + b, data = df_here, lower = c(-Inf, -Inf, 0), algorithm = "port")
    return(c(coef(rates.nls.offset), modelr::rsquare(rates.nls.offset, data = df_here), logLik(rates.nls.offset)))
  }, 
  error=function(cond) {
    return(NA)
  })
}

all_coefs <- data.frame(do.call("rbind", lapply(rownames(general_escapees), fit_nls_curve_for_gene)))
all_coefs$gene <- rownames(general_escapees)
all_coefs$halflifes <- log(1/2) / ( - all_coefs$k )
all_coefs$total <- total_counts[all_coefs$gene]
# 
all_coefs_offset <- data.frame(do.call("rbind", lapply(rownames(general_escapees), fit_nls_curve_for_gene_offset)))
all_coefs_offset$gene <- rownames(general_escapees)
all_coefs_offset$halflifes <- log(1/2) / ( - all_coefs_offset$k )
all_coefs_offset$total <- total_counts[all_coefs_offset$gene]

all_coefs_offset %>% ggplot() + geom_histogram(aes(x = halflifes), bins = 100) 
all_coefs_offset %>% ggplot(aes(x = reorder(gene, halflifes), y = halflifes)) + geom_point() + ylim(c(-1, 20)) + coord_flip()

all_coefs_offset$category <- rowData(data_here[all_coefs$gene, ])$EscapeAnnotation

all_coefs %>%
  ggplot(aes(x = halflifes, y = V3)) + geom_point() + xlab("Halflifes") + ylab("Regression R^2")


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
  
merged_df %>% 
  dplyr::filter(!(r2 > 0.2 & r2_off > 0.2)) %>% dim()

merged_df <- merged_df %>% 
  dplyr::filter(r2 > 0.2 & r2_off > 0.2)

merged_df %>%
  ggplot(aes(x = r2)) + geom_histogram()

merged_df %>%
  ggplot(aes(x = r2_off)) + geom_histogram()

merged_df %>%
  ggplot(aes(x = r2, r2_off, col = category)) + geom_point() + xlim(c(0,1 )) + ylim(c(0, 1)) + geom_abline() +
    ggrepel::geom_text_repel(aes(label = gene))

merged_df %>% dplyr::arrange( - b_off)

gene = "Armcx5"

df_here <- data.frame(
  active = as.numeric(counts_active(data_here)[gene, ]) + 1, 
  inactive = as.numeric(counts_inactive(data_here)[gene, ]) + 1, 
  timepoint = as.numeric(data_here$ndDox), 
  experiment = data_here$Experiment
) %>% mutate(total = active + inactive) %>% mutate(rate = inactive / (active + inactive))
#bin.glm <- glm(cbind(inactive, total) ~ timepoint, family = binomial(link = log), data = df_here)
rates.nls <- lm(log(rate) ~ timepoint, data = df_here)
#bin.glm <- lm(inactive / total ~ timepoint, data = df_here)

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
  annotate("text", label = paste0("R^2: ", round(all_coefs_filt[all_coefs_filt$gene == gene, ]$V3, digits = 4)), x = 2, y = 1, size = 8) + 
  annotate("text", label = paste0("t1/2: ", round(all_coefs_filt[all_coefs_filt$gene == gene, ]$halflifes, digits = 4)), x = 2, y = 0.9, size = 8)

# look at BICs between models: 

merged_df %>%
  ggplot(aes(x = BIC_native, BIC_offset, col = category)) + geom_point() + geom_abline() +
    ggrepel::geom_text_repel(aes(label = gene))

merged_df %>%
  ggplot(aes(x = halflifes_off, b_off, col = BIC_native > BIC_offset)) + geom_point() + xlim(c(0, 15)) +  ylim(c(0, 0.3)) + ggrepel::geom_text_repel(aes(label = gene))



