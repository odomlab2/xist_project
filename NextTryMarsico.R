## Ma-ma-ma-maaaaarsico

# Re-run script: 
# Load coefficients + data
# look at coefficient + data distributions
# run RF model (also try regularized mv regression against non-binarized target)
# plot CV-AUC
# For each feature plot for each CV iteration the MDA distribution
# Look at feature distributions

theme_paper <- function(){
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        text = element_text(size = 20))
}


clean_names <- function(.data, unique = FALSE) {
  # thanks, https://www.r-bloggers.com/2019/07/clean-consistent-column-names/
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


# Load the fitted coefficients
setwd("~/Desktop/Projects/XChromosome_Antonia/")

coefficients <- readRDS("./ProcessedData/LM_coefficients_E6.rds")
all_coefs <- coefficients

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
constitutive_escapees <- unique(escapee_annotation[escapee_annotation$`final status` == "E", ]$symbol)
facultative_escapees <- unique(escapee_annotation[escapee_annotation$`final status` == "V", ]$symbol)

# load their feature-set
load("~/Desktop/Projects/XChromosome_Antonia/Data/Supplemental_Code/data/feature_matrix/promoter_matrix_reannotated_normRAdjusted_all_chrX_genes.RData")

# overlap between our gene list and theirs
both_genes <- intersect(rownames(data_set), rownames(all_coefs))
length(both_genes)

data_set_fit <- data_set[both_genes, ]
data_set_fit$is_constitutive <- 0
data_set_fit$is_facultative <- 0
data_set_fit[rownames(data_set_fit) %in% constitutive_escapees, ]$is_constitutive <- 1
data_set_fit[rownames(data_set_fit) %in% facultative_escapees, ]$is_facultative <- 1
colnames(data_set_fit) <- clean_names(colnames(data_set_fit))

# Check distribution of features
# ...

#data_set_fit$distance_to_LADs <- log10(data_set_fit$distance_to_LADs + 1)
#data_set_fit$distance_to_Xist <- log10(data_set_fit$distance_to_Xist + 1)
#data_set_fit$distance_to_TAD_border <- log10(data_set_fit$distance_to_TAD_border + 1)
#data_set_fit$distance_to_LINE <- log10(data_set_fit$distance_to_LINE + 1)
data_set_fit_scaled <- data.frame(apply(data_set_fit, 2, function(x){x <- as.numeric(x); return((x - mean(x)) / sd(x))}))
coefs_fit_here <- all_coefs[both_genes, ]
rownames(data_set_fit_scaled) <- rownames(coefs_fit_here)

library(randomForest)
require(pROC)

rf_cv_auc <- function(data_fit, folds = 5){     
  n_chunks <- floor(nrow(data_fit) / folds)
  data_fit <- data_fit[sample(nrow(data_fit)), ]
  df_folds <- split(data_fit, (seq(nrow(data_fit))-1) %/% n_chunks)
  
  importances <- lapply(1:(length(df_folds) - 1), function(i){
    train_data <- do.call('rbind', df_folds[-i])
    test_data <- df_folds[[i]]
    rf <- randomForest(
      objective ~ .,
      data = train_data,
      importance = T, 
    )
    prediction <- stats::predict(rf, test_data, type = "vote")
    rf.roc <- roc(test_data$objective, prediction[,2])
    
    return(list(auc(rf.roc), cbind(importance(rf))))
  })
  
  aucs <- unlist(lapply(importances, function(x){x[[1]]}))
  importances <- do.call("rbind", lapply(importances, function(x) x[[2]] %>% data.frame() %>% rownames_to_column("FeatureName")))
  
  return(list(aucs, importances))
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

coefs_fit_here$ratio_3d <- ((coefs_fit_here$ndDox3 + 0.1) / (coefs_fit_here$washout_dayYES_3 + 0.1))
coefs_fit_here$ratio_7d <- ((coefs_fit_here$ndDox7 + 0.1) / (coefs_fit_here$washout_dayYES_7 + 0.1))

coefs_fit_here$ratio_3d <- - ((coefs_fit_here$ndDox3) + (coefs_fit_here$washout_dayYES_3))
coefs_fit_here$ratio_7d <- - ((coefs_fit_here$ndDox7) + (coefs_fit_here$washout_dayYES_7))

# 
coefs_fit_here %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(x = ndDox7, y = washout_dayYES_7, size = ratio_7d)) + geom_point() + geom_abline(slope = -1) + 
  theme_paper() + xlab("Silencing effect (7d)") + ylab("Washout effect (7d)") + labs(size = "Reversibility \n(Silencing - Washout)")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures/Fig_Marsico_prelim/visualization_reversibility.pdf")

data_set_fit_scaled$objective <- coefs_fit_here$ndDox7
train_data <- data_set_fit_scaled
train_data$objective <- factor(train_data$objective < mean(train_data$objective))
rownames(train_data) <- rownames(coefs_fit_here)

rf <- randomForest(
  objective ~ .,
  data = train_data,
  importance = T, 
)

# objective distribution: 
data.frame(X = data_set_fit_scaled$objective) %>%
  ggplot(aes(x = - X)) + geom_histogram(aes(fill = X > mean(data_set_fit_scaled$objective) )) + scale_fill_manual(values = c("black", "grey")) + 
  xlab("Silencing strength") + theme_paper() + labs(fill = "high / low")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures/Fig_Marsico_prelim/stratification_silencing.pdf")


set.seed(1234)
aucs_here <- rf_cv_auc(train_data)
aucs_here[[1]] %>%
  data.frame() %>%
  ggplot(aes(x = "AUC", y = .)) + geom_boxplot(fill = "darkred", col = "black") + ylim(0.5, 1) + theme_classic() + coord_flip() + 
  geom_jitter(size = 10, width = 0.1) + theme(text = element_text(size = 20)) + xlab("") + theme_paper()
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures/Fig_Marsico_prelim/aucs_silencing.pdf")

importance <- data.frame(importance(rf))
importances <- aucs_here[[2]]

importances %>%
  data.frame() %>%
  group_by(FeatureName) %>%
  mutate(AverageMDA = mean(MeanDecreaseAccuracy)) %>%
  arrange(-AverageMDA) %>%
  head(n = 5 * 50)

importances %>%
  data.frame() %>%
  group_by(FeatureName) %>%
  mutate(AverageMDA = mean(MeanDecreaseAccuracy)) %>%
  arrange(-AverageMDA) %>%
  #head(n = 5 * 50) %>%
  ggplot(aes(x = reorder(FeatureName, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) + stat_summary() + coord_flip()

importances %>%
  data.frame() %>%
  group_by(FeatureName) %>%
  mutate(AverageMDA = mean(MeanDecreaseAccuracy)) %>%
  arrange(-AverageMDA) %>%
  head(n = 5 * 30) %>%
  ggplot(aes(x = reorder(FeatureName, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) + stat_summary() + coord_flip() + theme_paper() + xlab("")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures/Fig_Marsico_prelim/mdas_silencing.pdf")

data_set_fit_plot <- data.frame(apply(data_set_fit, 2, as.numeric))
data_set_fit_plot$objective <- train_data$objective

feature_check <- "distance_to_la_ds"
data_set_fit_plot %>%
  data.frame() %>%
  pivot_longer(-objective) %>%
  add_column(MDA = importance[.$name, ]$MeanDecreaseAccuracy) %>%
  dplyr::filter(name %in% c(feature_check)) %>%
  ggplot(aes(x = reorder(name, MDA), y = value, fill = objective)) + geom_boxplot() + ggpubr::stat_compare_means() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) + theme_paper() + xlab("")
  ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures/Fig_Marsico_prelim/silencing_lads.pdf")

train_data_here <- train_data
train_data_here$objective <- coefs_fit_here$ndDox7
train_data_here$feature_plot <- train_data_here[,feature_check]
train_data_here %>%
  ggplot(aes(x = - objective, y = feature_plot)) + geom_point() + geom_smooth(method = "lm") + ylab(feature_check) + 
  theme_paper() + xlab(" - Dox effect (high means strong effect)")

train_data_here <- data_set_fit
train_data_here$objective <- coefs_fit_here$ndDox7
train_data_here$feature_plot <- train_data_here[,feature_check]
train_data_here %>%
  ggplot(aes(x = - objective, y = feature_plot)) + geom_point() + geom_smooth(method = "lm") + ylab(feature_check) + 
  theme_paper() + xlab(" - Dox effect (high means strong effect)")

feature_check <- "hdac3_gse116480_minus_500_500"
train_data_here <- data_set_fit
train_data_here$objective <- coefs_fit_here$ndDox7
train_data_here$feature_plot <- train_data_here[,feature_check]
train_data_here %>%
  ggplot(aes(x = - objective, y = feature_plot)) + geom_point() + geom_smooth(method = "lm") + ylab(feature_check) + 
  theme_paper() + xlab(" - Dox effect (high means strong effect)")

feature_check <- "is_constitutive"
train_data_here <- data_set_fit
train_data_here$objective <- coefs_fit_here$ndDox7
train_data_here$feature_plot <- train_data_here[,feature_check]
train_data_here %>%
  ggplot(aes(x = - objective, y = feature_plot)) + geom_point() + geom_smooth(method = "lm") + ylab(feature_check) + 
  theme_paper() + xlab(" - Dox effect (high means strong effect)")

#### predict based on reversibility
data_set_fit_scaled$objective <- coefs_fit_here$ratio_7d
train_data <- data_set_fit_scaled
train_data$objective <- factor(train_data$objective < mean(train_data$objective))
rownames(train_data) <- rownames(coefs_fit_here)

rf <- randomForest(
  objective ~ .,
  data = train_data,
  importance = T, 
)

# objective distribution: 
data.frame(X = data_set_fit_scaled$objective) %>%
  ggplot(aes(x = X)) + geom_histogram(aes(fill = X > mean(data_set_fit_scaled$objective) )) + scale_fill_manual(values = c("black", "grey")) + 
  xlab("Irreversibility") + theme_paper() + labs("high / low")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures/Fig_Marsico_prelim/stratification_reversibility.pdf")

set.seed(1234)
aucs_here <- rf_cv_auc(train_data)
aucs_here[[1]] %>%
  data.frame() %>%
  ggplot(aes(x = "AUC", y = .)) + geom_boxplot(fill = "darkred", col = "black") + ylim(0.5, 1) + theme_classic() + coord_flip() + 
  geom_jitter(size = 10, width = 0.1) + theme(text = element_text(size = 20)) + theme_paper() + xlab("")
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures/Fig_Marsico_prelim/aucs_irreversibility.pdf")

importance <- data.frame(importance(rf))
importances <- aucs_here[[2]]

importances %>%
  data.frame() %>%
  group_by(FeatureName) %>%
  mutate(AverageMDA = mean(MeanDecreaseAccuracy)) %>%
  arrange(-AverageMDA) %>%
  head(n = 5 * 30)

importances %>%
  data.frame() %>%
  group_by(FeatureName) %>%
  mutate(AverageMDA = mean(MeanDecreaseAccuracy)) %>%
  arrange(-AverageMDA) %>%
  ggplot(aes(x = reorder(FeatureName, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) + stat_summary() + coord_flip()

importances %>%
  data.frame() %>%
  group_by(FeatureName) %>%
  mutate(AverageMDA = mean(MeanDecreaseAccuracy)) %>%
  arrange(-AverageMDA) %>%
  head(n = 5 * 30) %>%
  ggplot(aes(x = reorder(FeatureName, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) + stat_summary() + coord_flip() + theme_paper() 
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures/Fig_Marsico_prelim/mdas_irreversibility.pdf")


data_set_fit_plot <- data.frame(apply(data_set_fit, 2, as.numeric))
data_set_fit_plot$objective <- train_data$objective

feature_check <- "taf1_gse30959_minus_500_500"
data_set_fit_plot %>%
  data.frame() %>%
  pivot_longer(-objective) %>%
  add_column(MDA = importance[.$name, ]$MeanDecreaseAccuracy) %>%
  dplyr::filter(name %in% c(feature_check)) %>%
  ggplot(aes(x = reorder(name, MDA), y = value, fill = objective)) + geom_boxplot() + ggpubr::stat_compare_means() + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) + theme_paper()
ggsave("~/Desktop/Projects/XChromosome_Antonia/Figures/Fig_Marsico_prelim/taf1_irreversibility.pdf")

feature_check <- "h3k4me1_gse47949_minus_1000_1000"
train_data_here <- train_data
train_data_here$objective <- coefs_fit_here$ratio_7d
train_data_here$feature_plot <- train_data_here[,feature_check]
train_data_here %>%
  ggplot(aes(x = objective, y = feature_plot)) + geom_point() + geom_smooth(method = "lm") + ylab(feature_check) + 
  theme_paper() + xlab("Difference Xist / Washout (high means less reversible)")

feature_check <- "hdac1_gse27841_minus_250_250"
train_data_here <- train_data
train_data_here$objective <- coefs_fit_here$ratio_7d
train_data_here$feature_plot <- train_data_here[,feature_check]
train_data_here %>%
  ggplot(aes(x = objective, y = feature_plot)) + geom_point() + geom_smooth(method = "lm") + ylab(feature_check) + 
  theme_paper() + xlab("Difference Xist / Washout (high means less reversible)")

feature_check <- "cpg_content"
train_data_here <- data_set_fit
train_data_here$objective <- coefs_fit_here$ratio_7d
train_data_here$feature_plot <- train_data_here[,feature_check]
train_data_here %>%
  ggplot(aes(x = objective, y = feature_plot)) + geom_point() + geom_smooth(method = "lm") + ylab(feature_check) + 
  theme_paper() + xlab("Difference Xist / Washout (high means less reversible)")

feature_check <- "is_constitutive"
train_data_here <- data_set_fit
train_data_here$objective <- coefs_fit_here$ratio_7d
train_data_here$feature_plot <- train_data_here[,feature_check]
train_data_here %>%
  ggplot(aes(x = objective, y = feature_plot)) + geom_point() + geom_smooth(method = "lm") + ylab(feature_check) + 
  theme_paper() + xlab("Difference Xist / Washout (high means less reversible)")



