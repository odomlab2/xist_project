### Marsico analysis

### given a vector of coefficients, check if they correlate with any genomic / epigenomic features 
### from https://genome.cshlp.org/content/29/7/1087.full.pdf+html

# their escapee-annotation
load("~/Desktop/Projects/XChromosome_Antonia/Data/Supplemental_Code/data/annotation_files/escapees/escapees.RData")

# their features # called data-set
load("~/Desktop/Projects/XChromosome_Antonia/Data/Supplemental_Code/data/feature_matrix/promoter_matrix_reannotated_normRAdjusted_all_chrX_genes.RData")

# overlap between our gene list and theirs
both_genes <- intersect(rownames(data_set), rownames(all_coefs))

data_set_fit <- data_set[both_genes, ]
data_set_fit <- data.frame(apply(data_set_fit, 2, as.numeric)) %>%
  `rownames<-`(both_genes)

# look at entire dataset -- there's a lot of binary / non-quantitative features
# data_set_fit %>%
#   pivot_longer(-c()) %>%
#   ggplot(aes(x = value)) + facet_wrap(~name, scale = "free") + geom_histogram()

data_set_fit_scaled <- data.frame(apply(data_set_fit, 2, function(x){x <- as.numeric(x); return(as.numeric(x - mean(x)) / sd(x))}))
coefs_fit_here <- all_coefs[both_genes, ]

# data_set_fit_scaled %>%
#   pivot_longer(-c()) %>%
#   ggplot(aes(x = value)) + facet_wrap(~name, scale = "free") + geom_histogram()

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

coefs_fit_here %>% 
  ggplot(aes(x = ndDox7, y = washoutYES)) + geom_point(aes(size = abs(ndDox7 / washoutYES))) + geom_abline(slope = -1)


## check for one vector

data_fit <- data_set_fit_scaled

coef_vector <- coefs_fit_here$ndDox7
correlations <- apply(data_fit, 2, function(x){cor(coef_vector, as.numeric(x), method = "spearman")})

names_save <- colnames(data_fit)
names(names_save) <- paste0("Feature", 1:ncol(data_fit))
colnames(data_fit) <- paste0("Feature", 1:ncol(data_fit))

data_fit <- cbind(data_fit, "objective" = coef_vector)
boruta_rf <- Boruta::Boruta(objective ~ . , data_fit)

important_features <- Boruta::getSelectedAttributes(boruta_rf)
names_save[important_features]

rf <- randomForest(
  objective ~ .,
  data=data_fit,
  importance = T
)
rf

importance_df <- data.frame(rf$importance)
importance_df$cor <- correlations
rownames(importance_df) <- names_save

head(importance_df[order(importance_df$X.IncMSE, decreasing = T), ])

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
  feature = data_set_fit[,"distance_to_LADs"]
) %>%
  ggplot(aes(x = feature, y = ratio)) + geom_point() + geom_smooth(method = "lm") + 
  ggrepel::geom_text_repel(aes(label = gene))

rownames(data_set_fit_scaled) <- rownames(data_set_fit)
pheatmap::pheatmap(data_set_fit_scaled[,names_save[important_features]], fontsize = 4)

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