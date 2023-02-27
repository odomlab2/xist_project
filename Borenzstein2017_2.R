library(tidyverse)

data <- readxl::read_excel("~/Desktop/Projects/XChromosome_Antonia/Data/Bronzstein2017/GSE80810_AllelicRatio.xls")
data <- data[data$Chromosomes == "chrX", ]

data <- data %>%
  column_to_rownames("Genes") %>%
  select(-Chromosomes)

stage <- unlist(lapply(str_split(colnames(data), "_"), function(x){x[[1]]}))
cross <- unlist(lapply(str_split(colnames(data), "_"), function(x){x[[2]]}))

data.frame(
  stage = factor(stage, levels = c("Oo", "2C", "4C", "8C", "16C", "32C", "64C")), 
  cross = cross, 
  allele_ratio = as.numeric(data["Ddx3x", ])
) %>%
  ggplot(aes(x = stage, y = allele_ratio, col = cross)) + geom_jitter(width = 0.1, height = 0.01)

means <- apply(data[,stage %in% c("32C", "64C")], 1, function(x){mean(as.numeric(x), na.rm = T)})
sum_not_na <- apply(data[,stage %in% c("32C", "64C")], 1, function(x){sum(!is.na(as.numeric(x)))})

data.frame(
  genes = rownames(data),
  means, 
  sum_not_na
) %>%
  ggplot(aes(x = sum_not_na, means)) + geom_point() + 
    ggrepel::geom_text_repel(aes(label = genes))

genes_check <- data.frame(
  genes = rownames(data),
  means, 
  sum_not_na
) %>%
  filter(means < 0.5 & means > 0.1 & sum_not_na > 20)

data_here <- t(apply(data[,stage %in% c("32C", "64C")], 1, as.numeric))
colnames(data_here) <- colnames(data[,stage %in% c("32C", "64C")])

cor_matrix <- cor(data_here[rownames(genes_check), ], use = "pairwise.complete.obs", method = "spearman")
clustering_cells <- hclust(as.dist(cor_matrix))
plot(clustering_cells)

ComplexHeatmap::Heatmap(data_here[rownames(genes_check), clustering_cells$order], cluster_rows = F, cluster_columns = F)
ComplexHeatmap::Heatmap(data_here[rownames(genes_check), clustering_cells$order])

