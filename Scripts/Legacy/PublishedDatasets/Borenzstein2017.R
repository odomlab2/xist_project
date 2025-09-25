

### Read data
setwd("~/Desktop/Projects/XChromosome_Antonia/")

library(tidyverse)
library(ggplot2)
library(scran)
library(scater)

data <- readxl::read_xls("~/Desktop/GSE80810_AllelicRatio.xls")
data_x <- data[data$Chromosomes == "chrX", ] %>%
  dplyr::select(-Chromosomes) %>%
  column_to_rownames("Genes")

cell_names <- colnames(data_x)
stage <- unlist(lapply(str_split(colnames(data_x), "_"), function(x){x[[1]]}))
cross <- unlist(lapply(str_split(colnames(data_x), "_"), function(x){x[[2]]}))
cell_number <- unlist(lapply(str_split(colnames(data_x), "_"), function(x){x[[3]]}))

data_x <- t(data.frame(apply(data_x, 1, as.numeric)))
colnames(data_x) <- cell_names

data.frame(
  ase = colMeans(data_x, na.rm = T), 
  stage = factor(stage, levels = c("Oo", "2C", "4C", "8C", "16C", "32C", "64C")), 
  cross = cross
) %>%
  ggplot(aes(x = stage, y = ase, col = cross)) + ggbeeswarm::geom_quasirandom()

data.frame(
  ase = data_x["Kdm6a", ], 
  stage = factor(stage, levels = c("Oo", "2C", "4C", "8C", "16C", "32C", "64C")), 
  cross = cross
) %>%
  ggplot(aes(x = stage, y = ase, fill = cross)) + geom_boxplot()






