library(tidyverse)
library(SingleCellExperiment)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../../")

source("./Scripts/General/reuse_functions.R")

counts_active <- function(sce){
  assays(sce)["counts_active"][[1]]
}

counts_inactive <- function(sce){
  assays(sce)["counts_inactive"][[1]]
}

#### generate invivo dataset 1

all.files <- list.files("../../invivo_data/processed_data/quant/", pattern = ".counts$", full.names = T)
all.data <- lapply(all.files, function(x){
  test <- read_delim(x, delim = "\t", skip = 1)
  #if(test[1, 1] == "Geneid"){
  #  test = test[-1, ]
  #}
  test[,c(1, 7)]
})

data <- do.call("cbind", lapply(all.data, function(x){x[,2]}))
rownames(data) <- all.data[[1]][,1]$Geneid

# add chromosome information to genes
library("EnsDb.Mmusculus.v79")
gene_info <- data.frame(genes(EnsDb.Mmusculus.v79))
rownames(gene_info) <- gene_info$gene_id
gene_info <- gene_info[!duplicated(gene_info$symbol), ]

# we're dropping duplicate ensembl ids here, which might not be ideal, 
# but there are no affected expressed genes on the X
genes_keep <- intersect(rownames(data), gene_info$gene_id)
gene_info <- gene_info[genes_keep, ]

gene_info <- makeGRangesFromDataFrame(gene_info, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = T)

# read in metadata
metadata <- data.frame(
  row.names = paste0("Sample", c("1", "8", "10", "11", "12")), 
  Sample = factor(paste0("Sample", c("1", "8", "10", "11", "12"))), 
  ShortSampleName = factor(paste0("Sample", c("1", "8", "10", "11", "12"))), 
  Sex = c(rep("female", 4), "male"), 
  RttaCounts = c(0, 0, 0, 0, 0), 
  RttaGenotype = c(rep("Rtta-", 4), "Rtta+")
)

data.all <- data[genes_keep, grepl("Gall", colnames(data))]
data.b6 <- data[genes_keep,grepl("C57BL-6J", colnames(data))]
data.jf1 <- data[genes_keep,grepl("JF1_MSJ", colnames(data))]

colnames(data.all) <- gsub("_Gall_sorted.bam", "", colnames(data.all))
colnames(data.b6) <- gsub("_C57BL-6J_sorted.bam", "", colnames(data.b6))
colnames(data.jf1) <- gsub("_JF1_MSJ_sorted.bam", "", colnames(data.jf1))

rownames(data.all) <- gene_info[genes_keep, ]$symbol
rownames(data.b6) <- gene_info[genes_keep, ]$symbol
rownames(data.jf1) <- gene_info[genes_keep, ]$symbol

names(gene_info) <- gene_info$symbol

data.all <- data.all[,rownames(metadata)]
data.b6 <- data.b6[,rownames(metadata)]
data.jf1 <- data.jf1[,rownames(metadata)]

# parse into SCE object (it's not single cell data, but still useful)
sce_dataset1 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.jf1, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata), rowData = gene_info)

sce_dataset1$Experiment <- "Experiment1"


#### generate invivo dataset 2

all.files <- list.files("../../invivo_data/processed_data_rtta_exp2/quant/", pattern = ".counts$", full.names = T)
all.data <- lapply(all.files, function(x){
  test <- read_delim(x, delim = "\t", skip = 1)
  #if(test[1, 1] == "Geneid"){
  #  test = test[-1, ]
  #}
  test[,c(1, 7)]
})

data <- do.call("cbind", lapply(all.data, function(x){x[,2]}))
rownames(data) <- all.data[[1]][,1]$Geneid

# add chromosome information to genes
library("EnsDb.Mmusculus.v79")
gene_info <- data.frame(genes(EnsDb.Mmusculus.v79))
rownames(gene_info) <- gene_info$gene_id
gene_info <- gene_info[!duplicated(gene_info$symbol), ]

# we're dropping duplicate ensembl ids here, which might not be ideal, 
# but there are no affected expressed genes on the X
genes_keep <- intersect(rownames(data), gene_info$gene_id)
gene_info <- gene_info[genes_keep, ]

gene_info <- makeGRangesFromDataFrame(gene_info, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = T)

# read in metadata
metadata <- data.frame(
  row.names = gsub("_Gall_sorted.bam", "", colnames(data)[grepl("_Gall_sorted.bam", colnames(data))]),  
  Sample = gsub("_Gall_sorted.bam", "", colnames(data)[grepl("_Gall_sorted.bam", colnames(data))]),
  Sex = rep("Unknown", ncol(data) / 3)
) %>% 
  mutate(ShortSampleName = gsub("(AACMTNTM5_EmbryoJF1Rtta2_23s001680-1-1_Kneuss_)|(AACYJJ7M5_E45doxSept23_23s002684-1-1_loda_)", "", Sample)) %>%
  mutate(Experiment = ifelse(grepl("lane1s", Sample), "Experiment2", "Experiment3"))

data.all <- data[genes_keep, grepl("Gall", colnames(data))]
data.b6 <- data[genes_keep,grepl("C57BL-6J", colnames(data))]
data.jf1 <- data[genes_keep,grepl("JF1_MSJ", colnames(data))]

colnames(data.all) <- gsub("_Gall_sorted.bam", "", colnames(data.all))
colnames(data.b6) <- gsub("_C57BL-6J_sorted.bam", "", colnames(data.b6))
colnames(data.jf1) <- gsub("_JF1_MSJ_sorted.bam", "", colnames(data.jf1))

rownames(data.all) <- gene_info[genes_keep, ]$symbol
rownames(data.b6) <- gene_info[genes_keep, ]$symbol
rownames(data.jf1) <- gene_info[genes_keep, ]$symbol

names(gene_info) <- gene_info$symbol

data.all <- data.all[,rownames(metadata)]
data.b6 <- data.b6[,rownames(metadata)]
data.jf1 <- data.jf1[,rownames(metadata)]

# parse into SCE object (it's not single cell data, but still useful)
sce_dataset2 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.jf1, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata), rowData = gene_info)

rtta_counts <- read_csv("/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/rtta_alignment/output/exp2/merged.rtta.counts")
rtta_counts_2 <- read_csv("/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/rtta_alignment/output/exp3/merged.rtta.counts")
rtta_counts_3 <- read_csv("/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/rtta_alignment/output/exp4/merged.rtta.counts")

rtta_counts <- rbind(rtta_counts, rtta_counts_2, rtta_counts_3)

sce_dataset2$RttaCounts <- (rtta_counts %>% pull(RttaCounts, name = Sample))[colnames(sce_dataset2)]
sce_dataset2$RttaGenotype <- ifelse(sce_dataset2$RttaCounts > 50, "Rtta+", "Rtta-")
sce_dataset2$Sex <- "female"

#### Merge datasets

sce_merged <- cbind(sce_dataset1, sce_dataset2)

#### 

escapee_annotation <- readRDS("./ProcessedData/metadata_escape.rds")
rowData(sce_merged)$EscapeAnnotation <- escapee_annotation$final_status_2[match(rowData(sce_merged)$gene_id, escapee_annotation$ENSEMBL_v102)]

#### save

saveRDS(sce_merged, "/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/allele-specific_expression/xist_project/ProcessedData/merged_dataset_invivo.rds")

