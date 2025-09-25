### Preprocessing and generation of expression objects

### ### ### ### ### ### ### ### ### ### ### 
# read astrocyte dataset
### ### ### ### ### ### ### ### ### ### ### 

library(tidyverse)
library(SingleCellExperiment)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../../../")

data <- read.csv("~/Desktop/PhD/Projects/XChromosome_Antonia/Data/November_2024/table_raw_counts.txt", sep = "\t")

# isolate genes information
genes <- data[,1:2]
data <- data[,-c(1:2)]

rownames(data) <- genes$Geneid

# add chromosome information to genes
library("EnsDb.Mmusculus.v79")
gene_info <- data.frame(genes(EnsDb.Mmusculus.v79))
rownames(gene_info) <- gene_info$gene_id
gene_info <- gene_info[!duplicated(gene_info$symbol), ]

genes_keep <- intersect(genes$Geneid, gene_info$gene_id)
gene_info <- gene_info[genes_keep, ]
genes <- genes[genes$Geneid %in% genes_keep, ]

genes <- gene_info[genes$Geneid, ]
genes <- makeGRangesFromDataFrame(genes, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = T)

# read in metadata
metadata <- readxl::read_excel("~/Desktop/PhD/Projects/XChromosome_Antonia/Data/November_2024/metadata_astro_1.xlsx")

metadata$Experiment = "November2024"

data <- data[,metadata$SampleName]

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "All")
metadata_reduced$SampleName <- gsub("_Gall_sorted.bam", "", metadata_reduced$SampleName)

data.all <- data[genes_keep, metadata$Allele == "All"]
data.b6 <- data[genes_keep, metadata$Allele == "C57B6"]
data.cast <- data[genes_keep, metadata$Allele == "CAST"]

colnames(data.all) <- gsub("_Gall_sorted.bam", "", colnames(data.all))
colnames(data.b6) <- gsub("_C57BL.6J_sorted.bam", "", colnames(data.b6))
colnames(data.cast) <- gsub("_CAST.EiJ_sorted.bam", "", colnames(data.cast))

data.all <- data.all[genes_keep, metadata_reduced$SampleName]
data.b6 <- data.b6[genes_keep,metadata_reduced$SampleName]
data.cast <- data.cast[genes_keep,metadata_reduced$SampleName]

rownames(data.all) <- genes[genes_keep, ]$symbol
rownames(data.b6) <- genes[genes_keep, ]$symbol
rownames(data.cast) <- genes[genes_keep, ]$symbol

names(genes) <- genes$symbol

# parse into SCE object
sce_dataset1 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)

### merge datasets
dataset_complete <- sce_dataset1

### add escapee annotation

# Read known information about escapee-status. 
# This is based on a literature survey and classifies genes into constitutive, facultative and non-escaping
# apart from Xist (which isn't really an escapee bc it's only from the Xi)
escapee_annotation <- readxl::read_excel("~/Desktop/PhD/Projects/XChromosome_Project/ProcessedData/ListOfEscapeeFromEdithLab.xlsx") %>%
  dplyr::select(c("ENSEMBL_v102", "final status"))
library("EnsDb.Mmusculus.v79")
symbols <- mapIds(EnsDb.Mmusculus.v79, keys = escapee_annotation$ENSEMBL_v102, keytype = "GENEID", column="SYMBOL")
escapee_annotation$symbol <- symbols
convert_status_names <- setNames(c("constitutive", "facultative", "variable"), c("E", "V", "S"))
convert_status_names2 <- setNames(c("constitutive", "facultative", "silenced / variable"), c("E", "V", "S"))
escapee_annotation$our_status <- unlist(convert_status_names[escapee_annotation$`final status`])
escapee_annotation$our_status2 <- unlist(convert_status_names2[escapee_annotation$`final status`])
escapee_annotation <- escapee_annotation[!duplicated(escapee_annotation$symbol), ]

add_to_sce <- setNames(rep(NA, nrow(dataset_complete)), rownames(dataset_complete))
add_to_sce[escapee_annotation$symbol] <- escapee_annotation$our_status2

rowData(dataset_complete) <- cbind(rowData(dataset_complete), "EscapeAnnotation" = add_to_sce[-length(add_to_sce)])

### Plot this annotation
escapee_annotation <- readxl::read_excel("~/Desktop/PhD/Projects/XChromosome_Project/ProcessedData/ListOfEscapeeFromEdithLab.xlsx")
symbols <- mapIds(EnsDb.Mmusculus.v79, keys = escapee_annotation$ENSEMBL_v102, keytype = "GENEID", column="SYMBOL")
escapee_annotation$symbol <- symbols

x_genes <- genes(EnsDb.Mmusculus.v79)

escapee_annotation <- escapee_annotation[escapee_annotation$ENSEMBL_v102 %in% names(x_genes), ]
escapee_annotation$position <- start(x_genes[escapee_annotation$ENSEMBL_v102, ])

escapee_annotation_status <- escapee_annotation[,grepl("Status|symbol|position", colnames(escapee_annotation))]

escape_colors <- setNames(c("grey", "darkgreen", "orange"), nm = c("silenced / variable", "facultative", "constitutive"))
escape_colors_other_name <- setNames(c("lightblue", "grey", "darkgreen", "orange"), nm = c("NA", "S", "V", "E"))

escapee_annotation <- escapee_annotation %>%
  mutate("escape_status" = c("NA" = "not measured", "S" = "silenced / variable", "V" = "facultative", "E" = "constitutive")[.$`final status`]) 

data <- dataset_complete

dataset_complete$Condition <- gsub("_R1|_R2", "", dataset_complete$CorrectedSampleName)
dataset_complete$ConditionClean <- gsub("_R1|_R2", "", dataset_complete$CorrectedSampleName)

###

saveRDS(dataset_complete, "./ProcessedData/merged_dataset_astrocytes.rds")
write.csv(colData(dataset_complete), "./ProcessedData/merged_dataset_astrocytes_metadata.csv")

### 
