### Preprocessing and generation of expression objects

### ### ### ### ### ### ### ### ### ### ### 
# read dataset1
### ### ### ### ### ### ### ### ### ### ### 

library(tidyverse)
library(SingleCellExperiment)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../../../")

source("./xist_project/Scripts/General/reuse_functions.R")

data <- read.csv("./Data/RNASeq/December_2021/table_raw_counts.txt", sep = "\t")

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
metadata <- readxl::read_excel("./Data/RNASeq/December_2021/RNAseq-samples-CL30_31.xlsx")
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Allele", "ndTreatment", "ndAux", "timeAux", "ndDox", "timeDox", 
                        "washout", "timeWO", "WOwithAux", "Replicates")

metadata$Guide = NA
metadata$Experiment = "December2021"

metadata <- metadata[,!colnames(metadata) %in% c("DayOfExperiment", "Replicates")]
metadata$Replicate = "Rep1"

metadata <- metadata %>% 
  mutate(SampleName = gsub(('\"'), "", gsub("^[^CL]*", "", SampleName)))

data <- data[,metadata$SampleName]

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "allelic reads + unassigned reads")
metadata_reduced$SampleName <- gsub("_Gall_sorted.bam", "", metadata_reduced$SampleName)

data.all <- data[genes_keep, metadata$Allele == "allelic reads + unassigned reads"]
data.b6 <- data[genes_keep, metadata$Allele == "C57BL.6J"]
data.cast <- data[genes_keep, metadata$Allele == "CAST.EiJ"]

colnames(data.all) <- gsub("_Gall_sorted.bam", "", colnames(data.all))
colnames(data.b6) <- gsub("_C57BL.6J_sorted.bam", "", colnames(data.b6))
colnames(data.cast) <- gsub("_CAST.EiJ_sorted.bam", "", colnames(data.cast))

data.all <- data.all[,metadata_reduced$SampleName]
data.b6 <- data.b6[,metadata_reduced$SampleName]
data.cast <- data.cast[,metadata_reduced$SampleName]

rownames(data.all) <- genes[genes_keep, ]$symbol
rownames(data.b6) <- genes[genes_keep, ]$symbol
rownames(data.cast) <- genes[genes_keep, ]$symbol

names(genes) <- genes$symbol

# parse into SCE object
sce_dataset1 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)

### ### ### ### ### ### ### ### ### ### ### 
# read dataset2
### ### ### ### ### ### ### ### ### ### ### 

data <- read.csv("./Data/RNASeq/October_2021/dataE6_Dox_rep1_table_raw_counts.txt", sep = "\t")

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
metadata <- readxl::read_excel("./Data/RNASeq/October_2021/Metadata_RNASeq_October2021_E6rep1.xlsx")
metadata <- metadata[-c(1:2), -1]
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Addition", "Guide", "Replicate", "ndTreatment", "ndAux", "timeAux", 
                        "ndDox", "timeDox", "washout", "timeWO", "WOwithAux", "reps", "Allele")

metadata <- metadata %>% 
  mutate(SampleName = gsub(('\"'), "", gsub("^[^CL]*", "", SampleName)))

metadata$Experiment = "October2021"
metadata$Replicate = "Rep0"

metadata <- metadata[,colnames(colData(sce_dataset1))]
metadata$SampleName <- gsub("-", ".", metadata$SampleName)
colnames(data) <- gsub("RNA_NP", "", colnames(data))
data <- data[,metadata$SampleName]

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "allelic reads + unassigned reads")
metadata_reduced$SampleName <- gsub("_Gall_sorted.bam", "", metadata_reduced$SampleName)

data.all <- data[genes_keep, metadata$Allele == "allelic reads + unassigned reads"]
data.b6 <- data[genes_keep, metadata$Allele == "C57BL.6J"]
data.cast <- data[genes_keep, metadata$Allele == "CAST.EiJ"]

colnames(data.all) <- gsub("_Gall_sorted.bam", "", colnames(data.all))
colnames(data.b6) <- gsub("_C57BL.6J_sorted.bam", "", colnames(data.b6))
colnames(data.cast) <- gsub("_CAST.EiJ_sorted.bam", "", colnames(data.cast))

data.all <- data.all[,metadata_reduced$SampleName]
data.b6 <- data.b6[,metadata_reduced$SampleName]
data.cast <- data.cast[,metadata_reduced$SampleName]

rownames(data.all) <- genes[genes_keep, ]$symbol
rownames(data.b6) <- genes[genes_keep, ]$symbol
rownames(data.cast) <- genes[genes_keep, ]$symbol

names(genes) <- genes$symbol

# parse into SCE object 
sce_dataset2 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)

### ### ### ### ### ### ### ### ### ### ### 
# read dataset3
### ### ### ### ### ### ### ### ### ### ### 

data <- read.csv("./Data/RNASeq/March_2022/Xist_data_raw_counts_Mar2022.txt", sep = "\t")

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
metadata <- readxl::read_excel("./Data/RNASeq/March_2022/Metadata_RNASeq_March2022_E6rep2andKd.xlsx")
metadata <- metadata[-c(1:2), -1]
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Addition", "Guide", "Replicate", "ndTreatment", "ndAux", "timeAux", 
                        "ndDox", "timeDox", "washout", "timeWO", "WOwithAux", "reps", "Allele")

metadata <- metadata %>% 
  mutate(SampleName = gsub(('\"'), "", gsub("^[^CL]*", "", SampleName)))

metadata$Experiment = "March2022"
metadata$Replicate = "Rep1"
metadata <- metadata[,colnames(colData(sce_dataset1))]
metadata$SampleName <- gsub("-", ".", metadata$SampleName)
colnames(data) <- gsub("RNA_NP", "", colnames(data))
data <- data[,metadata$SampleName]

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "allelic reads + unassigned reads")
metadata_reduced$SampleName <- gsub("_Gall_sorted.bam", "", metadata_reduced$SampleName)

data.all <- data[genes_keep, metadata$Allele == "allelic reads + unassigned reads"]
data.b6 <- data[genes_keep, metadata$Allele == "C57BL.6J"]
data.cast <- data[genes_keep, metadata$Allele == "CAST.EiJ"]

colnames(data.all) <- gsub("_Gall_sorted.bam", "", colnames(data.all))
colnames(data.b6) <- gsub("_C57BL.6J_sorted.bam", "", colnames(data.b6))
colnames(data.cast) <- gsub("_CAST.EiJ_sorted.bam", "", colnames(data.cast))

data.all <- data.all[,metadata_reduced$SampleName]
data.b6 <- data.b6[,metadata_reduced$SampleName]
data.cast <- data.cast[,metadata_reduced$SampleName]

rownames(data.all) <- genes[genes_keep, ]$symbol
rownames(data.b6) <- genes[genes_keep, ]$symbol
rownames(data.cast) <- genes[genes_keep, ]$symbol

names(genes) <- genes$symbol

# parse into SCE object
sce_dataset3 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)

### ### ### ### ### ### ### ### ### ### ### 
# read dataset4
### ### ### ### ### ### ### ### ### ### ### 

data_raw <- read.csv("./Data/RNASeq/May_2022_clones/RNA-seq_NPC_total_counts.txt", sep = "\t") %>%
  column_to_rownames("Geneid")

data <- read.csv("./Data/RNASeq/May_2022_clones/RNA-seq_NPC_allelic_counts.txt", sep = "\t")

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
metadata <- readxl::read_excel("./Data/RNASeq/March_2022/Metadata_RNASeq_March2022_E6rep2andKd.xlsx")
metadata <- metadata[-c(1:2), -1]
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Addition", "Guide", "Replicate", "ndTreatment", "ndAux", "timeAux", 
                        "ndDox", "timeDox", "washout", "timeWO", "WOwithAux", "reps", "Allele")

metadata <- metadata %>% 
  mutate(SampleName = gsub(('\"'), "", gsub("^[^CL]*", "", SampleName)))

metadata <- data.frame(
  SampleName = colnames(data), 
  CorrectedSampleName = gsub("RNA_NPC_|_C57BL.6J|_CAST.EiJ", "",  colnames(data)), 
  Clone = gsub("RNA_NPC_|_C57BL.6J|_CAST.EiJ|_rep.", "",  colnames(data)), 
  Allele = unlist(lapply(strsplit(colnames(data), "_"), function(x){x[[4]]})), 
  ndTreatment = NA, ndAux = NA, timeAux = NA, ndDox = NA, timeDox = NA, washout = NA, timeWO = NA, 
  WOwithAux = NA, Guide = NA, 
  row.names = colnames(data)
) %>%
  mutate(timeWO = ifelse(is.na(timeWO), "NA", timeWO))

metadata$Experiment = "May2022"
metadata$Replicate = "RepX"

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "C57BL.6J")
metadata_reduced$SampleName <- gsub("_C57BL.6J", "", metadata_reduced$SampleName)

data.all <- data_raw[genes_keep, ]
data.b6 <- data[genes_keep, metadata$Allele == "C57BL.6J"]
data.cast <- data[genes_keep, metadata$Allele == "CAST.EiJ"]

colnames(data.all) <- gsub("_Gall", "", colnames(data.all))
colnames(data.b6) <- gsub("_C57BL.6J", "", colnames(data.b6))
colnames(data.cast) <- gsub("_CAST.EiJ", "", colnames(data.cast))

data.all <- data.all[,metadata_reduced$SampleName]
data.b6 <- data.b6[,metadata_reduced$SampleName]
data.cast <- data.cast[,metadata_reduced$SampleName]

rownames(data.all) <- genes[genes_keep, ]$symbol
rownames(data.b6) <- genes[genes_keep, ]$symbol
rownames(data.cast) <- genes[genes_keep, ]$symbol

names(genes) <- genes$symbol

# we need to adjust cell lines here
b6_cell_lines <- c("B1", "C5", "H4")
cast_cell_lines <- c("CL30", "E6", "JTG")

# these will be the inactive ones
data.inactive <- cbind(data.cast[,grepl("B1|C5|H4", colnames(data.cast))], data.b6[,grepl("CL30|E6|JTG", colnames(data.b6))])
data.active <- cbind(data.b6[,grepl("B1|C5|H4", colnames(data.b6))], data.cast[,grepl("CL30|E6|JTG", colnames(data.cast))])

data.inactive <- data.inactive[,colnames(data.all)]
data.active <- data.active[,colnames(data.all)]

rownames(metadata_reduced) <- metadata_reduced$SampleName

# parse into SCE object
sce_dataset4 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.active, "counts_inactive" = data.inactive), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)

### ### ### ### ### ### ### ### ### ### ### 
# read dataset5
### ### ### ### ### ### ### ### ### ### ### 

data_raw <- read.csv("./Data/RNASeq/September_2022/table_raw_counts.txt", sep = "\t")

# isolate genes information
genes <- data_raw[,1:2]
rownames(data_raw) <- genes$Geneid

data <- data_raw[,-c(1:2)]

rownames(data_raw) <- genes$Geneid

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
metadata <- readxl::read_excel("./Data/RNASeq/September_2022/Metadata_Xistlong_Rep1.xlsx")
metadata <- metadata[-c(1:2), -1]
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Addition", "Guide", "Replicate", "ndTreatment", "ndAux", "timeAux", 
                        "ndDox", "timeDox", "washout", "timeWO", "WOwithAux", "reps", "Allele")

metadata <- metadata %>% 
  mutate(SampleName = gsub("-", ".", metadata$SampleName))

metadata$Experiment = "September2022"
metadata$Replicate = "Rep2"

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "C57BL.6J") %>%
  mutate(SampleName = gsub("_C57BL.6J", "", SampleName)) %>%
  dplyr::select(-c("reps", "Addition"))

data.all <- data[genes_keep, ]
data.b6 <- data[genes_keep, metadata$Allele == "C57BL.6J"]
data.cast <- data[genes_keep, metadata$Allele == "CAST.EiJ"]

colnames(data.all) <- gsub("_Gall", "", colnames(data))
colnames(data.b6) <- gsub("_C57BL.6J", "", colnames(data.b6))
colnames(data.cast) <- gsub("_CAST.EiJ", "", colnames(data.cast))

data.all <- data.all[,metadata_reduced$SampleName]
data.b6 <- data.b6[,metadata_reduced$SampleName]
data.cast <- data.cast[,metadata_reduced$SampleName]

rownames(data.all) <- genes[genes_keep, ]$symbol
rownames(data.b6) <- genes[genes_keep, ]$symbol
rownames(data.cast) <- genes[genes_keep, ]$symbol

names(genes) <- genes$symbol

# parse into SCE object
sce_dataset5 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6),
                                     colData = DataFrame(metadata_reduced), rowData = genes)

### ### ### ### ### ### ### ### ### ### ### 
# read dataset6
### ### ### ### ### ### ### ### ### ### ### 

data_raw <- read.csv("./Data/RNASeq/October_2022/table_raw_counts.txt", sep = "\t")

# isolate genes information
genes <- data_raw[,1:2]
rownames(data_raw) <- genes$Geneid

data <- data_raw[,-c(1:2)]

rownames(data_raw) <- genes$Geneid

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
metadata <- readxl::read_excel("./Data/RNASeq/October_2022/Metadata_RNASeq_Sep2022.xlsx")
metadata <- metadata[-c(1:2), -1]
metadata <- metadata[,colnames(metadata) != "no. of days of WO"] # inconsistent with previous annotations, recover this from period of WO
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Addition", "Guide", "Replicate", "ndTreatment", "ndAux", "timeAux", 
                        "ndDox", "timeDox", "washout", "timeWO", "WOwithAux", "reps", "Allele")

metadata <- metadata %>% 
  mutate(SampleName = gsub("-", ".", metadata$SampleName))

# Exclude unrelated samples
data <- data[ , !is.na(metadata$Clone)]
metadata <- metadata[!is.na(metadata$Clone), ]

metadata$Experiment = "October2022"
metadata$Replicate <- ifelse(metadata$Clone == "E6", "Rep3", "Rep2")

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "C57BL.6J") %>%
  mutate(SampleName = gsub("_C57BL.6J", "", SampleName)) %>%
  dplyr::select(-c("reps", "Addition")) 

data.all <- data[genes_keep, ]
data.b6 <- data[genes_keep, metadata$Allele == "C57BL.6J"]
data.cast <- data[genes_keep, metadata$Allele == "CAST.EiJ"]

colnames(data.all) <- gsub("_Gall", "", colnames(data))
colnames(data.b6) <- gsub("_C57BL.6J", "", colnames(data.b6))
colnames(data.cast) <- gsub("_CAST.EiJ", "", colnames(data.cast))

data.all <- data.all[,metadata_reduced$SampleName]
data.b6 <- data.b6[,metadata_reduced$SampleName]
data.cast <- data.cast[,metadata_reduced$SampleName]

rownames(data.all) <- genes[genes_keep, ]$symbol
rownames(data.b6) <- genes[genes_keep, ]$symbol
rownames(data.cast) <- genes[genes_keep, ]$symbol

names(genes) <- genes$symbol

# parse into SCE object
sce_dataset6 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)

### ### ### ### ### ### ### ### ### ### ### 
# read dataset7
### ### ### ### ### ### ### ### ### ### ### 
data_raw <- read.csv("./Data/RNASeq/April_2023/quant/table_raw_counts.txt", sep = "\t")

# isolate genes information
genes <- data_raw[,1:2]
rownames(data_raw) <- genes$Geneid

data <- data_raw[,-c(1:2)]

rownames(data_raw) <- genes$Geneid

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
metadata <- readxl::read_excel("./Data/RNASeq/April_2023/Metadata_RNASeq_April2023.xlsx")
metadata <- metadata[-c(1:2), -1]
metadata <- metadata[,colnames(metadata) != "no. of days of WO"] # inconsistent with previous annotations, recover this from period of WO
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Addition", "Guide", "Replicate", "ndTreatment", "ndAux", "timeAux", 
                        "ndDox", "timeDox", "washout", "timeWO", "WOwithAux", "reps", "Allele")

metadata <- metadata %>% 
  mutate(SampleName = gsub("-", ".", metadata$SampleName))

# Exclude unrelated samples
data <- data[ , !is.na(metadata$Clone)]
metadata <- metadata[!is.na(metadata$Clone), ]

metadata$Experiment = "April2023"
metadata$Replicate <- ifelse(metadata$Clone == "E6", "Rep4", "Rep2")

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "C57BL.6J") %>%
  mutate(SampleName = gsub("_C57BL.6J", "", SampleName)) %>%
  dplyr::select(-c("reps", "Addition")) 

data.all <- data[genes_keep, ]
data.b6 <- data[genes_keep, metadata$Allele == "C57BL.6J"]
data.cast <- data[genes_keep, metadata$Allele == "CAST.EiJ"]

colnames(data.all) <- gsub("_Gall", "", colnames(data))
colnames(data.b6) <- gsub("_C57BL.6J", "", colnames(data.b6))
colnames(data.cast) <- gsub("_CAST.EiJ", "", colnames(data.cast))

data.all <- data.all[,metadata_reduced$SampleName]
data.b6 <- data.b6[,metadata_reduced$SampleName]
data.cast <- data.cast[,metadata_reduced$SampleName]

rownames(data.all) <- genes[genes_keep, ]$symbol
rownames(data.b6) <- genes[genes_keep, ]$symbol
rownames(data.cast) <- genes[genes_keep, ]$symbol

names(genes) <- genes$symbol

# parse into SCE object
sce_dataset7 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)


### ### ### ### ### ### ### ### ### ### ### 
# read dataset8
### ### ### ### ### ### ### ### ### ### ### 
data_raw <- read.csv("./Data/RNASeq/October_2023/table_raw_counts.txt", sep = "\t")

# isolate genes information
genes <- data_raw[,1:2]
rownames(data_raw) <- genes$Geneid

data <- data_raw[,-c(1:2)]

rownames(data_raw) <- genes$Geneid

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
metadata <- readxl::read_excel("./Data/RNASeq/October_2023/Metadata_RNASeq_NextSeq2000_09_CTCFNPCdegronandlastXistsampels.xlsx")
metadata <- metadata[,c(-1, -16)]
metadata <- metadata[,colnames(metadata) != "no. of days of WO"] # inconsistent with previous annotations, recover this from period of WO
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Addition", "Guide", "Replicate", "ndTreatment", "ndAux", "timeAux", 
                        "ndDox", "timeDox", "washout", "timeWO", "WOwithAux", "reps", "Allele")

metadata <- metadata %>% 
  mutate(SampleName = gsub("-", ".", metadata$SampleName))

metadata$Experiment = "October2023"
metadata$Replicate <- ifelse(metadata$Clone == "E6", "Rep5", "Rep2")

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "C57BL.6J") %>%
  mutate(SampleName = gsub("_C57BL.6J", "", SampleName)) %>%
  dplyr::select(-c("reps", "Addition")) 

data <- data[ , metadata$SampleName]

data.all <- data[genes_keep, metadata$Allele == "Gall"]
data.b6 <- data[genes_keep, metadata$Allele == "C57BL.6J"]
data.cast <- data[genes_keep, metadata$Allele == "CAST.EiJ"]

colnames(data.all) <- gsub("_Gall", "", colnames(data.all))
colnames(data.b6) <- gsub("_C57BL.6J", "", colnames(data.b6))
colnames(data.cast) <- gsub("_CAST.EiJ", "", colnames(data.cast))

data.all <- data.all[,metadata_reduced$SampleName]
data.b6 <- data.b6[, metadata_reduced$SampleName]
data.cast <- data.cast[, metadata_reduced$SampleName]

rownames(data.all) <- genes[genes_keep, ]$symbol
rownames(data.b6) <- genes[genes_keep, ]$symbol
rownames(data.cast) <- genes[genes_keep, ]$symbol

names(genes) <- genes$symbol

# parse into SCE object
sce_dataset8 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)


### ### ### ### ### ### ### ### ### ### ### 
# read dataset9
### ### ### ### ### ### ### ### ### ### ### 

data_raw <- read.csv("./Data/RNASeq/May_2025/table_raw_counts.txt", sep = "\t")
data_raw_reseq <- read.csv("./Data/RNASeq/May_2025/table_raw_counts_reseq.txt", sep = "\t")
data_raw <- cbind(data_raw, data_raw_reseq[ , 3:8])

# isolate genes information
genes <- data_raw[,1:2]
rownames(data_raw) <- genes$Geneid

data <- data_raw[,-c(1:2)]

rownames(data_raw) <- genes$Geneid

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
metadata <- readxl::read_excel("./Data/RNASeq/May_2025/MetaDataSheet_Rebuttal_RNAseq_052025_Jasper.xlsx")
metadata <- metadata[, -1]
metadata <- metadata[,colnames(metadata) != "no. of days of WO"] # inconsistent with previous annotations, recover this from period of WO
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Addition", "Guide", "Replicate", "ndTreatment", "ndAux", "timeAux", 
                        "ndDox", "timeDox", "washout", "timeWO", "WOwithAux", "reps", "Allele")

metadata$Experiment = "May2025"

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "C57BL.6J") %>%
  dplyr::select(-c("reps", "Addition"))

data.all <- data[genes_keep, grepl("Gall", colnames(data))]
data.b6 <- data[genes_keep, grepl("C57BL", colnames(data))]
data.cast <- data[genes_keep, grepl("CAST", colnames(data))]

colnames(data.all) <- gsub("_Gall_sorted.bam", "", colnames(data.all))
colnames(data.b6) <- gsub("_C57BL.6J_sorted.bam", "", colnames(data.b6))
colnames(data.cast) <- gsub("_CAST.EiJ_sorted.bam", "", colnames(data.cast))

data.all <- data.all[,metadata_reduced$SampleName]
data.b6 <- data.b6[,metadata_reduced$SampleName]
data.cast <- data.cast[,metadata_reduced$SampleName]

rownames(data.all) <- genes[genes_keep, ]$symbol
rownames(data.b6) <- genes[genes_keep, ]$symbol
rownames(data.cast) <- genes[genes_keep, ]$symbol

names(genes) <- genes$symbol

# parse into SCE object
sce_dataset9 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)

# there is one wrong sample, remove that here:
sample_out <- c("CL307dDox7dWOR1", "E614dDox14dWOR2")
sce_dataset9 <- sce_dataset9[ , !sce_dataset9$SampleName %in% sample_out ]

### merge datasets
dataset_complete <- do.call("cbind", list(sce_dataset1, sce_dataset2, sce_dataset3, sce_dataset4, sce_dataset5, sce_dataset6, sce_dataset7, 
                                          sce_dataset8, sce_dataset9))

### add escapee annotation

# Read known information about escapee-status. 
# This is based on a literature survey and classifies genes into constitutive, facultative and non-escaping apart from Xist

escapee_annotation <- readr::read_csv("./ProcessedData/ListOfEscapeeFromEdithLabUpdated.csv")
library("EnsDb.Mmusculus.v79")
symbols <- mapIds(EnsDb.Mmusculus.v79, keys = escapee_annotation$ENSEMBL_v102, keytype = "GENEID", column="SYMBOL")
escapee_annotation$symbol <- symbols
convert_status_names <- setNames(c("constitutive", "facultative", "variable"), c("E", "V", "S"))
convert_status_names2 <- setNames(c("constitutive", "facultative", "silenced / variable"), c("E", "V", "S"))

to_add_escaping <- ifelse(is.na(as.numeric(escapee_annotation$Yang2022_Status == "E")), 0, as.numeric(escapee_annotation$Yang2022_Status == "E")) + 
  ifelse(is.na(as.numeric(escapee_annotation$Bowness22_Status == "E")), 0, as.numeric(escapee_annotation$Bowness22_Status == "E"))
to_add_silenced <- ifelse(is.na(as.numeric(escapee_annotation$Yang2022_Status == "S")), 0, as.numeric(escapee_annotation$Yang2022_Status == "S")) + 
  ifelse(is.na(as.numeric(escapee_annotation$Bowness22_Status == "S")), 0, as.numeric(escapee_annotation$Bowness22_Status == "S"))

escapee_annotation$`nb of times escapee 2` <- escapee_annotation$`nb of times escapee` + to_add_escaping
escapee_annotation$`nb of times silenced 2` <- escapee_annotation$`nb of times silenced` + to_add_silenced

escapee_annotation <- escapee_annotation %>%
  mutate(detected = `nb of times escapee 2` + `nb of times silenced 2`) %>%
  mutate(ratio = `nb of times escapee 2` / detected) %>% 
  mutate(final_status_2 = case_when(
    ratio > .5 & detected > 3 ~ "Constitutive", 
    (ratio <= .5 & ratio > 0) | (ratio > 0 & detected <= 3) ~ "Facultative", 
    .default = "NPC-specific"
  ))

escapee_annotation <- escapee_annotation[!duplicated(escapee_annotation$symbol), ]

add_to_sce <- setNames(rep(NA, nrow(dataset_complete)), rownames(dataset_complete))
add_to_sce[escapee_annotation$symbol] <- escapee_annotation$final_status_2
add_to_sce[is.na(add_to_sce)] <- "NPC-specific"

rowData(dataset_complete) <- cbind(rowData(dataset_complete), "EscapeAnnotation" = as.character(add_to_sce[-length(add_to_sce)]))

### Plot this annotation
escapee_annotation <- readr::read_csv("./ProcessedData/ListOfEscapeeFromEdithLabUpdated.csv")
symbols <- mapIds(EnsDb.Mmusculus.v79, keys = escapee_annotation$ENSEMBL_v102, keytype = "GENEID", column="SYMBOL")
escapee_annotation$symbol <- symbols

# check number of papers in the meta-analysis
escapee_annotation[ , grepl("Status", colnames(escapee_annotation))] %>% ncol()

x_genes <- genes(EnsDb.Mmusculus.v79)

escapee_annotation <- escapee_annotation[escapee_annotation$ENSEMBL_v102 %in% names(x_genes), ]
escapee_annotation$position <- start(x_genes[escapee_annotation$ENSEMBL_v102, ])

escapee_annotation_status <- escapee_annotation[,grepl("Status|symbol|position", colnames(escapee_annotation))]

escapee_annotation_status %>%
  pivot_longer(-c(symbol, position)) %>%
  mutate(value = c("E" = "Escaping Gene", "S" = "Silent Gene", "NA" = "Not measured")[.$value]) %>%
  #mutate(value = factor(value, levels = c("Silent Gene", "Escapeing Gene", "Not measured"))) %>%
  ggplot(aes(x = reorder(symbol, position), y = name, fill = value)) + geom_tile() + 
    theme_paper() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
  scale_fill_manual(values = c("red", "grey", "black")) + 
  ylab("Study") + xlab("Genes ordered by chromosomal coordinate") + 
  labs(fill = "")
ggsave("./Plots//FigS1/escape_annotation_heatmap.pdf", width = 8, height = 8)

escape_colors <- setNames(c("grey", "darkgreen", "orange"), nm = c("silenced / variable", "facultative", "constitutive"))
escape_colors_other_name <- setNames(c("lightblue", "grey", "darkgreen", "orange"), nm = c("NA", "S", "V", "E"))
escape_colors_other_name_2 <- setNames(c("lightblue", "grey", "darkgreen", "orange"), nm = c("NA", "NPC-specific", "Facultative", "Constitutive"))

escapee_annotation <- escapee_annotation %>%
  mutate("escape_status" = c("NA" = "not measured", "S" = "silenced / variable", "V" = "facultative", "E" = "constitutive")[.$`final status`]) 

escapee_annotation %>%
  ggplot(aes(x = `nb of times escapee`, y = `nb of times silenced`, col = escape_status)) + geom_jitter() + 
    theme_paper() + 
    ggrepel::geom_text_repel(data = escapee_annotation[escapee_annotation$`final status` == "E", ], aes(label = symbol, col = escape_status)) + 
    xlab("Number of studies in which gene escapes") + ylab("Number of studies in which gene is silenced") + 
    scale_color_manual(values = escape_colors)
ggsave("./Plots/FigS1/escape_annotation_scatter.pdf", width = 8, height = 8)

to_add_escaping <- ifelse(is.na(as.numeric(escapee_annotation$Yang2022_Status == "E")), 0, as.numeric(escapee_annotation$Yang2022_Status == "E")) + 
  ifelse(is.na(as.numeric(escapee_annotation$Bowness22_Status == "E")), 0, as.numeric(escapee_annotation$Bowness22_Status == "E"))
to_add_silenced <- ifelse(is.na(as.numeric(escapee_annotation$Yang2022_Status == "S")), 0, as.numeric(escapee_annotation$Yang2022_Status == "S")) + 
  ifelse(is.na(as.numeric(escapee_annotation$Bowness22_Status == "S")), 0, as.numeric(escapee_annotation$Bowness22_Status == "S"))

escapee_annotation$`nb of times escapee 2` <- escapee_annotation$`nb of times escapee` + to_add_escaping
escapee_annotation$`nb of times silenced 2` <- escapee_annotation$`nb of times silenced` + to_add_silenced

escapee_annotation <- escapee_annotation %>%
  mutate(detected = `nb of times escapee 2` + `nb of times silenced 2`) %>%
  mutate(ratio = `nb of times escapee 2` / detected) %>% 
  mutate(final_status_2 = case_when(
    ratio > .5 & detected > 3 ~ "Constitutive", 
    (ratio <= .5 & ratio > 0) | (ratio > 0 & detected <= 3) ~ "Facultative", 
    .default = "NPC-specific"
  ))

escapee_annotation %>%
  mutate(detected = `nb of times escapee 2` + `nb of times silenced 2`) %>%
  mutate(ratio = `nb of times escapee 2` / detected) %>% 
  {
    ggplot(data = ., aes(x = detected, y = ratio, fill = final_status_2)) + 
      geom_jitter(pch = 21, size = 3, width = .1, height = .005) + 
      theme_paper() + 
      ggrepel::geom_text_repel(data = .[.$`final status` == "E", ], aes(label = symbol, col = final_status_2)) + 
      xlab("Number of studies in which gene is analyzed") + ylab("Fraction of studies in which gene escapes") + 
      scale_fill_manual(values = escape_colors_other_name_2) + 
      scale_color_manual(values = escape_colors_other_name_2) + 
      annotate(geom = "rect", fill = "orange", xmin = 3, xmax = 16.5, ymin = 0.5, ymax = 1.05, alpha = .2) + 
      annotate(geom = "rect", fill = "darkgreen", xmin = 0, xmax = 16.5, ymin = 0, ymax = 0.5, alpha = .2) + 
      annotate(geom = "rect", fill = "darkgreen", xmin = 0, xmax = 3, ymin = 0.5, ymax = 1.05, alpha = .2)
  }

plot_df_s1h <- escapee_annotation %>%
  mutate(detected = `nb of times escapee 2` + `nb of times silenced 2`) %>%
  mutate(ratio = `nb of times escapee 2` / detected) %>%
  dplyr::select(c(ENSEMBL_v102, symbol, detected, ratio))

saveRDS(plot_df_s1h, "./ProcessedData/source_data/source_data_s1h.rds")

escapee_annotation %>% dplyr::filter(detected > 0) %>% pull(final_status_2) %>% table()

### Add number of days washout as variable
dataset_complete$ndWashout <- unlist(lapply(dataset_complete$timeWO, function(x){
  if (is.na(x) | x == "NO" | x == "NA" | x == "0"){
    return(0)
  } else if(x == "7") {
    return(7)
  } else {
    nd_days <- as.numeric(stringr::str_split(x, ":")[[1]][[2]]) - as.numeric(stringr::str_split(x, ":")[[1]][[1]])
  }
}))

### create condition annotations
dataset_complete$Condition <- paste0("Aux_", dataset_complete$ndAux, "_Dox_", dataset_complete$ndDox, 
                                     "_WO_", dataset_complete$ndWashout, "_WOAux", dataset_complete$WOwithAux)

### Subset on main dataset for paper, remove May2022 (other clones) and crispr experiment (not necessary)
data <- dataset_complete[,!dataset_complete$Guide %in% c("g13", "P3") & dataset_complete$Experiment != "May2022"]

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
  "Aux_0_Dox_14_WO_21_WOAuxNO" = "Dox (14d) - washout (21d)",
  
  "Aux_0_Dox_21_WO_0_WOAuxNO" = "Dox (21d)",
  "Aux_0_Dox_21_WO_7_WOAuxNO" = "Dox (21d) - washout (7d)",
  "Aux_0_Dox_21_WO_14_WOAuxNO" = "Dox (21d) - washout (14d)",
  "Aux_0_Dox_21_WO_21_WOAuxNO" = "Dox (21d) - washout (21d)",
  
  "Aux_2_Dox_0_WO_0_WOAuxNO" = "Aux (2d)", 
  "Aux_2_Dox_0_WO_4_WOAuxNO" = "Aux (2d) - washout (4d)", 
  "Aux_5_Dox_0_WO_0_WOAuxNO" = "Aux (5d)", 
  "Aux_5_Dox_3_WO_0_WOAuxNO" = "Dox (3d), Aux (5d)", 
  "Aux_5_Dox_3_WO_4_WOAuxNO" = "Dox (3d) - washout (4d), Aux (5d)", 
  "Aux_9_Dox_0_WO_0_WOAuxNO" = "Aux (9d)", 
  "Aux_9_Dox_7_WO_0_WOAuxNO" = "Dox (7d), Aux (9d)", 
  "Aux_9_Dox_7_WO_4_WOAuxNO" = "Dox (7d) - washout (4d), Aux (9d)", 
  "Aux_9_Dox_7_WO_4_WOAuxYES" = "Dox (7d) - washout (4d) (with Aux), Aux (9d)",
  "Aux_21_Dox_21_WO_0_WOAuxNO" = "Dox (21d), Aux (21d)"
)

# We use this grouping to subset the data
data$ConditionClean <- unlist(name_conversion_condition[data$Condition])

#
colData(data) <- colData(data)[ , c(
  "SampleName", "CorrectedSampleName", "Clone", "Allele", "ndTreatment", "ndAux", "timeAux", "ndDox", "timeDox", "washout", "timeWO", "WOwithAux",
  "Guide", "Experiment", "ndWashout", "Condition", "ConditionClean", "Replicate"
)]

###

saveRDS(data, "./ProcessedData/merged_dataset_paper.rds")
write.csv(colData(data), "./ProcessedData/merged_dataset_paper_metadata.csv")

saveRDS(dataset_complete, "./ProcessedData/merged_dataset.rds")
write.csv(colData(dataset_complete), "./ProcessedData/merged_dataset_metadata.csv")

### define geneset to use across analysis
# 10 average reads across samples based on the pre-revision dataset, to keep things consistent
data <- readRDS("./ProcessedData/merged_dataset_paper.rds")
genes_keep <- (rowSums(assays(data[,(data$Experiment != "May2025")])[["counts_inactive"]]) + 
                 rowSums(assays(data[,(data$Experiment != "May2025")])[["counts_active"]])) / 
  ncol(data[,(data$Experiment != "May2025")]) > 10

# remove genes with high allelic ratios on the X
data_here <- data[genes_keep, data$ConditionClean == "Control" & data$Experiment != "May2025"]
data_here <- data_here[seqnames(data_here) == "X", ]

# get d-scores across genes + samples
ratios <- assays(data_here)[["counts_inactive"]] / (assays(data_here)[["counts_inactive"]] + assays(data_here)[["counts_active"]])

# Genes with average ASE > 0.8 in the control sample are likely mapping artefacts:
genes_out <- names(rowMeans(ratios, na.rm = T)[rowMeans(ratios, na.rm = T) > 0.8 & rownames(ratios) != "Xist" ])
genes_out_na <- rownames( ratios[!apply(ratios, 1, function(x){!any(is.na(x))}), ] )
genes_out <- c(genes_out, genes_out_na)

genes_keep[genes_out] <- FALSE
saveRDS(genes_keep, "./ProcessedData/genes_keep.rds")

# data_x <- data[ seqnames(rowRanges(data)) == "X", ]
