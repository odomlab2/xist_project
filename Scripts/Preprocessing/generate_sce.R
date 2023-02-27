### Preprocessing and generation of expression objects

### ### ### ### ### ### ### ### ### ### ### 
# read dataset1
### ### ### ### ### ### ### ### ### ### ### 

library(tidyverse)
library(SingleCellExperiment)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../../")

source("./Scripts/General/reuse_functions.R")

data <- read.csv("./Data/December_2021/table_raw_counts.txt", sep = "\t")

# isolate genes information
genes <- data[,1:2]
data <- data[,-c(1:2)]

rownames(data) <- genes$Geneid

# add chromosome information to genes
library("EnsDb.Mmusculus.v79")
gene_info <- data.frame(genes(EnsDb.Mmusculus.v79))
rownames(gene_info) <- gene_info$gene_id
gene_info <- gene_info[!duplicated(gene_info$symbol), ]

# we're dropping duplicate ensembl ids here, which might not be ideal, 
# but there are no affected expressed genes on the X
genes_keep <- intersect(genes$Geneid, gene_info$gene_id)
gene_info <- gene_info[genes_keep, ]
genes <- genes[genes$Geneid %in% genes_keep, ]

genes <- gene_info[genes$Geneid, ]
genes <- makeGRangesFromDataFrame(genes, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = T)

# read in metadata
metadata <- readxl::read_excel("./Data/December_2021/RNAseq-samples-CL30_31.xlsx")
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Allele", "ndTreatment", "ndAux", "timeAux", "ndDox", "timeDox", 
                        "washout", "timeWO", "WOwithAux", "Replicates")

metadata$Guide = NA
metadata$Experiment = "December2021"

metadata <- metadata[,!colnames(metadata) %in% c("DayOfExperiment", "Replicates")]

metadata <- metadata %>% 
  mutate(SampleName = gsub(('\"'), "", gsub("^[^CL]*", "", SampleName)))
# metadata <- metadata[match(colnames(data), metadata$SampleName), ]

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

# parse into SCE object (it's not single cell data, but still useful)
sce_dataset1 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)

### ### ### ### ### ### ### ### ### ### ### 
# read dataset2
### ### ### ### ### ### ### ### ### ### ### 

data <- read.csv("./Data/October_2021/dataE6_Dox_rep1_table_raw_counts.txt", sep = "\t")

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
metadata <- readxl::read_excel("./Data/October_2021/Metadata_RNASeq_October2021_E6rep1.xlsx")
metadata <- metadata[-c(1:2), -1]
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Addition", "Guide", "Replicate", "ndTreatment", "ndAux", "timeAux", 
                        "ndDox", "timeDox", "washout", "timeWO", "WOwithAux", "reps", "Allele")

metadata <- metadata %>% 
  mutate(SampleName = gsub(('\"'), "", gsub("^[^CL]*", "", SampleName)))

metadata$Experiment = "October2021"

metadata <- metadata[,colnames(colData(sce_dataset1))]
metadata$SampleName <- gsub("-", ".", metadata$SampleName)
colnames(data) <- gsub("RNA_NP", "", colnames(data))
data <- data[,metadata$SampleName]

# metadata <- metadata[match(colnames(data), metadata$SampleName), ]
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

data <- read.csv("./Data/March_2022/Xist_data_raw_counts_Mar2022.txt", sep = "\t")

# isolate genes information
genes <- data[,1:2]
data <- data[,-c(1:2)]

rownames(data) <- genes$Geneid

### IMPORTANT: The first two samples are swapped. We correct this here. 
# 
# colnames_save <- colnames(data)
# 
# data <- data[,c(4:6, 1:3, 7:ncol(data))]
# colnames(data) <- colnames_save

# double check reference sample
# x_genes <- gene_info[gene_info$seqnames == "X", ]
# x_genes <- rownames(data) [rownames(data) %in% rownames(x_genes)]
# 
# data_x_here <- data[rownames(data) %in% x_genes, ]
# 
# data_b6_here <- data_x_here$RNA_NPC.C5E6_NoDox_r2_C57BL.6J_sorted.bam
# data_cast_here <- data_x_here$RNA_NPC.C5E6_NoDox_r2_CAST.EiJ_sorted.bam
# 
# data_b6_here_x <- data_b6_here[data_b6_here + data_cast_here > 10]
# data_cast_here_x <- data_cast_here[data_b6_here + data_cast_here > 10]
# 
# data.frame(
#   total = data_b6_here_x + data_cast_here_x, 
#   ase = data_b6_here_x / (data_b6_here_x + data_cast_here_x)
# ) %>%
#   ggplot(aes(x = total, y = ase)) + geom_point() + scale_x_log10()

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
metadata <- readxl::read_excel("./Data/March_2022/Metadata_RNASeq_March2022_E6rep2andKd.xlsx")
metadata <- metadata[-c(1:2), -1]
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Addition", "Guide", "Replicate", "ndTreatment", "ndAux", "timeAux", 
                        "ndDox", "timeDox", "washout", "timeWO", "WOwithAux", "reps", "Allele")

metadata <- metadata %>% 
  mutate(SampleName = gsub(('\"'), "", gsub("^[^CL]*", "", SampleName)))

metadata$Experiment = "March2022"
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

# parse into SCE object (it's not single cell data, but still useful)
sce_dataset3 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)

### ### ### ### ### ### ### ### ### ### ### 
# read dataset4
### ### ### ### ### ### ### ### ### ### ### 

data_raw <- read.csv("./Data/May_2022_clones/RNA-seq_NPC_total_counts.txt", sep = "\t") %>%
  column_to_rownames("Geneid")

data <- read.csv("./Data/May_2022_clones/RNA-seq_NPC_allelic_counts.txt", sep = "\t")

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
metadata <- readxl::read_excel("./Data/March_2022/Metadata_RNASeq_March2022_E6rep2andKd.xlsx")
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
)

metadata$Experiment = "May2022"

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

# parse into SCE object (it's not single cell data, but still useful)
sce_dataset4 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.active, "counts_inactive" = data.inactive), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)

### ### ### ### ### ### ### ### ### ### ### 
# read dataset5
### ### ### ### ### ### ### ### ### ### ### 

data_raw <- read.csv("./Data/September_2022/table_raw_counts.txt", sep = "\t")

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
metadata <- readxl::read_excel("./Data/September_2022/Metadata_Xistlong_Rep1.xlsx")
metadata <- metadata[-c(1:2), -1]
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Addition", "Guide", "Replicate", "ndTreatment", "ndAux", "timeAux", 
                        "ndDox", "timeDox", "washout", "timeWO", "WOwithAux", "reps", "Allele")

metadata <- metadata %>% 
  mutate(SampleName = gsub("-", ".", metadata$SampleName))

metadata$Experiment = "September2022"

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "C57BL.6J") %>%
  mutate(SampleName = gsub("_C57BL.6J", "", SampleName)) %>%
  dplyr::select(-c("reps", "Addition", "Replicate"))

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

# parse into SCE object (it's not single cell data, but still useful)
sce_dataset5 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6),
                                     colData = DataFrame(metadata_reduced), rowData = genes)

### ### ### ### ### ### ### ### ### ### ### 
# read dataset6
### ### ### ### ### ### ### ### ### ### ### 

data_raw <- read.csv("./Data/October_2022/table_raw_counts.txt", sep = "\t")

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
metadata <- readxl::read_excel("./Data/October_2022/Metadata_RNASeq_Sep2022.xlsx")
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

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "C57BL.6J") %>%
  mutate(SampleName = gsub("_C57BL.6J", "", SampleName)) %>%
  dplyr::select(-c("reps", "Addition", "Replicate")) 

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

# parse into SCE object (it's not single cell data, but still useful)
sce_dataset6 <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)

### merge datasets
#dataset_complete <- cbind(cbind(cbind(cbind(sce_dataset1, sce_dataset2), sce_dataset3), sce_dataset4), sce_dataset5)
dataset_complete <- do.call("cbind", list(sce_dataset1, sce_dataset2, sce_dataset3, sce_dataset4, sce_dataset5, sce_dataset6))


### add escapee annotation

# Read known information about escapee-status. 
# This is based on a literature survey and classifies genes into constitutive, facultative and non-escaping
# apart from Xist (which isn't really an escapee bc it's only from the Xi)
escapee_annotation <- readxl::read_excel("./ProcessedData/ListOfEscapeeFromEdithLab.xlsx") %>%
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

### Plot this annotation - 
# For the paper we need: 
# Cite all studies, show annotation per study, show how they are defined

escapee_annotation <- readxl::read_excel("./ProcessedData/ListOfEscapeeFromEdithLab.xlsx")
symbols <- mapIds(EnsDb.Mmusculus.v79, keys = escapee_annotation$ENSEMBL_v102, keytype = "GENEID", column="SYMBOL")
escapee_annotation$symbol <- symbols

x_genes <- genes(EnsDb.Mmusculus.v79)

# a small number of genes does not seem to be in Ens79... 
table(escapee_annotation$ENSEMBL_v102 %in% names(x_genes)) # 79 / 2547

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

escapee_annotation <- escapee_annotation %>%
  mutate("escape_status" = c("NA" = "not measured", "S" = "silenced / variable", "V" = "facultative", "E" = "constitutive")[.$`final status`]) 

escapee_annotation%>%
  ggplot(aes(x = `nb of times escapee`, y = `nb of times silenced`, col = escape_status)) + geom_jitter() + 
    theme_paper() + 
    ggrepel::geom_text_repel(data = escapee_annotation[escapee_annotation$`final status` == "E", ], aes(label = symbol, col = escape_status)) + 
    xlab("Number of studies in which gene escapes") + ylab("Number of studies in which gene is silenced") + 
    scale_color_manual(values = escape_colors)
ggsave("./Plots/FigS1/escape_annotation_scatter.pdf", width = 8, height = 8)

### Add number of days washout as variable
dataset_complete$ndWashout <- unlist(lapply(dataset_complete$timeWO, function(x){
  if (is.na(x) | x == "NO" | x == "NA"){
    return(0)
  } else {
    nd_days <- as.numeric(stringr::str_split(x, ":")[[1]][[2]]) - as.numeric(stringr::str_split(x, ":")[[1]][[1]])
  }
}))

dataset_complete$Condition <- paste0("Aux_", dataset_complete$ndAux, "_Dox_", dataset_complete$ndDox, 
                                     "_WO_", dataset_complete$ndWashout, "_WOAux", dataset_complete$WOwithAux)

### Subset on main dataset for paper, remove May2022 (other clones) and crispr experiment (not really necessary)

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

# We use this grouping to subset the data
data$ConditionClean <- unlist(name_conversion_condition[data$Condition])

###

saveRDS(dataset_complete, "./ProcessedData/merged_dataset.rds")
write.csv(colData(dataset_complete), "./ProcessedData/merged_dataset_metadata.csv")

saveRDS(data, "./ProcessedData/merged_dataset_paper.rds")
write.csv(colData(data), "./ProcessedData/merged_dataset_paper_metadata.csv")


### 