### Preprocessing and generation of expression objects

### ### ### ### ### ### ### ### ### ### ### 
# read dataset1
### ### ### ### ### ### ### ### ### ### ### 

setwd("~/Desktop/Projects/XChromosome_Antonia/")

data <- read.csv("~/Desktop/Projects/XChromosome_Antonia/Data/December_2021/table_raw_counts.txt", sep = "\t")

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
metadata <- readxl::read_excel("~/Desktop/Projects/XChromosome_Antonia/Data/December_2021/RNAseq-samples-CL30_31.xlsx")
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Allele", "ndTreatment", "ndAux", "timeAux", "ndDox", "timeDox", 
                        "washout", "timeWO", "WOwithAux", "Replicates")

metadata$Guide = NA
metadata$Experiment = "December2021"

metadata <- metadata[,!colnames(metadata) %in% c("DayOfExperiment", "Replicates")]

metadata <- metadata %>% 
  mutate(SampleName = gsub(('\"'), "", gsub("^[^CL]*", "", SampleName))) %>%
  add_column("Condition" = paste0("Aux_", .$ndAux, "_Dox_", .$ndDox, "_WO_", .$washout, "_WOAux", .$WOwithAux))

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

data <- read.csv("~/Desktop/Projects/XChromosome_Antonia/Data/October_2021/dataE6_Dox_rep1_table_raw_counts.txt", sep = "\t")

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
metadata <- readxl::read_excel("~/Desktop/Projects/XChromosome_Antonia/Data/October_2021/Metadata_RNASeq_October2021_E6rep1.xlsx")
metadata <- metadata[-c(1:2), -1]
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Addition", "Guide", "Replicate", "ndTreatment", "ndAux", "timeAux", 
                        "ndDox", "timeDox", "washout", "timeWO", "WOwithAux", "reps", "Allele")

metadata <- metadata %>% 
  mutate(SampleName = gsub(('\"'), "", gsub("^[^CL]*", "", SampleName))) %>%
  add_column("Condition" = paste0("Aux_", .$ndAux, "_Dox_", .$ndDox, "_WO_", .$washout, "_WOAux", .$WOwithAux))

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

data <- read.csv("~/Desktop/Projects/XChromosome_Antonia/Data/March_2022/Xist_data_raw_counts_Mar2022.txt", sep = "\t")

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
metadata <- readxl::read_excel("~/Desktop/Projects/XChromosome_Antonia/Data/March_2022/Metadata_RNASeq_March2022_E6rep2andKd.xlsx")
metadata <- metadata[-c(1:2), -1]
colnames(metadata) <- c("SampleName", "CorrectedSampleName", "Clone", "Addition", "Guide", "Replicate", "ndTreatment", "ndAux", "timeAux", 
                        "ndDox", "timeDox", "washout", "timeWO", "WOwithAux", "reps", "Allele")

metadata <- metadata %>% 
  mutate(SampleName = gsub(('\"'), "", gsub("^[^CL]*", "", SampleName))) %>%
  add_column("Condition" = paste0("Aux_", .$ndAux, "_Dox_", .$ndDox, "_WO_", .$washout, "_WOAux", .$WOwithAux))

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

### merge datasets
dataset_complete <- cbind(cbind(sce_dataset1, sce_dataset2), sce_dataset3)

saveRDS(dataset_complete, "./ProcessedData/merged_dataset.rds")

