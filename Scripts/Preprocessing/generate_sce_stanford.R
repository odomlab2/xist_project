### Preprocessing and generation of expression objects

### ### ### ### ### ### ### ### ### ### ### 
# read dataset1
### ### ### ### ### ### ### ### ### ### ### 

library(tidyverse)
library(SingleCellExperiment)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("../../../")

source("./xist_project/Scripts/General/reuse_functions.R")

# data <- read.csv("./Data/RNASeq/multiclone_geo_upload_processed/multiclone_total_counts.txt", sep = "\t") %>%
#   column_to_rownames("Gene")
# data_genome1 <- read.csv("./Data/RNASeq/multiclone_geo_upload_processed/multiclone_genome1_counts.txt", sep = "\t") %>%
#   column_to_rownames("Gene")
# data_genome2 <- read.csv("./Data/RNASeq/multiclone_geo_upload_processed/multiclone_genome2_counts.txt", sep = "\t") %>%
#   column_to_rownames("Gene")

data <- read.csv("./Data/RNASeq/multiclone_geo_upload_processed/table_raw_counts.txt", sep = "\t")

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

metadata <- data.frame(
  "Clone" = str_extract(colnames(data), "CL_F[0-9]*"), 
  "Allele" = str_extract(colnames(data), "Gall|(CAST\\.EiJ)|(129S1\\.SvImJ)")
)
metadata

metadata_reduced <- metadata %>%
  dplyr::filter(Allele == "Gall")
metadata_reduced$Sample <- gsub("_Gall_sorted.bam", "", metadata_reduced$Clone)
rownames(metadata_reduced) <- metadata_reduced$Clone

data.all <- data[genes_keep, metadata$Allele == "Gall"]
data.b6 <- data[genes_keep, metadata$Allele == "129S1.SvImJ"]
data.cast <- data[genes_keep, metadata$Allele == "CAST.EiJ"]

colnames(data.all) <- gsub("_Gall_sorted.bam", "", colnames(data.all))
colnames(data.b6) <- gsub("_129S1.SvImJ_sorted.bam", "", colnames(data.b6))
colnames(data.cast) <- gsub("_CAST.EiJ_sorted.bam", "", colnames(data.cast))

rownames(data.all) <- genes[genes_keep, ]$symbol
rownames(data.b6) <- genes[genes_keep, ]$symbol
rownames(data.cast) <- genes[genes_keep, ]$symbol

data.all <- data.all[,metadata_reduced$Sample]
data.b6 <- data.b6[,metadata_reduced$Sample]
data.cast <- data.cast[,metadata_reduced$Sample]

names(genes) <- genes$symbol

# parse into SCE object
sce_dataset <- SingleCellExperiment(assays = list("counts" = data.all, "counts_active" = data.cast, "counts_inactive" = data.b6), 
                                     colData = DataFrame(metadata_reduced), rowData = genes)

### add escapee annotation

# Read known information about escapee-status. 
# This is based on a literature survey and classifies genes into constitutive, facultative and non-escaping apart from Xist

# escapee_annotation <- readr::read_csv("./ProcessedData/ListOfEscapeeFromEdithLabUpdated.csv")
# library("EnsDb.Mmusculus.v79")
# symbols <- mapIds(EnsDb.Mmusculus.v79, keys = escapee_annotation$ENSEMBL_v102, keytype = "GENEID", column="SYMBOL")
# escapee_annotation$symbol <- symbols
# convert_status_names <- setNames(c("constitutive", "facultative", "variable"), c("E", "V", "S"))
# convert_status_names2 <- setNames(c("constitutive", "facultative", "silenced / variable"), c("E", "V", "S"))
# 
# to_add_escaping <- ifelse(is.na(as.numeric(escapee_annotation$Yang2022_Status == "E")), 0, as.numeric(escapee_annotation$Yang2022_Status == "E")) + 
#   ifelse(is.na(as.numeric(escapee_annotation$Bowness22_Status == "E")), 0, as.numeric(escapee_annotation$Bowness22_Status == "E"))
# to_add_silenced <- ifelse(is.na(as.numeric(escapee_annotation$Yang2022_Status == "S")), 0, as.numeric(escapee_annotation$Yang2022_Status == "S")) + 
#   ifelse(is.na(as.numeric(escapee_annotation$Bowness22_Status == "S")), 0, as.numeric(escapee_annotation$Bowness22_Status == "S"))
# 
# escapee_annotation$`nb of times escapee 2` <- escapee_annotation$`nb of times escapee` + to_add_escaping
# escapee_annotation$`nb of times silenced 2` <- escapee_annotation$`nb of times silenced` + to_add_silenced
# 
# escapee_annotation <- escapee_annotation %>%
#   mutate(detected = `nb of times escapee 2` + `nb of times silenced 2`) %>%
#   mutate(ratio = `nb of times escapee 2` / detected) %>% 
#   mutate(final_status_2 = case_when(
#     ratio > .5 & detected > 3 ~ "Constitutive", 
#     (ratio <= .5 & ratio > 0) | (ratio > 0 & detected <= 3) ~ "Facultative", 
#     .default = "NPC-specific"
#   ))
# 
# escapee_annotation <- escapee_annotation[!duplicated(escapee_annotation$symbol), ]
# 
# add_to_sce <- setNames(rep(NA, nrow(dataset_complete)), rownames(dataset_complete))
# add_to_sce[escapee_annotation$symbol] <- escapee_annotation$final_status_2
# add_to_sce[is.na(add_to_sce)] <- "NPC-specific"
# 
# rowData(dataset_complete) <- cbind(rowData(dataset_complete), "EscapeAnnotation" = as.character(add_to_sce[-length(add_to_sce)]))
# 
# ### Plot this annotation
# escapee_annotation <- readr::read_csv("./ProcessedData/ListOfEscapeeFromEdithLabUpdated.csv")
# symbols <- mapIds(EnsDb.Mmusculus.v79, keys = escapee_annotation$ENSEMBL_v102, keytype = "GENEID", column="SYMBOL")
# escapee_annotation$symbol <- symbols
# 
# x_genes <- genes(EnsDb.Mmusculus.v79)
# 
# escapee_annotation <- escapee_annotation[escapee_annotation$ENSEMBL_v102 %in% names(x_genes), ]
# escapee_annotation$position <- start(x_genes[escapee_annotation$ENSEMBL_v102, ])
# 
# escapee_annotation_status <- escapee_annotation[,grepl("Status|symbol|position", colnames(escapee_annotation))]
# 
# escapee_annotation_status %>%
#   pivot_longer(-c(symbol, position)) %>%
#   mutate(value = c("E" = "Escaping Gene", "S" = "Silent Gene", "NA" = "Not measured")[.$value]) %>%
#   #mutate(value = factor(value, levels = c("Silent Gene", "Escapeing Gene", "Not measured"))) %>%
#   ggplot(aes(x = reorder(symbol, position), y = name, fill = value)) + geom_tile() + 
#     theme_paper() + 
#     theme(axis.text.x=element_blank(),
#           axis.ticks.x=element_blank()) + 
#   scale_fill_manual(values = c("red", "grey", "black")) + 
#   ylab("Study") + xlab("Genes ordered by chromosomal coordinate") + 
#   labs(fill = "")
# 
# escape_colors <- setNames(c("grey", "darkgreen", "orange"), nm = c("silenced / variable", "facultative", "constitutive"))
# escape_colors_other_name <- setNames(c("lightblue", "grey", "darkgreen", "orange"), nm = c("NA", "S", "V", "E"))
# escape_colors_other_name_2 <- setNames(c("lightblue", "grey", "darkgreen", "orange"), nm = c("NA", "NPC-specific", "Facultative", "Constitutive"))
# 
# escapee_annotation <- escapee_annotation %>%
#   mutate("escape_status" = c("NA" = "not measured", "S" = "silenced / variable", "V" = "facultative", "E" = "constitutive")[.$`final status`]) 
# 
# escapee_annotation %>%
#   ggplot(aes(x = `nb of times escapee`, y = `nb of times silenced`, col = escape_status)) + geom_jitter() + 
#     theme_paper() + 
#     ggrepel::geom_text_repel(data = escapee_annotation[escapee_annotation$`final status` == "E", ], aes(label = symbol, col = escape_status)) + 
#     xlab("Number of studies in which gene escapes") + ylab("Number of studies in which gene is silenced") + 
#     scale_color_manual(values = escape_colors)
# 
# to_add_escaping <- ifelse(is.na(as.numeric(escapee_annotation$Yang2022_Status == "E")), 0, as.numeric(escapee_annotation$Yang2022_Status == "E")) + 
#   ifelse(is.na(as.numeric(escapee_annotation$Bowness22_Status == "E")), 0, as.numeric(escapee_annotation$Bowness22_Status == "E"))
# to_add_silenced <- ifelse(is.na(as.numeric(escapee_annotation$Yang2022_Status == "S")), 0, as.numeric(escapee_annotation$Yang2022_Status == "S")) + 
#   ifelse(is.na(as.numeric(escapee_annotation$Bowness22_Status == "S")), 0, as.numeric(escapee_annotation$Bowness22_Status == "S"))
# 
# escapee_annotation$`nb of times escapee 2` <- escapee_annotation$`nb of times escapee` + to_add_escaping
# escapee_annotation$`nb of times silenced 2` <- escapee_annotation$`nb of times silenced` + to_add_silenced
# 
# escapee_annotation <- escapee_annotation %>%
#   mutate(detected = `nb of times escapee 2` + `nb of times silenced 2`) %>%
#   mutate(ratio = `nb of times escapee 2` / detected) %>% 
#   mutate(final_status_2 = case_when(
#     ratio > .5 & detected > 3 ~ "Constitutive", 
#     (ratio <= .5 & ratio > 0) | (ratio > 0 & detected <= 3) ~ "Facultative", 
#     .default = "NPC-specific"
#   ))
# 
# escapee_annotation %>%
#   mutate(detected = `nb of times escapee 2` + `nb of times silenced 2`) %>%
#   mutate(ratio = `nb of times escapee 2` / detected) %>% 
#   {
#     ggplot(data = ., aes(x = detected, y = ratio, fill = final_status_2)) + 
#       geom_jitter(pch = 21, size = 3, width = .1, height = .005) + 
#       theme_paper() + 
#       ggrepel::geom_text_repel(data = .[.$`final status` == "E", ], aes(label = symbol, col = final_status_2)) + 
#       xlab("Number of studies in which gene is analyzed") + ylab("Fraction of studies in which gene escapes") + 
#       scale_fill_manual(values = escape_colors_other_name_2) + 
#       scale_color_manual(values = escape_colors_other_name_2) + 
#       annotate(geom = "rect", fill = "orange", xmin = 3, xmax = 16.5, ymin = 0.5, ymax = 1.05, alpha = .2) + 
#       annotate(geom = "rect", fill = "darkgreen", xmin = 0, xmax = 16.5, ymin = 0, ymax = 0.5, alpha = .2) + 
#       annotate(geom = "rect", fill = "darkgreen", xmin = 0, xmax = 3, ymin = 0.5, ymax = 1.05, alpha = .2)
#   }
# 
# escapee_annotation %>% dplyr::filter(detected > 0) %>% pull(final_status_2) %>% table()
# 
# ### Add number of days washout as variable
# dataset_complete$ndWashout <- unlist(lapply(dataset_complete$timeWO, function(x){
#   if (is.na(x) | x == "NO" | x == "NA" | x == "0"){
#     return(0)
#   } else if(x == "7") {
#     return(7)
#   } else {
#     nd_days <- as.numeric(stringr::str_split(x, ":")[[1]][[2]]) - as.numeric(stringr::str_split(x, ":")[[1]][[1]])
#   }
# }))
# 
# ### create condition annotations
# dataset_complete$Condition <- paste0("Aux_", dataset_complete$ndAux, "_Dox_", dataset_complete$ndDox, 
#                                      "_WO_", dataset_complete$ndWashout, "_WOAux", dataset_complete$WOwithAux)
# 
# ### Subset on main dataset for paper, remove May2022 (other clones) and crispr experiment (not necessary)
# data <- dataset_complete[,!dataset_complete$Guide %in% c("g13", "P3") & dataset_complete$Experiment != "May2022"]
# 
# name_conversion_condition <- list(
#   "Aux_0_Dox_0_WO_0_WOAuxNO" = "Control", 
#   
#   "Aux_0_Dox_3_WO_0_WOAuxNO" = "Dox (3d)", 
#   "Aux_0_Dox_3_WO_4_WOAuxNO" = "Dox (3d) - washout (4d)", 
#   
#   "Aux_0_Dox_7_WO_0_WOAuxNO" = "Dox (7d)", 
#   "Aux_0_Dox_7_WO_4_WOAuxNO" = "Dox (7d) - washout (4d)", 
#   "Aux_0_Dox_7_WO_7_WOAuxNO" = "Dox (7d) - washout (7d)",
#   "Aux_0_Dox_7_WO_14_WOAuxNO" = "Dox (7d) - washout (14d)",
#   "Aux_0_Dox_7_WO_4_WOAuxYES" = "Dox (7d) - washout (with Aux)",
#   
#   "Aux_0_Dox_14_WO_0_WOAuxNO" = "Dox (14d)",
#   "Aux_0_Dox_14_WO_7_WOAuxNO" = "Dox (14d) - washout (7d)",
#   "Aux_0_Dox_14_WO_14_WOAuxNO" = "Dox (14d) - washout (14d)",
#   "Aux_0_Dox_14_WO_21_WOAuxNO" = "Dox (14d) - washout (21d)",
#   
#   "Aux_0_Dox_21_WO_0_WOAuxNO" = "Dox (21d)",
#   "Aux_0_Dox_21_WO_7_WOAuxNO" = "Dox (21d) - washout (7d)",
#   "Aux_0_Dox_21_WO_14_WOAuxNO" = "Dox (21d) - washout (14d)",
#   "Aux_0_Dox_21_WO_21_WOAuxNO" = "Dox (21d) - washout (21d)",
#   
#   "Aux_2_Dox_0_WO_0_WOAuxNO" = "Aux (2d)", 
#   "Aux_2_Dox_0_WO_4_WOAuxNO" = "Aux (2d) - washout (4d)", 
#   "Aux_5_Dox_0_WO_0_WOAuxNO" = "Aux (5d)", 
#   "Aux_5_Dox_3_WO_0_WOAuxNO" = "Dox (3d), Aux (5d)", 
#   "Aux_5_Dox_3_WO_4_WOAuxNO" = "Dox (3d) - washout (4d), Aux (5d)", 
#   "Aux_9_Dox_0_WO_0_WOAuxNO" = "Aux (9d)", 
#   "Aux_9_Dox_7_WO_0_WOAuxNO" = "Dox (7d), Aux (9d)", 
#   "Aux_9_Dox_7_WO_4_WOAuxNO" = "Dox (7d) - washout (4d), Aux (9d)", 
#   "Aux_9_Dox_7_WO_4_WOAuxYES" = "Dox (7d) - washout (4d) (with Aux), Aux (9d)",
#   "Aux_21_Dox_21_WO_0_WOAuxNO" = "Dox (21d), Aux (21d)"
# )
# 
# # We use this grouping to subset the data
# data$ConditionClean <- unlist(name_conversion_condition[data$Condition])
# 
# #
# colData(data) <- colData(data)[ , c(
#   "SampleName", "CorrectedSampleName", "Clone", "Allele", "ndTreatment", "ndAux", "timeAux", "ndDox", "timeDox", "washout", "timeWO", "WOwithAux",
#   "Guide", "Experiment", "ndWashout", "Condition", "ConditionClean", "Replicate"
# )]
# 
# ###

saveRDS(sce_dataset, "./ProcessedData/merged_dataset_multiclone.rds")
