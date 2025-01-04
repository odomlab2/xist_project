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
  Sex = c(rep("female", 4), "male")
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

###

data <- sce_dataset1

data_before_filtering <-  data[seqnames(data) == "X", ]

data <- data[(rowSums(counts_inactive(data)) + rowSums(counts_active(data))) / ncol(data) > 10, ]

data_with_autosomes <- data
data <- data[seqnames(data) == "X", ]

test <- counts(sce_dataset1[seqnames(rowRanges(sce_dataset1)) == "Y", ])

pdf("./Plots/Fig5/y_expression.pdf")
pheatmap(log10(test[rowSums(test) > 0, ] + 1), cluster_rows = F, cluster_cols = F)
dev.off()


# look at X/Y coverage
data.frame(
  Ycov = colSums(counts(sce_dataset1[seqnames(rowRanges(sce_dataset1)) == "Y", ])) / colSums(counts(sce_dataset1)), 
  Xcov = colSums(counts(sce_dataset1[seqnames(rowRanges(sce_dataset1)) == "X", ])) / colSums(counts(sce_dataset1))
) %>%
  ggplot(aes(Xcov, Ycov)) + geom_point()


sample_show = 1
colnames(data_before_filtering)[[sample_show]]

data.frame(
  counts_total = counts(data_before_filtering[,sample_show])[[1]], 
  counts_allelic = (counts_active(data_before_filtering[,sample_show]) +  counts_inactive(data_before_filtering[,sample_show]))[[1]]
) %>%
  ggplot(aes(x = counts_total + 1, y = counts_allelic + 1)) + geom_point() + scale_x_log10() + scale_y_log10() + 
  theme_paper() + geom_vline(xintercept = 10, linetype = "dashed") + xlab("Total reads mapped") + ylab("Allele-specific reads mapped") + 
  geom_abline(linetype = "dashed")
ggsave("./Plots/Fig5/mapped_reads_sample1.pdf")

data.frame(
  chromosome = as.character(seqnames(rowRanges(data_with_autosomes))) == "X", 
  total = (counts_active(data_with_autosomes[,sample_show]) + counts_inactive(data_with_autosomes[,sample_show]))[[1]], 
  ase =  (counts_active(data_with_autosomes[,sample_show]) / (counts_active(data_with_autosomes[,sample_show]) + counts_inactive(data_with_autosomes[,sample_show])))[[1]]
) %>%
  ggplot(aes(x = ase, fill = chromosome)) + geom_density() + theme_paper() + ylab("Density") + xlab("B6 / (B6 + JF1)")
ggsave("./Plots/Fig5/density_ase_c30.pdf")


full_df <- do.call("rbind", lapply(1:ncol(data), function(i){
  sample_show = i
  data.frame(
    chromosome = as.character(seqnames(rowRanges(data_with_autosomes))) == "X", 
    total = (counts_active(data_with_autosomes[,sample_show]) + counts_inactive(data_with_autosomes[,sample_show]))[[1]], 
    ase = (counts_active(data_with_autosomes[,sample_show]) / (counts_active(data_with_autosomes[,sample_show]) + counts_inactive(data_with_autosomes[,sample_show])))[[1]], 
    sample = data$Sample[[i]]
  )
}))

full_df %>% 
  mutate(sample = factor(sample, levels = paste0("Sample", c(1, 8, 10, 11, 12)))) %>%
  ggplot(aes(x = ase, fill = chromosome)) + geom_density() + theme_paper() + ylab("Density") + xlab("B6 / (B6 + JF1)") + facet_wrap(~sample)
ggsave("./Plots/Fig5/density_all_clones.pdf")


### 
data_here <- data
ratios <- counts_inactive(data_here) / (counts_inactive(data_here) + counts_active(data_here))
# ratios <- ratios[!rownames(ratios) %in% genes_out, ]

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

add_to_sce <- setNames(rep(NA, nrow(data)), rownames(data))
add_to_sce[escapee_annotation$symbol] <- escapee_annotation$our_status2

# rowData(data) <- cbind(rowData(data), "EscapeAnnotation" = add_to_sce[-length(add_to_sce)])

# For Xist-expression, we need a scaling factor per sample to account for differences in sequencing depth

sample_nr = 1
data.frame(
  total = rowSums(counts_inactive(data[,sample_nr])) + rowSums(counts_active(data[,sample_nr])),
  ase = rowSums(counts_inactive(data[,sample_nr])) / (rowSums(counts_inactive(data[,sample_nr])) +
                                                                   rowSums(counts_active(data[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle(paste0("Sample", sample_nr)) +
  xlab("Total expression (allelic reads)") + ylab("ASE (d-score)")
ggsave("./Plots/Fig5/sample1_dscore_plot.pdf")

sample_nr = 2
data.frame(
  total = rowSums(counts_inactive(data[,sample_nr])) + rowSums(counts_active(data[,sample_nr])),
  ase = rowSums(counts_inactive(data[,sample_nr])) / (rowSums(counts_inactive(data[,sample_nr])) +
                                                        rowSums(counts_active(data[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle(colnames(data)[[sample_nr]]) +
  xlab("Total expression (allelic reads)") + ylab("ASE (d-score)")
ggsave("./Plots/Fig5/sample8_dscore_plot.pdf")

sample_nr = 3
data.frame(
  total = rowSums(counts_inactive(data[,sample_nr])) + rowSums(counts_active(data[,sample_nr])),
  ase = rowSums(counts_inactive(data[,sample_nr])) / (rowSums(counts_inactive(data[,sample_nr])) +
                                                        rowSums(counts_active(data[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle(colnames(data)[[sample_nr]]) +
  xlab("Total expression (allelic reads)") + ylab("ASE (d-score)")
ggsave("./Plots/Fig5/sample10_dscore_plot.pdf")

sample_nr = 4
data.frame(
  total = rowSums(counts_inactive(data[,sample_nr])) + rowSums(counts_active(data[,sample_nr])),
  ase = rowSums(counts_inactive(data[,sample_nr])) / (rowSums(counts_inactive(data[,sample_nr])) +
                                                        rowSums(counts_active(data[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle(colnames(data)[[sample_nr]]) +
  xlab("Total expression (allelic reads)") + ylab("ASE (d-score)")
ggsave("./Plots/Fig5/sample11_dscore_plot.pdf")

sample_nr = 5
data.frame(
  total = rowSums(counts_inactive(data[,sample_nr])) + rowSums(counts_active(data[,sample_nr])),
  ase = rowSums(counts_inactive(data[,sample_nr])) / (rowSums(counts_inactive(data[,sample_nr])) +
                                                        rowSums(counts_active(data[,sample_nr])))
) %>%
  rownames_to_column("Gene") %>%
  ggplot(aes(total, ase)) + geom_point() + ggrepel::geom_text_repel(aes(label = Gene)) +
  theme_paper() + scale_x_log10() + geom_hline(yintercept = 0.1, linetype = "dashed") + ggtitle(colnames(data)[[sample_nr]]) +
  xlab("Total expression (allelic reads)") + ylab("ASE (d-score)")
ggsave("./Plots/Fig5/sample12_dscore_plot.pdf")



# check that Xist-overexpression works by plotting Xist counts per million (CPM) per sample
size_factors <- colSums(counts(data_with_autosomes)) / colSums(counts(data_with_autosomes))[[1]]

data.frame(
  Sample = factor(colnames(data), levels = paste0("Sample", c(1, 8, 10, 11, 12))), 
  Expression = as.numeric(counts(data_here["Xist", ]) / colSums(counts(data_here)) * 1e6)
) %>%
  ggplot(aes(x = Sample, y = Expression)) + 
    stat_summary(stat = "mean", geom = "bar", fill = "grey", col = "black") + 
    coord_flip() + theme_paper() + 
    xlab("") + ylab("Expression (CPM)") + geom_jitter(size = 3, width = 0.1) + 
    # scale_y_continuous(expand = c(0, 0), limits = c(0, 320000)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    ggtitle("Xist expression")
ggsave("./Plots/Fig5/xist_expression.pdf")

# Distribution d-scores
ratios %>%
  t() %>% data.frame() %>%
  add_column(Sample = factor(data$Sample, levels = paste0("Sample", c(1, 8, 10, 11, 12)))) %>%
  pivot_longer(-c("Sample")) %>%
  # dplyr::filter(!is.na(escape_status)) %>%
  ggplot(aes(x = Sample, y = value)) + 
    geom_boxplot(col = "black", outlier.color = "grey") + 
    theme_paper() + ylab("ASE (d-score)") + 
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) + 
    # scale_fill_manual(values = c("grey", "darkgreen", "orange"), labels = c("Silenced / Variable", "Facultative", "Constitutive")) + 
    labs(fill = "Escape Category") + xlab("")
ggsave("./Plots/Fig5/d_score_distribution.pdf")



