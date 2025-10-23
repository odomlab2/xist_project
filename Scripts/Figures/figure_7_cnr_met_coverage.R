library(GenomicRanges)
library(Rsamtools)
library(bamsignals)

setwd(paste0(dirname(rstudioapi::getSourceEditorContext()$path)))

blacklist <- rtracklayer::import("../../../ProcessedData/mm10_blacklist_280617.bed")

# get files
data_dir <- "../processed_files_antonia/"
all_files <- list.files(data_dir, pattern = "sorted_rmdup.bam")

plan_Gall <- data.frame(all_files) %>%
  mutate(Allele = str_extract(all_files, "(C57BL-6J)|(CAST-EiJ)")) %>% 
  mutate(Allele = ifelse(Allele == "C57BL-6J", "Xi_count", "Xa_count")) %>%
  mutate(sName = unlist(lapply(str_split(all_files, "_"), function(x){x[[1]]}))) %>%
  mutate(mark = unlist(lapply(str_split(all_files, "_"), function(x){x[[2]]}))) %>%
  pivot_wider(names_from = Allele, values_from = all_files) %>% 
  mutate(Xi_count = paste0("../processed_files_antonia/", Xi_count)) %>%
  mutate(Xa_count = paste0("../processed_files_antonia/", Xa_count))

pe.bam.all <- apply(plan_Gall[, c("Xi_count", "Xa_count")], 1, function(row) {
  paste(as.character(row))
})

ref.chr_all <- paste0("chr", c(1:19, "X", "Y"))

condition_rename = c(
  "NoDox" = "No Dox",
  "7dDox" = "7d Dox",
  "7dD7dWO" = "7d Dox - Washout",
  "21dDox" = "21d Dox",
  "21dD7dWO" = "21d Dox - Washout"
)

# 
get_interval_counts <- function(file, chr, start, end, n_windows = 20){
  
  bfl <- BamFile(file)
  param2 <- ScanBamParam(which=GRanges(seqnames = paste0(chr, ":", start, "-", end)))
  pparam <- PileupParam(distinguish_nucleotides=FALSE, distinguish_strands=FALSE)
  coverage <- pileup(bfl, scanBamParam=param2, pileupParam=pparam) %>% dplyr::select(-which_label) %>% 
    full_join(data.frame(pos = start:end, count = 0)) %>%
    arrange(pos) %>%
    mutate(window = cut(pos, n_windows, labels = F)) %>%
    group_by(window) %>%
    summarize(count = mean(count)) %>%
    add_column(chr = chr)
  
  coverage
}

get_counts_per_gene <- function(gene, sample_list, sample_to_file, gene_annotation, size_factors, n_windows = 20, offset = 1e4){
  
  chr = paste0("chr", seqnames(gene_annotation[gene_annotation$gene_name == gene]))
  start = start(gene_annotation[gene_annotation$gene_name == gene]) - offset
  end = end(gene_annotation[gene_annotation$gene_name == gene]) + offset
  
  print(end - start)
  
  all_counts <- lapply(sample_list, function(sample){
    
    file = sample_to_file[[sample]]
    sf = size_factors[[sample]]
    
    get_interval_counts(file, chr, start, end, n_windows) %>%
      add_column("Gene" = gene) %>%
      add_column("Sample" = sample) %>%
      mutate(count_norm = count / sf)
  })
  
  do.call("rbind", all_counts)
}

library(EnsDb.Mmusculus.v79)
gene_info <- genes(EnsDb.Mmusculus.v79)
gene_info <- gene_info[!duplicated(gene_info$symbol), ]
basic_chr <- c(1:19, "X", "Y")
gene_info <- gene_info[seqnames(gene_info) %in% basic_chr, ]

tss_ranges <- GRanges(
  seqnames = seqnames(gene_info),
  ranges = IRanges(
    start = ifelse(strand(gene_info) == "+", start(gene_info) - 500, end(gene_info) - 1000),
    end = ifelse(strand(gene_info) == "+", start(gene_info) + 1000, end(gene_info) + 500),
    names = names(gene_info)
  ),
  strand = strand(gene_info)
)
tss_ranges$gene_name = gene_info$gene_name

range_dist = 5e4
adjacent_ranges_left <- GRanges(
  seqnames = seqnames(gene_info),
  ranges = IRanges(
    start = ifelse(strand(gene_info) == "+", start(gene_info) - 500 - range_dist, start(gene_info) - range_dist),
    end = ifelse(strand(gene_info) == "+", start(gene_info) - 500, start(gene_info)),
    names = names(gene_info)
  ),
  strand = strand(gene_info)
)
adjacent_ranges_left$gene_name = gene_info$gene_name

adjacent_ranges_right <- GRanges(
  seqnames = seqnames(gene_info),
  ranges = IRanges(
    start = ifelse(strand(gene_info) == "+", end(gene_info), end(gene_info) + 500),
    end = ifelse(strand(gene_info) == "+", end(gene_info) + range_dist, end(gene_info) + 500 +  range_dist),
    names = names(gene_info)
  ),
  strand = strand(gene_info)
)
adjacent_ranges_right$gene_name = gene_info$gene_name

gene_info_x <- gene_info[seqnames(gene_info) == "X", ]

# 
swapping_order <- c(1, 2, 13, 4, 5, 6, 17, 8, 9, 10, 11, 12, 3, 14, 15, 16, 7, 18, 19, 20)
# swapping_order <- 1:20
escape_results <- read_csv("../../../ProcessedData/excel_files/figure3_d_score_differences.csv") %>% dplyr::select(-c("...1"))
escape_results %>% 
  mutate(gene_category = case_when(
    - `Dox (21 days), Washout (0 days)_dif` < .2 & - `Dox (21 days), Washout (7 days)_dif` < .2 ~ "Unresponsive", 
    - `Dox (21 days), Washout (0 days)_dif` > .2 & - `Dox (21 days), Washout (7 days)_dif` > .2 ~ "Responsive, Irreversible",
    - `Dox (21 days), Washout (0 days)_dif` > .2 & - `Dox (21 days), Washout (7 days)_dif` < .2 ~ "Responsive, Reversible", 
    .default = "test")) %>%
  dplyr::rename(Symbol = name) -> gene_annotation_reversibility

cnr_metadata <- readRDS("../../../ProcessedData/cnr_metadata.rds") %>%
  mutate(file_name_b6 = paste0("../processed_files_antonia/", "E6", Condition, "_", HistoneMark, "_", "C57BL-6J", "_sorted_rmdup.bam")[swapping_order]) %>%
  mutate(file_name_cast = paste0("../processed_files_antonia/", "E6", Condition, "_", HistoneMark, "_", "CAST-EiJ", "_sorted_rmdup.bam")[swapping_order]) %>%
  mutate(Condition = gsub("r[2|3]", "", Condition)) %>%
  mutate(Condition = ifelse(Condition == "21dD", "21dDox", Condition))
cnr_metadata <- readRDS("../../../ProcessedData/cnr_metadata.rds") %>%
  mutate(file_name_b6 = paste0("../processed_files_antonia/", "E6", Condition, "_", HistoneMark, "_", "C57BL-6J", "_sorted_rmdup.bam")[swapping_order]) %>%
  mutate(file_name_cast = paste0("../processed_files_antonia/", "E6", Condition, "_", HistoneMark, "_", "CAST-EiJ", "_sorted_rmdup.bam")[swapping_order]) %>%
  mutate(Condition = gsub("r[2|3]", "", Condition)) %>%
  mutate(Condition = ifelse(Condition == "21dD", "21dDox", Condition))
rownames(cnr_metadata) <- rownames(cnr_metadata)[swapping_order]
  
sample_to_file_xa <- cnr_metadata %>% rownames_to_column("SampleName") %>% pull(file_name_cast, name = SampleName)
sample_to_file_xi <- cnr_metadata %>% rownames_to_column("SampleName") %>% pull(file_name_b6, name = SampleName)

sample_to_sf <- cnr_metadata %>% rownames_to_column("SampleName") %>% pull(sizeFactor, name = SampleName)

# 

samples_here <- cnr_metadata %>% dplyr::filter(HistoneMark == "H3K27me3") %>% rownames_to_column("SampleName") %>% pull(SampleName)

gene_counts_test <- get_counts_per_gene("Kdm6a", samples_here, sample_to_file_xi, gene_info_x, sample_to_sf, n_windows = 100, offset = 1e5)

gene_counts_test %>%
  mutate(Condition = cnr_metadata[Sample, ]$Condition) %>%
  mutate(Condition = factor(condition_rename[Condition], levels = as.character(condition_rename))) %>%
  ggplot(aes(x = window, y = count_norm, group = Sample)) + geom_ribbon(aes(ymax = count_norm, ymin = 0), alpha = 0.3, fill = "brown") +
    facet_wrap(~Condition, ncol = 1) + theme_bw(base_size = 17) + xlab("") + ylab("Normalized Coverage")

samples_here_here <- c("E6NoDoxr2_H3K27me3", "E6NoDoxr3_H3K27me3")
me3_data_shelfs <- do.call("rbind", lapply(gene_annotation_reversibility_cut$Symbol, function(gene){
  
  strand = as.character(strand(gene_info_x[gene_info_x$gene_name == gene]))
  
  data_here_left_shelf <- get_counts_per_gene(gene, samples_here_here, sample_to_file_xi, adjacent_ranges_left, sample_to_sf, n_windows = 100, offset = 0) %>%
    add_column("Context" = "Left_Shelf")
  data_here_promoter = get_counts_per_gene(gene, samples_here_here, sample_to_file_xi, tss_ranges, sample_to_sf, n_windows = 10, offset = 0) %>%
    add_column("Context" = "Promoter") %>%
    mutate(window = window + max(data_here_left_shelf$window))
  data_here_gb = get_counts_per_gene(gene, samples_here_here, sample_to_file_xi, gene_info_x, sample_to_sf, n_windows = 100, offset = 0) %>%
    add_column("Context" = "Genebody") %>%
    mutate(window = window + max(data_here_promoter$window))
  data_here_right_shelf <- get_counts_per_gene(gene, samples_here_here, sample_to_file_xi, adjacent_ranges_left, sample_to_sf, n_windows = 100, offset = 0) %>%
    add_column("Context" = "Left_Shelf") %>%
    mutate(window = window + max(data_here_gb$window))
  
  # if negative strand, invert genebody results and promoter
  if (strand == "-"){
    data_here_gb %>% 
      group_by(Gene, Sample) %>%
      mutate(count_norm = rev(count_norm)) -> data_here_gb
    data_here_promoter %>% 
      group_by(Gene, Sample) %>%
      mutate(count_norm = rev(count_norm)) -> data_here_promoter
    data_here_left_shelf %>% 
      group_by(Gene, Sample) %>%
      mutate(count_norm = rev(count_norm)) -> data_here_left_shelf
    data_here_right_shelf %>% 
      group_by(Gene, Sample) %>%
      mutate(count_norm = rev(count_norm)) -> data_here_right_shelf
  }
  
  data_done <- rbind(data_here_promoter, data_here_gb, data_here_left_shelf, data_here_right_shelf)
  data_done
  
}))

me3_data_shelfs %>%
  mutate(Condition = cnr_metadata[Sample, ]$Condition) %>%
  group_by(window, Gene, Sample) %>%
  summarize(count_norm = mean(count_norm)) %>%
  ggplot(aes(x = window, y = count_norm, group = Sample)) + geom_ribbon(aes(ymax = count_norm, ymin = 0), alpha = 0.3, fill = "brown") +
    theme_bw(base_size = 17) + xlab("") + ylab("Normalized Coverage")

## get counts for all escapees and make metaplot 
gene_annotation_reversibility_cut <- gene_annotation_reversibility %>% 
  dplyr::filter(gene_category %in% c("Responsive, Irreversible", "Responsive, Reversible"))

# me3_data <- do.call("rbind", lapply(gene_annotation_reversibility_cut$Symbol, function(gene){
#   
#   strand = as.character(strand(gene_info_x[gene_info_x$gene_name == gene]))
#   
#   data_here_promoter = get_counts_per_gene(gene, samples_here, sample_to_file_xi, tss_ranges, sample_to_sf, n_windows = 10, offset = 0) %>%
#     add_column("Context" = "Promoter")
#   data_here_gb = get_counts_per_gene(gene, samples_here, sample_to_file_xi, gene_info_x, sample_to_sf, n_windows = 100, offset = 0) %>%
#     add_column("Context" = "Genebody") %>%
#     mutate(window = window + max(data_here_promoter$window))
#   
#   if (strand == "-"){
#     data_here_gb %>% 
#       group_by(Gene, Sample) %>%
#       mutate(count_norm = rev(count_norm)) -> data_here_gb
#     data_here_promoter %>% 
#       group_by(Gene, Sample) %>%
#       mutate(count_norm = rev(count_norm)) -> data_here_promoter
#   }
#   
#   data_done <- rbind(data_here_promoter, data_here_gb) %>%
#     add_column("Reversibility" = gene_annotation_reversibility_cut[gene_annotation_reversibility_cut$Symbol == gene, ]$gene_category)
#   data_done
#   
# }))
# 
# saveRDS(me3_data, "../me3_gene_coverage.rds")

me3_data <- readRDS("../../../ProcessedData/me3_gene_coverage.rds")

me3_data %>%
  mutate(Condition = cnr_metadata[Sample, ]$Condition) %>%
  group_by(Gene, Condition, Reversibility) %>%
  summarize(count_norm = mean(count_norm)) %>%
  mutate(Reversibility = factor(Reversibility, levels = c("Responsive, Reversible", "Responsive, Irreversible"))) %>%
  mutate(Condition = factor(Condition, levels = c("NoDox", "7dDox", "7dD7dWO", "21dDox", "21dD7dWO"))) %>%
  ggplot(aes(x = Condition, y = count_norm, col = Reversibility)) + stat_summary(size = 1.5, position = position_dodge(.3)) + 
    scale_color_manual(values = c("darkgreen", "darkgrey"))

me3_data %>%
  mutate(Condition = cnr_metadata[Sample, ]$Condition) %>%
  mutate(Replicate = cnr_metadata[Sample, ]$Replicate) %>%
  group_by(Gene, Condition, Reversibility, Replicate) %>%
  summarize(count_norm = mean(count_norm)) %>%
  group_by(Gene, Reversibility, Replicate) %>%
  mutate(lfc = log2((count_norm + 0.1) / (mean(count_norm[Condition == "NoDox"]) + 0.1))) %>%
  group_by(Gene, Condition, Reversibility) %>%
  summarize(lfc = mean(lfc)) %>%
  mutate(Reversibility = factor(Reversibility, levels = c("Responsive, Reversible", "Responsive, Irreversible"))) %>%
  mutate(Condition = factor(Condition, levels = c("NoDox", "7dDox", "7dD7dWO", "21dDox", "21dD7dWO"))) %>%
  ggplot(aes(x = Condition, y = lfc, col = Reversibility)) + stat_summary(size = 1.5, position = position_dodge(.3)) + 
    scale_color_manual(values = c("darkgreen", "darkgrey"))
#ggsave("../plots/me3_coverage_foldchanges_summary.pdf", width = 10, height = 10)

me3_data %>%
  dplyr::select(-c("count", "chr")) %>%
  mutate(Condition = cnr_metadata[Sample, ]$Condition) %>%
  mutate(Replicate = cnr_metadata[Sample, ]$Replicate) %>%
  group_by(window, Gene, Context, Reversibility, Replicate) %>%
  mutate(lfc = log2((count_norm + 0.1) / (mean(count_norm[Condition == "NoDox"]) + 0.1))) %>%
  ungroup() %>%
  group_by(window, Context, Reversibility, Condition) %>%
  summarize(lfc = mean(lfc)) %>%
  mutate(Condition = factor(condition_rename[Condition], levels = as.character(condition_rename))) -> df_plot_7e_1

df_plot_7e_1 %>%
  ggplot(aes(x = window, y = lfc, fill = Reversibility)) + geom_ribbon(aes(ymax = lfc, ymin = 0), alpha = 0.5) +
    facet_wrap(~Condition, ncol = 1) + theme_bw(base_size = 17) + xlab("Genomic Coordinate") + ylab("Coverage (Log Fold Change)") + 
    geom_vline(xintercept = 10) + scale_fill_manual(values = c("darkgrey", "darkgreen")) + ylim(c(-0.1, 3.5))
#ggsave("../plots/me3_coverage_foldchanges.pdf", width = 10, height = 10)

### same for ubi
samples_here <- cnr_metadata %>% dplyr::filter(HistoneMark == "H2AK119ubi") %>% rownames_to_column("SampleName") %>% pull(SampleName)

# ubi_data <- do.call("rbind", lapply(gene_annotation_reversibility_cut$Symbol, function(gene){
#   
#   strand = as.character(strand(gene_info_x[gene_info_x$gene_name == gene]))
#   
#   data_here_promoter = get_counts_per_gene(gene, samples_here, sample_to_file_xi, tss_ranges, sample_to_sf, n_windows = 10, offset = 0) %>%
#     add_column("Context" = "Promoter")
#   data_here_gb = get_counts_per_gene(gene, samples_here, sample_to_file_xi, gene_info_x, sample_to_sf, n_windows = 100, offset = 0) %>%
#     add_column("Context" = "Genebody") %>%
#     mutate(window = window + max(data_here_promoter$window))
#   
#   if (strand == "-"){
#     data_here_gb %>% 
#       group_by(Gene, Sample) %>%
#       mutate(count_norm = rev(count_norm)) -> data_here_gb
#     data_here_promoter %>% 
#       group_by(Gene, Sample) %>%
#       mutate(count_norm = rev(count_norm)) -> data_here_promoter
#   }
#   
#   data_done <- rbind(data_here_promoter, data_here_gb) %>%
#     add_column("Reversibility" = gene_annotation_reversibility_cut[gene_annotation_reversibility_cut$Symbol == gene, ]$gene_category)
#   data_done
# }))
# 
# saveRDS(ubi_data, "../ubi_gene_coverage.rds")

ubi_data <- readRDS("../../../ProcessedData/ubi_gene_coverage.rds")

ubi_data %>%
  mutate(Condition = cnr_metadata[Sample, ]$Condition) %>%
  group_by(Gene, Condition, Reversibility) %>%
  summarize(count_norm = mean(count_norm)) %>%
  mutate(Reversibility = factor(Reversibility, levels = c("Responsive, Reversible", "Responsive, Irreversible"))) %>%
  mutate(Condition = factor(Condition, levels = c("NoDox", "7dDox", "7dD7dWO", "21dDox", "21dD7dWO"))) %>%
  ggplot(aes(x = Condition, y = count_norm, col = Reversibility)) + stat_summary(size = 1.5, position = position_dodge(.3)) + 
    scale_color_manual(values = c("darkgreen", "darkgrey"))

ubi_data %>%
  mutate(Condition = cnr_metadata[Sample, ]$Condition) %>%
  mutate(Replicate = cnr_metadata[Sample, ]$Replicate) %>%
  group_by(Gene, Condition, Reversibility, Replicate) %>%
  summarize(count_norm = mean(count_norm)) %>%
  group_by(Gene, Reversibility, Replicate) %>%
  mutate(lfc = log2((count_norm + 0.1) / (mean(count_norm[Condition == "NoDox"]) + 0.1))) %>%
  group_by(Gene, Condition, Reversibility) %>%
  summarize(lfc = mean(lfc)) %>%
  mutate(Reversibility = factor(Reversibility, levels = c("Responsive, Reversible", "Responsive, Irreversible"))) %>%
  mutate(Condition = factor(Condition, levels = c("NoDox", "7dDox", "7dD7dWO", "21dDox", "21dD7dWO"))) %>%
  ggplot(aes(x = Condition, y = lfc, col = Reversibility)) + stat_summary(size = 1.5, position = position_dodge(.3)) + 
    scale_color_manual(values = c("darkgreen", "darkgrey"))
#ggsave("../plots/ubi_coverage_foldchanges_summary.pdf", width = 10, height = 10)

ubi_data %>%
  dplyr::select(-c("count", "chr")) %>%
  mutate(Condition = cnr_metadata[Sample, ]$Condition) %>%
  mutate(Replicate = cnr_metadata[Sample, ]$Replicate) %>%
  group_by(window, Gene, Context, Reversibility, Replicate) %>%
  mutate(lfc = log2((count_norm + 0.1) / (mean(count_norm[Condition == "NoDox"]) + 0.1))) %>%
  ungroup() %>%
  group_by(window, Context, Reversibility, Condition) %>%
  summarize(lfc = mean(lfc)) %>%
  mutate(Condition = factor(condition_rename[Condition], levels = as.character(condition_rename))) -> df_plot_7e_2

df_plot_7e_2 %>%
  ggplot(aes(x = window, y = lfc, fill = Reversibility)) + geom_ribbon(aes(ymax = lfc, ymin = 0), alpha = 0.5) +
    facet_wrap(~Condition, ncol = 1) + theme_bw(base_size = 17) + xlab("Genomic Coordinate") + ylab("Coverage (Log Fold Change)") + 
    geom_vline(xintercept = 10) + scale_fill_manual(values = c("darkgrey", "darkgreen")) + ylim(c(-0.1, 3.5))
# ggsave("../plots/ubi_coverage_foldchanges.pdf", width = 10, height = 10)

rbind(cbind(df_plot_7e_1, "HistoneMark" = "meth"), 
      cbind(df_plot_7e_2, "HistoneMark" = "ubi")) %>% list() %>% saveRDS(., "../../../ProcessedData/source_data/source_data_5_coverage.rds")
