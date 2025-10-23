### format source data
# S4Vectors::rename

library(openxlsx)

setwd(paste0(dirname(rstudioapi::getSourceEditorContext()$path), "/../../../"))

#
data <- readRDS("./ProcessedData/merged_dataset_paper.rds")
colData(data) %>% data.frame() %>% dplyr::select(c("Experiment", "Replicate")) %>% unique() %>% pull(Replicate, name = Experiment) -> exp_to_clone

### read prepared data
source_data_0 <- readRDS("./ProcessedData/source_data/source_data_0.rds")
source_data_1 <- readRDS("./ProcessedData/source_data/source_data_1.rds")
source_data_2 <- readRDS("./ProcessedData/source_data/source_data_2.rds")
source_data_3 <- readRDS("./ProcessedData/source_data/source_data_3.rds")
source_data_3_deg <- readRDS("./ProcessedData/source_data/source_data_3_deg.rds")
source_data_4 <- readRDS("./ProcessedData/source_data/source_data_4.rds")
source_data_5 <- readRDS("./ProcessedData/source_data/source_data_5.rds")
source_data_5_add <- readRDS("./ProcessedData/source_data/source_data_5_add.rds")
source_data_5_coverage <- readRDS("./ProcessedData/source_data/source_data_5_coverage.rds")
source_data_5_met <- readRDS("./ProcessedData/source_data/source_data_5_met.rds")

### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
# Figure 1
### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
file_export_name <- "./ProcessedData/source_data/source_data_fig_1.xlsx"

export_fig1a <- source_data_0[[1]]

export_fig1d <- source_data_1[[1]] %>% 
  mutate(Clone = gsub("E6 - ", "", Clone)) %>% 
  mutate(Replicate = exp_to_clone[Clone]) %>%
  dplyr::select(c("Condition", "Replicate", "Expression")) %>% 
  rename(Expression = "`Xist Expression`") %>% 
  mutate(Condition = factor(Condition, levels = c("Control", "Dox (3 days)", "Dox (7 days)", "Dox (14 days)", "Dox (21 days)"))) %>%
  arrange(Condition)

export_fig1e <- source_data_1[[2]] %>% 
  mutate(Clone = gsub("E6 - ", "", Clone)) %>% 
  mutate(Replicate = exp_to_clone[Clone]) %>%
  rename(name = "Gene", value = "Allelic Ratio") %>%
  dplyr::select(c("Condition", "Replicate", "Gene", "Allelic Ratio", "position", "Clusterings"))

export_fig1f <- source_data_1[[3]] %>%
  mutate(Replicate = gsub("E6 - ", "", Clone)) %>%
  rename(name = "Gene", value = "Allelic Ratio") %>%
  dplyr::select(c("Condition", "Replicate", "Gene", "Allelic Ratio", "position", "Clusterings"))

export_fig1g <- source_data_1[[4]] %>%
  mutate(Replicate = gsub("E6 - ", "", Clone)) %>%
  rename(name = "Gene", value = "Allelic Ratio") %>%
  dplyr::select(c("Condition", "Replicate", "Gene", "Allelic Ratio", "escape_status")) %>%
  rename(escape_status = "Escape Status")

export_fig1i <- source_data_1[[5]] %>%
  dplyr::select(c("gene", "category", "halflifes_off", "b_off")) %>%
  rename(gene = "Gene", category = "Escape Status", halflifes_off = "Halflife fit (offset model)", b_off = "Offset fit (offset model)")

export_fig1j <- source_data_1[[6]] %>%
  rename(clustering = "Clusterings", halftime_mean = "Halflife fit mean (offset model)", sd_halftime = "Halflife fit standard dev (offset model)", 
         position = "Position", n_genes = "Number of genes", genes = "Genes")

export_fig1k <- source_data_1[[5]] %>%
  dplyr::select(c("gene", "category", "halflifes_off", "b_off")) %>%
  rename(gene = "Gene", category = "Escape Status", halflifes_off = "Halflife fit (offset model)", b_off = "Offset fit (offset model)")

list_of_datasets <- list("Figure 1a" = export_fig1a, 
                         "Figure 1d" = export_fig1d,
                         "Figure 1e" = export_fig1e,
                         "Figure 1f" = export_fig1f,
                         "Figure 1g" = export_fig1g,
                         "Figure 1i" = export_fig1i,
                         "Figure 1j" = export_fig1j,
                         "Figure 1k" = export_fig1k
)

openxlsx::write.xlsx(list_of_datasets, file = file_export_name)


### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
# Figure S1
### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###

file_export_name <- "./ProcessedData/source_data/source_data_fig_s1.xlsx"

export_figs1a <- source_data_0[[2]] %>%
  rename(counts_total = "Total Counts", counts_allelic = "Allele-specific Counts")

export_figs1b <- source_data_0[[3]] %>%
  mutate(chromosome = ifelse(chromosome, "X-chromosome", "Autosome")) %>%
  rename(chromosome = "Chromosome", total = "Total Counts", ase = "Allelic Ratio")

export_figs1c <- source_data_0[[4]] %>%
  rename(total = "Total Counts", ase = "Allelic Ratio")

export_figs1d <- source_data_0[[5]] %>%
  rename(name = "Gene", value = "Allelic Ratio", escape_in_clones = "Escapes in all clones") %>%
  dplyr::select(c("Gene", "Sample", "Allelic Ratio", "Escapes in all clones"))

export_figs1f <- source_data_1[[8]] %>%
  dplyr::select(c("name", "Clone", "Condition", "position", "value")) %>%
  rename(name = "Gene", value = "Allelic Ratio", position = "Position")

export_figs1g <- source_data_1[[9]] %>%
  dplyr::select(c("name", "Clone", "Condition", "escape_in_clones", "escape_group", "value")) %>%
  rename(name = "Gene", value = "Allelic Ratio", escape_in_clones = "Escapes in all clones", escape_group = "Escape Status")

export_figs1h <- readRDS("./ProcessedData/source_data/source_data_s1h.rds") %>%
  rename(ENSEMBL_v102 = "Ensembl ID", symbol = "Gene Symbol", detected = "Studies detected", ratio = "Fraction of studies in which gene escapes")

export_figs1i <- source_data_1[[10]] %>%
  dplyr::select(c("name", "Clone", "Condition", "escape_in_clones", "escape_group", "value")) %>%
  rename(name = "Gene", value = "Allelic Ratio", escape_in_clones = "Escapes in all clones", escape_group = "Escape Status")

list_of_datasets <- list("Figure S1a" = export_figs1a, 
                         "Figure S1b" = export_figs1b,
                         "Figure S1c" = export_figs1c,
                         "Figure S1d" = export_figs1d,
                         "Figure S1f" = export_figs1f,
                         "Figure S1g" = export_figs1g,
                         "Figure S1h" = export_figs1h,
                         "Figure S1i" = export_figs1i
)

openxlsx::write.xlsx(list_of_datasets, file = file_export_name)

### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
# Figure S2
### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###

file_export_name <- "./ProcessedData/source_data/source_data_fig_s2.xlsx"

export_figs2a <- source_data_1[[11]] %>%
  dplyr::select(c("name", "position", "Clone", "Condition", "value")) %>%
  rename(name = "Gene", position = "Position", "value" = "Allelic Ratio")

export_figs2b <- source_data_1[[12]] %>%
  dplyr::select(c("name", "Clone", "Condition", "value")) %>%
  rename(name = "Gene", "value" = "Allelic Ratio")

export_figs2c <- source_data_1[[13]] %>%
  dplyr::select(c("name", "position", "Clone", "Condition", "value")) %>%
  rename(name = "Gene", position = "Position", "value" = "Allelic Ratio")

export_figs2d <- source_data_1[[14]] %>%
  dplyr::select(c("name", "escape_status", "Clone", "Condition", "value")) %>%
  rename(name = "Gene", escape_status = "Escape Status", "value" = "Allelic Ratio")

export_figs2e <- source_data_1[[15]] %>%
  dplyr::select(c("Clone", "Condition", "Expression"))

export_figs2f <- source_data_1[[16]] %>%
  dplyr::select(c("name", "escape_status", "Clone", "Condition", "value")) %>%
  rename(name = "Gene", escape_status = "Escape Status", "value" = "Allelic Ratio")

export_figs2g <- source_data_1[[17]] %>%
  dplyr::select(c("name", "position", "Clone", "Condition", "value")) %>%
  rename(name = "Gene", position = "Position", "value" = "Allelic Ratio")

export_figs2h <- source_data_1[[18]] %>%
  dplyr::select(c("name", "position", "Clone", "Condition", "value")) %>%
  rename(name = "Gene", position = "Position", "value" = "Allelic Ratio")

export_figs2i <- source_data_1[[19]] %>%
  dplyr::select(c("gene", "EscapeAnnotation", "Clone", "Timepoint", "Control", "Treatment", "pvalue")) %>%
  rename(gene = "Gene", EscapeAnnotation = "Escape Status")

list_of_datasets <- list("Figure S2a" = export_figs2a, 
                         "Figure S2b" = export_figs2b,
                         "Figure S2c" = export_figs2c,
                         "Figure S2d" = export_figs2d,
                         "Figure S2e" = export_figs2e,
                         "Figure S2f" = export_figs2f,
                         "Figure S2g" = export_figs2g,
                         "Figure S2h" = export_figs2h,
                         "Figure S2i" = export_figs2i
)

openxlsx::write.xlsx(list_of_datasets, file = file_export_name)

### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
# Figure S3
### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###

file_export_name <- "./ProcessedData/source_data/source_data_fig_s3.xlsx"

export_figs3a <- source_data_1[[20]] %>%
  dplyr::select(c("gene", "category", "BIC_native", "BIC_offset")) %>%
  rename(gene = "Gene", category = "Escape Status")

export_figs3b <- source_data_1[[21]] %>%
  dplyr::select(c("gene", "experiment", "timepoint", "active", "inactive", "rates")) %>%
  mutate(experiment = exp_to_clone[experiment]) %>%
  rename(gene = "Gene", experiment = "Replicate", "timepoint" = "Timepoint", "active" = "Counts (active X)", "inactive" = "Counts (inactive X)", "rates" = "Allelic ratio")

export_figs3c <- source_data_1[[22]] %>%
  rownames_to_column("Gene") %>%
  dplyr::select(c("Gene", "EscapeCategory", "BasalExpressionLevel", "EscapeDifference")) %>%
  rename(EscapeCategory = "Escape Status", "BasalExpressionLevel" = "Expression level X total", "EscapeDifference" = "Difference in allelic ratio")

export_figs3d <- source_data_1[[23]] %>%
  rownames_to_column("Gene") %>%
  dplyr::select(c("Gene", "EscapeCategory", "BasalInactiveExpressionLevel", "EscapeDifference")) %>%
  rename(EscapeCategory = "Escape Status", "BasalInactiveExpressionLevel" = "Expression level X inactive", "EscapeDifference" = "Difference in allelic ratio")

export_figs3e <- source_data_1[[24]] %>%
  dplyr::select(c("condition", "halftime_diff")) %>%
  rename(condition = "Condition", halftime_diff = "Difference in halftime")

export_figs3f <- source_data_1[[25]] %>%
  dplyr::select(c("Gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"))

export_figs3g <- source_data_3[[1]] %>%
  dplyr::select(c("Gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"))

export_figs3h <- source_data_3[[2]] %>%
  dplyr::select(c("Gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"))

export_figs3i <- source_data_3_deg[[1]] %>%
  dplyr::select(c(condition, comparison, it, ratio)) %>%
  rename(condition = "Condition", comparison = "Comparison", it = "Iteration", ratio = "Ratio")

export_figs3j <- source_data_2[[1]] %>%
  dplyr::select(c("Gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"))

export_figs3k <- source_data_2[[2]] %>%
  dplyr::select(c(condition, it, ratio)) %>%
  rename(condition = "Condition", it = "Iteration", ratio = "Ratio")

export_figs3l <- source_data_3_deg[[2]] %>%
  rename(cor = "Correlation")

export_figs3m <- source_data_3_deg[[3]] %>%
  dplyr::select(c("Gene", "padj", "Comparison", "Direction"))


list_of_datasets <- list("Figure S3a" = export_figs3a, 
                         "Figure S3b" = export_figs3b,
                         "Figure S3c" = export_figs3c,
                         "Figure S3d" = export_figs3d,
                         "Figure S3e" = export_figs3e,
                         "Figure S3f" = export_figs3f,
                         "Figure S3g" = export_figs3g,
                         "Figure S3h" = export_figs3h,
                         "Figure S3i" = export_figs3i,
                         "Figure S3j" = export_figs3j,
                         "Figure S3k" = export_figs3k,
                         "Figure S3l" = export_figs3l,
                         "Figure S3m" = export_figs3m
)

openxlsx::write.xlsx(list_of_datasets, file = file_export_name)

### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
# Figure S4
### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###

file_export_name <- "./ProcessedData/source_data/source_data_fig_s4.xlsx"

export_figs4a <- source_data_1[[26]] %>%
  dplyr::select(c("Gene", "Status", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"))

export_figs4b <- source_data_1[[27]] %>%
  dplyr::select(c("Gene", "Status", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"))

export_figs4c <- source_data_1[[28]] %>%
  dplyr::select(c("Gene", "genomic_position", "Status", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")) %>%
  rename(genomic_position = "Position")

export_figs4d <- source_data_1[[29]] %>%
  dplyr::select(c("Gene", "genomic_position", "Status", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")) %>%
  rename(genomic_position = "Position")

export_figs4e <- source_data_2[[3]] %>%
  dplyr::select(c("Gene", "genomic_position", "Status", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")) %>%
  rename(genomic_position = "Position")

export_figs4f <- source_data_3[[3]] %>%
  dplyr::select(c("Gene", "genomic_position", "Status", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")) %>%
  rename(genomic_position = "Position")

export_figs4g <- source_data_3[[4]] %>%
  dplyr::select(c("Gene", "genomic_position", "Status", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")) %>%
  rename(genomic_position = "Position")


list_of_datasets <- list("Figure S4a" = export_figs4a, 
                         "Figure S4b" = export_figs4b,
                         "Figure S4c" = export_figs4c,
                         "Figure S4d" = export_figs4d,
                         "Figure S4e" = export_figs4e,
                         "Figure S4f" = export_figs4f,
                         "Figure S4g" = export_figs4g
)

openxlsx::write.xlsx(list_of_datasets, file = file_export_name)

### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
# Figure 2
### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###

file_export_name <- "./ProcessedData/source_data/source_data_fig_2.xlsx"

export_figs2c <- source_data_2[[4]]

export_figs2d <- source_data_2[[5]] %>% 
  dplyr::select(c("name", "position", "Condition", "Replicate", "value")) %>%
  rename(name = "Gene", position = "Position", value = "Allelic Ratio")

export_figs2e <- source_data_2[[6]] %>% 
  dplyr::select(c("name", "escape_status", "Condition...1", "Condition...6", "Astrocyte", "NPC")) %>%
  rename(name = "Gene", escape_status = "Escape Status", Condition...1 = "Condition x-axis", Condition...6 = "Condition y-axis", NPC = "Allelic Ratio NPC", Astrocyte = "Allelic Ratio Astrocyte")

list_of_datasets <- list("Figure 2c" = export_figs2c, 
                         "Figure 2d" = export_figs2d,
                         "Figure 2e" = export_figs2e
)

openxlsx::write.xlsx(list_of_datasets, file = file_export_name)

### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
# Figure S5
### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###

file_export_name <- "./ProcessedData/source_data/source_data_fig_s5.xlsx"

export_figs5b <- source_data_2[[7]] %>%
  mutate(Replicate = str_extract(Sample, "R1|R2")) %>%
  dplyr::select(c("name", "position", "Condition", "Replicate", "value")) %>%
  rename(name = "Gene", position = "Position", value = "Allelic Ratio")

export_figs5c <- source_data_2[[8]] %>%
  mutate(Replicate = str_extract(Sample, "R1|R2")) %>%
  dplyr::select(c("name", "Condition", "Replicate", "value", "escape_in_clones")) %>%
  rename(name = "Gene", value = "Allelic Ratio", escape_in_clones = "Escape in all samples")

export_figs5d <- source_data_2[[9]] %>%
  mutate(Replicate = str_extract(Sample, "R1|R2")) %>%
  dplyr::select(c("name", "escape_group", "Condition", "Replicate", "value", "escape_in_clones")) %>%
  rename(name = "Gene", escape_group = "Escape Status", value = "Allelic Ratio", escape_in_clones = "Escape in all samples")

export_figs5e <- source_data_2[[10]] %>%
  dplyr::select(c("name", "escape_status", "Condition", "Replicate", "value")) %>%
  rename(name = "Gene", escape_status = "Escape Status", value = "Allelic Ratio")

export_figs5f <- source_data_2[[11]] %>%
  dplyr::select(c("name", "escape_status", "Condition", "Replicate", "value")) %>%
  rename(name = "Gene", escape_status = "Escape Status", value = "Allelic Ratio")

export_figs5g <- source_data_2[[12]] %>%
  dplyr::select(c("name", "group", "value")) %>%
  rename(name = "Gene", group = "Group", value = "Allelic Ratio")

export_figs5h <- source_data_2[[13]] %>%
  dplyr::select(c("gene", "EscapeAnnotation", "Timepoint", "Control", "Treatment", "pvalue")) %>%
  rename(gene = "Gene", EscapeAnnotation = "Escape Status", Control = "Allelic Ratio Control", Treatment = "Allelic Ratio Treatment")

list_of_datasets <- list("Figure S5b" = export_figs5b, 
                         "Figure S5c" = export_figs5c,
                         "Figure S5d" = export_figs5d,
                         "Figure S5e" = export_figs5e,
                         "Figure S5f" = export_figs5f,
                         "Figure S5g" = export_figs5g,
                         "Figure S5h" = export_figs5h
)

openxlsx::write.xlsx(list_of_datasets, file = file_export_name)


### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
# Figure 3
### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###

file_export_name <- "./ProcessedData/source_data/source_data_fig_3.xlsx"

export_fig3b <- source_data_3[[5]] %>%
  dplyr::select(c("name", "Clone", "Condition2", "TimeSeries", "value")) %>%
  rename(name = "Gene", Condition2 = "Condition", TimeSeries = "Timepoint", value = "Allelic Ratio")

export_fig3c <- source_data_3[[6]] %>%
  dplyr::select(c("Clone", "Condition2", "Timepoint", "Expression")) %>%
  rename(Condition2 = "Condition")

export_fig3d <- source_data_3[[7]] %>%
  dplyr::select(c("name", "Clone", "Condition2", "Timepoint", "value")) %>%
  rename(name = "Gene", Condition2 = "Condition", value = "Allelic Ratio")

export_fig3e <- source_data_3[[8]] %>%
  dplyr::select(c("Gene", "escape_category", "Condition", "x_axis", "y_axis")) %>%
  rename(escape_category = "Escape Status", x_axis = "Allelic Ratio Control", y_axis = "Allelic Ratio Condition")

export_fig3f <- source_data_3[[9]] %>%
  dplyr::select(c("Gene", "escape_category", "Condition", "x_axis", "y_axis")) %>%
  rename(escape_category = "Escape Status", x_axis = "Allelic Ratio Control", y_axis = "Allelic Ratio Condition")

list_of_datasets <- list("Figure 3b" = export_figs3b, 
                         "Figure 3c" = export_figs3c,
                         "Figure 3d" = export_figs3d,
                         "Figure 3e" = export_figs3e,
                         "Figure 3f" = export_figs3f
)

openxlsx::write.xlsx(list_of_datasets, file = file_export_name)


### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
# Figure 6
### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###

file_export_name <- "./ProcessedData/source_data/source_data_fig_6.xlsx"

export_fig6b <- source_data_4[[1]] %>%
  mutate(Clone = gsub("E6 - ", "", Clone)) %>% 
  mutate(Replicate = exp_to_clone[Clone]) %>%
  dplyr::select(c("Condition", "Replicate", "TimeSeries", "Expression")) %>%
  rename(TimeSeries = "Timepoint")

export_fig6d <- source_data_4[[2]] %>%
  mutate(Clone = gsub("E6 - ", "", Clone)) %>% 
  mutate(Replicate = exp_to_clone[Clone]) %>%
  dplyr::select(c("name", "position", "Condition", "Replicate", "value")) %>%
  rename(name = "Gene", position = "Position", "value" = "Allelic Ratio")

export_fig6e <- source_data_4[[3]] %>%
  dplyr::select(c("Gene", "control", "nddox7_washout", "nddox14_washout", "nddox21_washout")) %>%
  mutate_at(vars(matches("dox")), list(ratio = ~ . / control)) %>% 
  dplyr::select(c("Gene", "nddox7_washout_ratio", "nddox14_washout_ratio", "nddox21_washout_ratio")) %>%
  pivot_longer(-c("Gene")) %>% 
  mutate(name = gsub("nd", "", gsub("_washout_ratio", "", name))) %>%
  mutate(name = gsub("dox", "Dox", name)) %>%
  mutate(name = factor(name, levels = paste0("Dox", c(7, 14, 21)))) %>%
  add_column(EscapeAnnotation = rowData(data[.$Gene, ])$EscapeAnnotation) %>%
  dplyr::select(c("Gene", "name", "EscapeAnnotation", "value")) %>%
  rename(name = "Timepoint", EscapeAnnotation = "Escape Status", "value" = "Allelic Ratio")

export_fig6f <- source_data_4[[4]] %>%
  dplyr::select(c("Gene", "name", "Escape", "Ratio", "reversibility_category_joint")) %>%
  rename(name = "Timepoint", reversibility_category_joint = "Reversibility Category")

export_fig6g <- source_data_4[[5]] %>%
  dplyr::select(c("genes", "region_category", "reversibility_index", "halflifes_off")) %>%
  rename(genes = "Gene", region_category = "Escape Group", reversibility_index = "Reversibility", halflifes_off = "Silencing halflife")


list_of_datasets <- list("Figure 6b" = export_fig6b, 
                         "Figure 6d" = export_fig6d, 
                         "Figure 6e" = export_fig6e, 
                         "Figure 6f" = export_fig6f, 
                         "Figure 6g" = export_fig6g
)

openxlsx::write.xlsx(list_of_datasets, file = file_export_name)


### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
# Figure S8
### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###

file_export_name <- "./ProcessedData/source_data/source_data_fig_s8.xlsx"

export_figs8a <- source_data_4[[6]] %>%
  dplyr::select(c("name", "escape_status", "TimeSeries", "Condition", "washout_ratio")) %>%
  rename(name = "Gene", escape_status = "Escape Status", TimeSeries = "Timepoint", washout_ratio = "Washout allelic ratio")

export_figs8b <- source_data_4[[7]] %>%
  dplyr::select(c("name", "TimeSeries", "Condition", "reversibility_category.x", "reversibility_category.y")) %>%
  rename(name = "Gene", TimeSeries = "Timepoint", reversibility_category.x = "Reversibility category (x)", reversibility_category.y = "Reversibility category (y)")

export_figs8c <- source_data_4[[8]] %>%
  dplyr::select(c("TimeSeries", "Condition", Expression)) %>%
  rename(TimeSeries = "Timepoint")

export_figs8d <- source_data_4[[9]] %>%
  dplyr::select(c("TimeSeries", "Condition", Expression)) %>%
  rename(TimeSeries = "Timepoint")

export_figs8e <- source_data_4[[10]] %>%
  dplyr::select(c("name", "escape_status", "TimeSeries", "Condition", "washout_ratio")) %>%
  rename(name = "Gene", escape_status = "Escape Status", TimeSeries = "Timepoint", washout_ratio = "Washout allelic ratio")

export_figs8f <- source_data_4[[11]] %>%
  dplyr::select(c("name", "escape_status", "TimeSeries", "Condition", "washout_ratio")) %>%
  rename(name = "Gene", escape_status = "Escape Status", TimeSeries = "Timepoint", washout_ratio = "Washout allelic ratio")

list_of_datasets <- list("Figure S8a" = export_figs8a, 
                         "Figure S8b" = export_figs8b, 
                         "Figure S8c" = export_figs8c, 
                         "Figure S8d" = export_figs8d, 
                         "Figure S8e" = export_figs8e, 
                         "Figure S8f" = export_figs8f
)

openxlsx::write.xlsx(list_of_datasets, file = file_export_name)


### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
# Figure 7
### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###

file_export_name <- "./ProcessedData/source_data/source_data_fig_7.xlsx"

export_figs7a <- source_data_5[[1]] %>%
  rename(H2AK119ubi = "H2AK119ubi Normalized Coverage", H3K27me3 = "H3K27me3 Normalized Coverage")

export_figs7b <- source_data_5[[2]] %>% 
  rename(value_norm = "Normalized Coverage")

export_figs7c <- source_data_5[[3]] %>% 
  rename(value_norm = "Normalized Coverage")

export_figs7d <- source_data_5[[4]] %>% 
  rename(value_norm = "Normalized Coverage", escape_status = "Escape Status", gene_category = "Reversibility Category")

export_figs7e <- source_data_5_coverage[[1]] %>%
  rename(lfc = "Coverage Log Fold Change", window = "Window")

export_figs7g <- source_data_5_met[[1]] %>%
  rename(coverage = "Number of Reads", beta = "Beta value")

export_figs7h <- source_data_5_met[[2]] %>%
  rename(genes = "Gene", reversibility_category = "Reversibility Category", beta = "Beta value", beta_ref = "Beta value control", beta_diff = "Beta value difference")

list_of_datasets <- list("Figure 7a" = export_figs7a, 
                         "Figure 7b" = export_figs7b, 
                         "Figure 7c" = export_figs7c, 
                         "Figure 7d" = export_figs7d, 
                         "Figure 7e" = export_figs7e, 
                         "Figure 7g" = export_figs7g,
                         "Figure 7h" = export_figs7h
)

openxlsx::write.xlsx(list_of_datasets, file = file_export_name)


### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###
# Figure S9
### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###### ### ###

file_export_name <- "./ProcessedData/source_data/source_data_fig_s9.xlsx"

export_figs9a <- source_data_5_add[[1]] %>%
  dplyr::select(-c("Chr")) %>%
  rename("Pos" = "Position", name2 = "Condition", replicate = "Replicate", value = "Coverage Fold Change")

export_figs9b <- source_data_5_add[[2]] %>%
  dplyr::select(-c("Chr")) %>%
  rename("Pos" = "Position", name2 = "Condition", replicate = "Replicate", value = "Coverage Fold Change")

export_figs9c <- source_data_5[[5]] %>%
  rename(value_norm = "Normalized Coverage")

export_figs9d <- source_data_5[[6]] %>%
  rename(value_norm = "Normalized Coverage")

export_figs9e <- source_data_5[[7]]

export_figs9f <- source_data_3_deg[[4]] %>%
  dplyr::select(c("Gene", "Comparison2", "log2FoldChange", "value")) %>% 
  rename(Comparison2 = "Condition", log2FoldChange = "log2FoldChange RNA", value = "log2FoldChange CnR")

export_figs9g <- source_data_5[[8]] %>%
  ungroup() %>%
  dplyr::select(c("Symbol", "gene_category", "NoDox", `21dDox`, `21dD7dWO`, "FC_dox", "FC_wo")) %>%
  rename(Symbol = "Gene", gene_category = "Reversibility Category", NoDox = "Coverage Control", 
                `21dDox` = "Coverage 21d Dox", `21dD7dWO` = "Coverage 21d Dox + Washout", FC_dox = "Fold Change Dox / Control", FC_wo = "Fold Change Washout / Control")

export_figs9h <- source_data_5[[9]] %>%
  ungroup() %>%
  dplyr::select(c("Symbol", "gene_category", "NoDox", `21dDox`, `21dD7dWO`, "FC_dox", "FC_wo")) %>%
  rename(Symbol = "Gene", gene_category = "Reversibility Category", NoDox = "Coverage Control", 
         `21dDox` = "Coverage 21d Dox", `21dD7dWO` = "Coverage 21d Dox + Washout", FC_dox = "Fold Change Dox / Control", FC_wo = "Fold Change Washout / Control")

export_figs9i <- source_data_5_met[[3]] %>%
  dplyr::select(c("chr", "start", "end", "gene", "is_expressed", "is_escaping", "Condition", "Activity", "coverage", "beta")) %>%
  rename(chr = "Chromosome", start = "Start", end = "End", gene = "Gene", "is_expressed" = "Expression Category", "is_escaping" = "Escape Category",
         coverage = "Read Coverage", beta = "Beta value")
  
list_of_datasets <- list("Figure S9a" = export_figs9a, 
                         "Figure S9b" = export_figs9b, 
                         "Figure S9c" = export_figs9c, 
                         "Figure S9d" = export_figs9d, 
                         "Figure S9e" = export_figs9e, 
                         "Figure S9f" = export_figs9f,
                         "Figure S9g" = export_figs9g,
                         "Figure S9h" = export_figs9h,
                         "Figure S9i" = export_figs9i
)

openxlsx::write.xlsx(list_of_datasets, file = file_export_name)

