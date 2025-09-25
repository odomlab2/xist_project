### escapees new studies

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(EnsDb.Mmusculus.v79)
all_genes <- genes(EnsDb.Mmusculus.v79)
symbol_to_ensembl <- setNames(all_genes$gene_id, nm = all_genes$symbol)

### load excel file
excel_sheet <- readxl::read_excel("../../../ProcessedData/ListOfEscapeeFromEdithLab.xlsx")
head(excel_sheet)
genes <- excel_sheet$ENSEMBL_v102
genes <- setNames(rep(NA, length(genes)), nm = genes)

# Xist exerts gene-specific silencing during XCI maintenance and 
# impacts lineage-specific cell differentiation and proliferation during hematopoiesis
data_here <- readxl::read_excel("../../../ProcessedData/41467_2022_32273_MOESM6_ESM.xlsx", skip = 1)
data_here$ensembl <- symbol_to_ensembl[data_here$`Gene name`]
data_here <- data_here %>% mutate(XCI_status = case_when(
  `XCI status` == "NA" ~ NA, 
  `XCI status` == "Subjective" ~ "S", 
  `XCI status` == "Escape" ~ "E"
))

data_here %>% 
  dplyr::filter(`XCI status` == "Escape") %>%
  pull(`Gene name`)

# add to big table
excel_sheet$Yang2022_Status <- as.character(data_here[match(excel_sheet$ENSEMBL_v102, data_here$ensembl), "XCI_status"]$XCI_status)

### MEF study
data_here <- readxl::read_excel("../../../ProcessedData/mmc2(1).xlsx", skip = 0)
data_here$ensembl <- symbol_to_ensembl[data_here$Geneid]
data_here <- data_here %>% mutate(XCI_status = case_when(
  SmcHD1_dependence == "MEF_escapee" ~ "E", 
  .default = "S"
))

table(data_here$XCI_status)
excel_sheet$Bowness22_Status <- as.character(data_here[match(excel_sheet$ENSEMBL_v102, data_here$ensembl), "XCI_status"]$XCI_status)

write.csv(excel_sheet, "../../../ProcessedData/NewEscapeeClas/ListOfEscapeeFromEdithLabUpdated.csv")

