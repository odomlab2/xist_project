### make sampleDescription.txt file for second invivo experiment

data.dir <- "/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/data/20231213_invivo_exp4/2023-12-11-AAFCTG2M5/"

read1.files <- list.files(data.dir, pattern = "1_sequence.txt.gz$") %>% sort()
read2.files <- gsub("1_sequence.txt.gz", "2_sequence.txt.gz", read1.files)

read1.files <- paste0(data.dir, read1.files)
read2.files <- paste0(data.dir, read2.files)

sample.name <- basename(gsub("_1_sequence.txt.gz$", "", read1.files))

df_out <- data.frame(
  "sample" = sample.name,
  "reads1" = read1.files, 
  "reads2" = read2.files
)

write.csv(df_out, "/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/sampleDescription.Exp4.txt", 
          row.names = F, quote = F)