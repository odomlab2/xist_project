### run rtta alignments and quantify 

library(bsub)

data.dir <- "/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/data/20231213_invivo_exp4/2023-12-11-AAFCTG2M5/"
samples <- gsub("_1_sequence.txt.gz", "", list.files(data.dir, pattern = "_1_sequence.txt.gz$"))

alignment_script <- "/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/scripts/bowtie_rtta_alignment_exp4.sh"

for (sample_here in samples) {
  
  print(sample_here)
  
  #bsub_chunk(name = paste0('map_to_rtta_', sample_here), variables = c("sample_here"), memory = 100, core = 4, hour = 36,{
  
  alignment_script <- "/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/scripts/bowtie_rtta_alignment_exp4.sh"
  
  system(paste0("bash ", alignment_script, " ", sample_here))
  
  #})
}

for (sample_here in samples) {
  
  print(sample_here)
  
  #bsub_chunk(name = paste0('count_rtta_', sample_here), variables = c("sample_here"), memory = 100, core = 4, hour = 36, {
  
  counting_script <- "/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/scripts/get_rtta_reads.sh"
  bam_file <- paste0("/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/rtta_alignment/output/exp4/", sample_here, ".rtta.sam")
  
  system(paste0("bash ", counting_script, " ", bam_file))
  
  #})
}

quant_files <- list.files("/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/rtta_alignment/output/exp4/", pattern = ".rtta.sam.rtta_reads", full.names = T)
rtta_quants <- as.numeric(unlist(lapply(quant_files, function(x){readLines(x)})))
names(rtta_quants) <- gsub(".rtta.sam.rtta_reads", "", basename(quant_files))

data.frame(
  Sample = names(rtta_quants), 
  RttaCounts = rtta_quants
) %>% write_csv("/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/rtta_alignment/output/exp4/merged.rtta.counts")


