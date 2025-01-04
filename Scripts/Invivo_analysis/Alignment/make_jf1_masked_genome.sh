source ~/miniconda3/bin/activate bulk_allelic_nextflow

snp_file=/omics/groups/OE0538/internal/users/panten/projects/genome_files/snp_files/mgp_REL2021_snps.vcf.gz
ref_file=/omics/groups/OE0538/internal/users/panten/projects/genome_files/genome_fastas/GRCm39/fasta_removechr/

SNPsplit_genome_preparation --vcf $snp_file --ref $ref_file --strain JF1_MsJ
