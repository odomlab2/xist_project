module load STAR

index_directory=/omics/groups/OE0538/internal/users/panten/projects/genome_files/star_indices/B6_masked_JF1_only_rtta/
fasta_files=/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/rtta_alignment/rtta.fa

STAR --runThreadN 6 \
    --runMode genomeGenerate \
    --genomeDir $index_directory \
    --genomeFastaFiles $fasta_files
#    --sjdbGTFfile /omics/groups/OE0538/internal/users/panten/projects/genome_files/genome_fastas/GRCm39/gtf/Mus_musculus.GRCm39.109.gtf \
#    --sjdbOverhang 99
