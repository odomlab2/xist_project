

module load STAR

sample_id=s1

DATA_DIR=/omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/data/20220307_invivo_exp1/2023-03-07-AACKYJCM5/

fastq_1=$DATA_DIR/AACKYJCM5_E35_dox_Feb23_23s000758-1-1_loda_lane1${sample_id}_1_sequence.txt.gz
fastq_2=$DATA_DIR/AACKYJCM5_E35_dox_Feb23_23s000758-1-1_loda_lane1${sample_id}_2_sequence.txt.gz

ls $fastq_1
ls $fastq_2

mkdir -p /omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/processed_data/rtta_star_align/${sample_id}

STAR --genomeDir /omics/groups/OE0538/internal/users/panten/projects/genome_files/star_indices/B6_masked_JF1_only_rtta/ \
    --runThreadN 6 \
    --readFilesIn $fastq_1 $fastq_2 \
    --outFileNamePrefix /omics/groups/OE0538/internal/users/panten/projects/AntoniaESCProject/invivo_data/processed_data/rtta_star_align/${sample_id}/run \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard  \
    --readFilesCommand zcat
