#!/bin/bash
set -e
set -o pipefail

# Run bowtie2
# bowtie2 -p ${THREADS} \
# -x ${ref} \
# --local \
# --no-unal \
# --sensitive-local \
# -N 1 \
# -1 ${processed_dir}/trimmed.R1.fastq.gz \
# -2 ${processed_dir}/trimmed.R2.fastq.gz \
# -S processed/temp_aligned_reads.sam

tmp_dir="temp"

# Set up date format for logging
date_format='+%Y%m%d_%H-%M-%S'
date=$(date "$date_format")

# Log file creation with the timestamp
log="${tmp_dir}/temp_log_${date}.log"
touch $log

# Function to log and echo simultaneously
log_and_echo() {
  echo "$@" | tee -a $log
}

# Function to execute a command, log its output and error, and check for success
execute_and_log() {
  start_time=$(date "$date_format")
  log_and_echo "Starting $1 at $start_time"
  
  if ! "$@" 2>&1 | tee -a $log; then
    log_and_echo "Error with $1"
    exit 1
  fi
  
  end_time=$(date "$date_format")
  log_and_echo "Completed $1 at $end_time"
}

# Define variables
sample="CHI-006"
ref="/mnt/d/JorritvU/refgenome/GRCh38/ensemble/Homo_sapiens.GRCh38.dna.primary_assembly"
processed_dir="processed/SCC-scNlaIII-EMC-${sample}"
scripts="/mnt/d/JorritvU/script/BuysDB"
THREADS=24

# Run bowtie2 (command commented out as an example)
# log_and_echo "Aligning starts at $(date "$date_format")"

# # Perform alignment and convert to BAM format while logging stderr
# bwa mem -t ${THREADS} ${ref}.fa ${processed_dir}/trimmed.R1.fastq.gz ${processed_dir}/trimmed.R2.fastq.gz 2>>$log | \
# samtools view -bS - > ${tmp_dir}/temp_aligned_reads.bam

# if [ $? -ne 0 ]; then
#   log_and_echo "Error during alignment."
#   exit 1
# fi

# log_and_echo "Alignment completed at $(date "$date_format")"

# execute_and_log java -jar /home/jorrit/picard.jar SortSam \
#   INPUT=${tmp_dir}/temp_aligned_reads.bam \
#   OUTPUT=${tmp_dir}/sorted.bam \
#   SORT_ORDER=coordinate 

## Add or replace read groups
# execute_and_log java -jar /home/jorrit/picard.jar AddOrReplaceReadGroups \
#   I=${tmp_dir}/sorted.bam \
#   O=${tmp_dir}/sorted.rg.bam \
#   RGID=${sample} \
#   RGLB=NlaIIIseq \
#   RGPL=ILLUMINA \
#   RGPU=BC \
#   RGSM=${sample}

#execute_and_log java -jar /home/jorrit/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${tmp_dir}/sorted.rg.bam O=${tmp_dir}/sorted.rg.dedup.bam METRICS_FILE=${tmp_dir}/sorted.rg.dedup.metrics.txt PROGRAM_RECORD_ID=MarkDuplicates PROGRAM_GROUP_VERSION=null PROGRAM_GROUP_NAME=MarkDuplicates
execute_and_log java -jar /home/jorrit/picard.jar BuildBamIndex -I ${tmp_dir}/sorted.rg.dedup.bam -O ${tmp_dir}/sorted.rg.dedup.bai

log_and_echo "Moving files.. $(date "$date_format")"
mv ${tmp_dir}/sorted.rg.dedup.bam ${processed_dir}/${sample}.bam
mv ${tmp_dir}/sorted.rg.dedup.bai ${processed_dir}/${sample}.bai
mv ${tmp_dir}/sorted.rg.dedup.metrics.txt ${processed_dir}/${sample}.metrics.txt
log_and_echo "Done moving files.. $(date "$date_format")"

# Additional steps can be added here with the execute_and_log function for logging and error checking
# execute_and_log python3 ${scripts}/universalBamTagger/bamtagmultiome.py ${processed_dir}/aligned_reads_sorted.bam -method nla -o ${processed_dir}/tagged/tagged.bam -mapfile ${ref}.mappability.stats.bgzf
