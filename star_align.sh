
#!/bin/bash

# BELOW ARE MORE VARIABLES!!
# Initialize Sample variable (this is modified by the flags used)
SAMPLE=""

# Manual parsing of command line arguments
for arg in "$@"; do
  shift
  case "$arg" in
    "--SAMPLE") set -- "$@" "-s" ;;
    *)        set -- "$@" "$arg"
  esac
done

# Extract options and their arguments into variables.
while getopts s: flag
do
  case "${flag}" in
    s) SAMPLE=${OPTARG};;
  esac
done

if [ -z "$SAMPLE" ]; then
  echo "No SAMPLE specified. Exiting.\nUse --SAMPLE or -s to provide SAMPLE name."
  exit 1
fi

ext_dir="/mnt/d/JorritvU/Tripolar/scRNA-seq/${SAMPLE}"  # the script extends this with /Fastq/*  --> Change the code below if this is unwanted behaviour
ref_gi59="refgenome/Grch38/ensemble/indexHuman59"  
ref_g38="refgenome/Grch38/ensemble/Homo_sapiens.GRCh38.gtf"
bc_file="refgenome/Grch38/ensemble/Homo_sapiens.GRCh38.gtf"



# Copy the merged fastq to current DIR
echo "Copying $ext_dir/Fastq/*.fastq.gz.." 
cp ${ext_dir}/Fastq/*.fastq.gz ./

# Set max amount open files
ulimit -n 65535

# Run STAR
STAR --runMode alignReads \
--runThreadN 32 \
--genomeDir $ref_gi59 \
--sjdbGTFfile $ref_g38 \
--sjdbOverhang 59 \
--readFilesCommand zcat \
--readFilesIn *R2_001.fastq.gz *R1_001.fastq.gz \
--soloType CB_UMI_Simple \
--soloCBstart 7 \
--soloCBlen 8 \
--soloUMIstart 1 \
--soloUMIlen 6 \
--soloBarcodeMate 2 \
--soloCBwhitelist $bc_file \
--outFilterScoreMinOverLread 0.5 \
--outFilterMatchNminOverLread 0.5 \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ./STAR/${SAMPLE} \
--soloOutFileNames _solo \
--soloStrand Unstranded \
--soloMultiMappers Unique EM \
--alignIntronMax 1 \
--soloUMIdedup Exact \
--clip5pNbases 36 0 \
--soloCellFilter EmptyDrops_CR \
--soloFeatures GeneFull

# Everything from --alignIntroMax -> down, is added based on experimentation

# Move files to harddisk in order to save space
mv ./STAR $ext_dir
rm *.fastq.gz
