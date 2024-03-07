#!/bin/bash

ulimit -Sn 100000
ulimit -n 100000


# Config settings
Monopogen=/home/jorrit/Monopogen  # Install Monopogen as usual, and replace path with location
apps=/home/jorrit/Monopogen/apps  # This is where Monopogen's tools are located (e.g. bcftools)
regions=/home/jorrit/Monopogen/resource/GRCh38.region.lst  # Chromosomal regions per chromosome, to digest the files in acceptable chunks
ref_genome=/mnt/d/JorritvU/refgenome/GRCh38/ensemble/ref_with_chr/Homo_sapiens.GRCh38.dna.primary_assembly.fa  # Ref genome with CHR prefixes
ref_panel=/mnt/d/JorritvU/refgenome/1KG3/  # References panel, downloaded from recommended link
Monopogen=/home/jorrit/Monopogen/src/Monopogen.py  # The Monopogen run file, runs on python
OUTPUT=$PWD/processed/SNV  # Output folder
THREADS=12

# If no index files present, create index
for bam in processed/singleBAM/*.bam; do
	filename=$(basename "$bam" .bam);
	if [ ! -f "processed/singleBAM/${filename}.bam.bai" ]; then
		echo "Index file for $bam not found. Creating index..."
        java -jar /home/jorrit/picard.jar BuildBamIndex -I "$bam" -O processed/singleBAM/$filename.bam.bai
	fi
done

# Creating sample list
echo "Creating sample list.."
for bam in processed/singleBAM/*.bam; do     
	filename=$(basename "$bam" .bam); 
	echo "$filename,$PWD/$bam"; 
done > processed/singleBAM/bam.lst

# Preprocess the BAM files in bam.lst
/home/jorrit/anaconda3/bin/python $Monopogen preProcess -b processed/singleBAM/bam.lst -a $apps -o $OUTPUT -t $THREADS

# Running germline variant calling
/home/jorrit/anaconda3/bin/python $Monopogen germline -a $apps -t $THREADS -r $regions -p $ref_panel -g $ref_genome -m 3 -s all -o $OUTPUT

