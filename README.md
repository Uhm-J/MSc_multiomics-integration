# MSc Multiomics-Integration Documentation

## Exploring Single Nucleotide Variants as shared features for unpaired Genome-Transcriptome integration

## Overview
This repository serves as a comprehensive archive for the codes utilized in the MSc Thesis focused on the exploration of multi-omics integration, with a special emphasis on Single Nucleotide Variants (SNVs). The study delves into the complex interplay between genomic and transcriptomic data, highlighting the potential and limitations of current methodologies. Despite the anticipated utility of genomic markers in integration efforts, significant challenges persist, primarily due to the static nature of genomic data against the dynamic transcriptome. This thesis concludes with a call for innovative computational techniques and a broader perspective on multi-omics integration to advance our understanding of biological systems.

![image](https://github.com/Uhm-J/MSc_multiomics-integration/assets/115709478/3459798b-fa4e-4a1a-8105-f02390e181ef)

## Methods

### 1. Aligning WGS
The process for preparing Whole Genome Sequencing (WGS) data involves multiple steps, beginning with demultiplexing and trimming, based on the [BuysDB](https://github.com/BuysDB/SingleCellMultiOmics) scripts. The snakefile is modified so that rule_all only needs the trimmed reads in order to run custom aligning. Following these preparatory steps, the `align.sh` script is executed for aligning the WGS data using BWA and Samtools, targeting the GRCh38.p14 human reference genome. Subsequent steps include sorting and adding read groups with Picard, marking duplicates to refine the dataset further, and creating an index to facilitate efficient data retrieval. 

### 2. Aligning and Counting RNA
`star_align.sh` facilitates the alignment and quantification of RNA sequencing data via STAR aligner, emphasizing the integration of single-cell RNA sequencing (scRNAseq) to dissect cellular heterogeneity. The script includes steps for gene expression quantification across the entire gene body, enabling a detailed examination of transcriptomic activity. It requires a parameter `--SAMPLE` or `-s`. Within the script there are some folder location to be modified.

### 3. Splitting BAM Files (RNA AND DNA)

The `split_bamfiles.py` and `split_bamfiles_rna.py` scripts segment BAM files into cell-specific files, enabling precise genomic and transcriptomic analyses at the single-cell level. This process is essential for accurately interpreting data based on cell barcodes and unique molecular identifiers (UMIs).

#### Dependencies:
- `pysam`: A python module for reading, manipulating, and writing genomic data sets.
- `os`: A module for interacting with the operating system.
- `argparse`: A parser for command-line options, arguments, and sub-commands.
- `re`: A module for regular expression matching operations.
- `tqdm`: A fast, extensible progress bar for loops and CLI.
- `subprocess`: A module for running new applications or programs through Python.
- `pandas`: A fast, powerful, flexible, and easy-to-use open-source data analysis and manipulation tool.

#### Usage:

##### split_bamfiles.py and split_bamfiles_rna.py
```bash
python split_bamfiles.py \
  --input <input_bam> \
  --output_dir <output_directory> \
  --sample <sample_name> \
  --annotation_file <annotation_excel_file> \
  --skip_rows <number_of_rows_to_skip> \
python split_bamfiles_rna.py \
  --input <input_bam> \
  --output_dir <output_directory> \
  --sample <sample_name> \
  --annotation_file <annotation_excel_file> \
  --skip_rows <number_of_rows_to_skip>
```

- `--input`: Path to the input BAM file.
- `--output_dir`: Directory where the split BAM files will be saved.
- `--sample`: Name of the sample being processed.
- `--annotation_file`: Excel file containing annotations for cell identifiers and phenotypes.
- `--skip_rows`: Number of initial rows to skip in the annotation file.

These scripts employ pysam for BAM file manipulation, utilizing argparse for command-line interface creation, and pandas for handling the annotations provided in Excel format. The `re` module is used for parsing cell identifiers and UMIs from read names, and `subprocess` runs external commands for adding or replacing read groups with Picard. Progress feedback is provided via `tqdm`. Admittedly, these files and parameters are highly customized for the infrastructure present in the lab. 

For both scripts, ensure that the required dependencies are installed, which can be done using pip:
```bash
pip install pysam argparse pandas tqdm
```

Additionally, `java` and `picard.jar` need to be available in your environment to execute Picard commands for read group manipulation.

### 4. Monopogen
Employing `Monopogen.sh`, Single Nucleotide Variants (SNVs) are called and processed as shared features for integration. This step highlights Monopogen's efficiency in detecting SNVs across scWGS and scRNA-seq datasets, setting the stage for comprehensive multi-omic data integration.

### 5. ConvertVCF to h5ad (RNA and DNA)
The Python scripts `convertVCF_h5ad.py` (DNA) and `convertVCF_h5ad_rna.py` (RNA) convert variant call format (VCF) files into h5ad format. This conversion facilitates the organization of genetic variants and metadata for integrated analysis, emphasizing the significance of SNVs in the context of cellular phenotypes and tumor heterogeneity.

### 6. Omic-Specific Clusters
#### 6.1 Aneufinder_caller.R
This script utilizes AneuFinder for CNV calling in scWGS data, contributing to the understanding of chromosomal instability (CIN) and its implications for cancer progression. The analysis underscores the need for high-quality data and rigorous QC to accurately profile CNVs and their impact on cellular phenotypes.

#### 6.2 DEG Limma + Seurat
Differential expression analysis, combining Limma with Seurat, focuses on identifying differentially expressed genes (DEGs) that might be influenced by CIN-related genomic alterations. This approach enriches the integration of scWGS and scRNAseq data by associating genomic features with transcriptomic profiles.

### 7. SNV Profile Check
The `SNV_profile.ipynb` notebook is dedicated to the examination and selection of SNVs for downstream integration. This step involves quality control, SNV calling, and the application of criteria to select SNVs that best represent the biological variations within the dataset.

### 8. Integration
`Integration_scRNA-DNA.ipynb` outlines the methodology for integrating scWGS and scRNAseq data using SIMBA, aiming to co-embed cells along with genomic features. This integration enables a nuanced exploration of the relationship between genotype, gene expression, and CNVs, offering new insights into tumor heterogeneity and the dynamics of cancer progression.
