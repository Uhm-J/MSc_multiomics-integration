# MSc Research Project

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
Employing `Monopogen.sh`, Single Nucleotide Variants (SNVs) are called and processed as shared features for integration. This step highlights [Monopogen](https://github.com/KChen-lab/Monopogen)'s efficiency in detecting SNVs across scWGS and scRNA-seq datasets, setting the stage for comprehensive multi-omic data integration. Within the script, there are a couple of config variables that need to be tuned to suit the infrastructure.

### 5. ConvertVCF to h5ad (RNA and DNA)

The Python scripts `convertVCF_h5ad.py` (DNA) and `convertVCF_h5ad_rna.py` (RNA) are designed to convert variant call format (VCF) files into h5ad format, facilitating the organization of genetic variants and metadata for integrated analysis. This process emphasizes the significance of Single Nucleotide Variants (SNVs) in the context of cellular phenotypes and tumor heterogeneity.

#### Dependencies:
- `argparse`: For parsing command-line options.
- `gzip`: To read `.gz` compressed VCF files.
- `pandas`: For data manipulation and analysis.
- `anndata`: For creating an h5ad file which is commonly used in single-cell genomics.
- `scipy`: Specifically `scipy.sparse`, for efficient storage of the genotype matrix.

#### Usage:

```bash
python convertVCF_h5ad.py --vcf path/to/your/file.vcf --out output_filename_without_extension --batch batch_name_optional

python convertVCF_h5ad_rna.py --vcf path/to/your/file.vcf --out output_filename_without_extension --batch batch_name_optional
```

- `--vcf`: Path to the input VCF file. The script supports `.vcf` and `.vcf.gz` formats.
- `--out`: The output filename for the h5ad file without the extension. The script automatically appends `.h5ad` to the specified filename.
- `--batch` (optional): A string identifier for the batch. This can be used to prepend a batch name to the cell identifiers.

These scripts perform the following key steps:
- Parse the VCF file to extract SNV information and genotype data.
- Preprocess the data to convert genotype strings to numerical format suitable for analysis.
- Create an AnnData object with sparse matrix representation of the genotype data, along with SNV and cell annotations.
- Write the AnnData object to an h5ad file for further analysis.

Please ensure that you have the required dependencies installed in your environment, which can be done using pip:

```bash
pip install argparse pandas anndata scipy
```

### 6. Omic-Specific Clusters
#### 6.1 Aneufinder_caller.R
This script utilizes AneuFinder for CNV calling in scWGS data, contributing to the understanding of chromosomal instability (CIN) and its implications for cancer progression. The analysis underscores the need for high-quality data and rigorous QC to accurately profile CNVs and their impact on cellular phenotypes.

#### 6.2 DEG Limma + Seurat
Differential expression analysis, combining Limma with Seurat, focuses on identifying differentially expressed genes (DEGs) that might be influenced by CIN-related genomic alterations. This approach enriches the integration of scWGS and scRNAseq data by associating genomic features with transcriptomic profiles.

### 7. SNV Profile Check
The `SNV_profile.ipynb` notebook is dedicated to the examination and selection of SNVs for downstream integration. 

#### SNV Selection Process

The SNV selection process is integral for identifying Single Nucleotide Variants (SNVs) that serve as a reliable foundation for the integration of genomic (DNA) and transcriptomic (RNA) datasets. This methodical approach aims to highlight SNVs that are consistent across these datasets, ensuring their biological relevance.

#### Steps in the Selection Process:

The SNV selection process is a critical step in ensuring the robustness and relevance of SNVs for integration between genomic (DNA) and transcriptomic (RNA) datasets. This selection is based on a comparative analysis of genotypes across these datasets, aiming to identify SNVs that exhibit consistent and biologically meaningful patterns.

#### Steps in the Selection Process:

1. **Genotype Proportion Calculation**: For each SNV, we calculate the proportion of each genotype (represented as 1 for homozygous reference, 2 for heterozygous variant, and 3 for homozygous variant) within both DNA and RNA datasets. This proportion is the count of a specific genotype for an SNV divided by the total counts of all genotypes for that SNV within the dataset.

2. **Ratio Calculation**: We then calculate the ratio of genotype proportions between DNA and RNA datasets for each genotype. This ratio helps in assessing the consistency of genotypic expression across genomic and transcriptomic levels.

3. **Applying Constraints**: To ensure the selection of SNVs with biologically relevant variations, we apply constraints on the ratio and proportions:
   - The ratio must fall within a specified interval around 1, indicating similar proportions of a genotype in both datasets.
   - Genotype proportions must exceed a minimum threshold, ensuring that only SNVs with significant presence are considered.
   - Proportions must also demonstrate sufficient variation within genotypes to exclude non-informative SNVs.

#### Criteria for Selection:

- **Allowed Variation**: A parameter (`a`) defines the allowed variation of proportions between datasets, ensuring that only SNVs with consistent representation across DNA and RNA are selected.
- **Minimum Threshold**: A threshold (`b`) filters out noise by disregarding SNVs with very low proportions, focusing on those with significant biological relevance.
- **Intra-genotype Variation**: A parameter (`c`) is set to ensure that selected SNVs show enough variation within their genotypes, excluding potential sequencing errors or artifacts.

#### Visualization and Integration:

Heatmaps are generated to visualize the selection process, showcasing the SNVs that meet the selection criteria versus those that are rejected. This visual representation aids in fine-tuning the selection parameters and understanding the distribution of SNVs across datasets. The final set of selected SNVs, satisfying all criteria, is then ready for downstream integration, contributing to a comprehensive multi-omic analysis.

This selective approach ensures that the integration of genomic and transcriptomic data is based on SNVs that are not only consistent across datasets but also biologically informative, enhancing the quality and interpretability of the integration results.


### 8. Integration
`Integration_scRNA-DNA.ipynb` outlines the methodology for integrating scWGS and scRNAseq data using [SIMBA](https://simba-bio.readthedocs.io/en/latest/), aiming to co-embed cells along with genomic features. This integration enables a nuanced exploration of the relationship between genotype, gene expression, and CNVs.
