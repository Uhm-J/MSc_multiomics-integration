#!/usr/bin/env Rscript

# Set the R library path with the installed packages
.libPaths("/mnt/d/JorritvU/environment/R/4.3.2/packages")

# Check if the necessary packages are installed
packages <- c("optparse", "BiocManager", "ggplot2", "stringr", "dplyr")
bioc.packages <- c("BSgenome.Hsapiens.UCSC.hg38", "AneuFinder")

# Install packages if not installed
for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# Install Bioconductor packages if not installed
for (package in bioc.packages) {
  if (!require(package, character.only = TRUE)) {
    BiocManager::install(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}


# Load necessary libraries
library(optparse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(AneuFinder)
library(ggplot2)
library(stringr)
library(dplyr)

# Setup optparse to handle command line arguments
option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL, help="Sample name"),
  make_option(c("-c", "--config"), type="character", default=NULL, help="Config path"),
  make_option(c("-m", "--model"), type="character", default="edivisive", help="Model type (edivisive, dnacopy, HMM)"),
  make_option("--run", type="character", default="all", help="Modules to run: all, callCNV, qc, plots, save")
)

# Parse command-line arguments
args = parse_args(OptionParser(option_list=option_list))

# Assign arguments to variables
sample = args$sample
config_path = args$config
model = args$model
modules_to_run <- strsplit(args$run, ",")[[1]]


# QC parameters
qc_params <- list(
  concordanceCutoff = 0.80,
  spikiness = 0.65,
  num.segments = 5, # Not used
  entropy = 0.5, # Not used
  bhattacharyya = 0.4,
  sos = 0.5, # Not used
  total.read.count = 0.5 # Not used
)

# Verify arguments
if (is.null(sample) || is.null(config_path)) {
  stop("Both sample name and config path are required.", call.=FALSE)
}

if (model != "edivisive" && model != "dnacopy" && model != "HMM") {
  stop("Model argument must be one of the following: edivisive, dnacopy, HMM", call.=FALSE)
}

should_run <- function(module) {
  "all" %in% modules_to_run || module %in% modules_to_run
}

# Print the setup details
cat("Running analysis for sample:", sample, "\n")
cat("Using configuration file:", config_path, "\n")
cat("Selected model for analysis:", model, "\n")
cat("Modules to run:", modules_to_run, "\n")

# Setup the environmental variables
basePath <- "/mnt/d/JorritvU/Tripolar/scDNA-seq/"
Input <- paste0(basePath, sample, "/processed/singleBAM/")
outputfolder <- paste0(basePath, sample, "/processed/CNVs/")
outputfolder_plots <- paste0(basePath, sample, "/processed/CNVs/test/")

# Function to retrieve the celltype of the sample
patterns <- c('Bipolar', 'Tripolar', 'Control')
FileLabel <- function(filelist, patterns){
  labels <- sapply(filelist, function(file) {
    for (pattern in patterns) {
      if (str_detect(file, pattern)) {
        return(pattern)
      }
    }
    return(NA) # Return NA or some default value if no pattern matches
  })
  return(labels)
}

# function to generate concordance plot between models
Concordance_plot <- function(data, sample_name, title){
  data <- as.data.frame(data)
  #pdf(file = paste0(outputfolder_plots,names(data),"_Concordance.pdf"))
  print(ggplot(data,aes(x=!!as.name(names(data)), y=!!as.name(names(data)))) + 
          geom_point() + geom_vline(xintercept = 0.9, color = 'red') + 
          geom_hline(yintercept = 0.9, color = 'red') + 
          geom_text(aes(0,0.9,label = 0.9, vjust = -1)) + 
          geom_text(aes(0.9,0,label = 0.9, vjust = -1)) + 
          labs(title = paste0("Compare ", names(data), " Methods"), subtitle = paste0(sample_name, ": ",title)) + 
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)))
  #dev.off()
  #cat("Concordance plot saved to:", plot_file, "\n")
}

# Running AneuFinder for the specified sample
if (should_run("callCNV")) {
  cat("Running AneuFinder for sample:", sample, "\n")
  Aneufinder(inputfolder = Input, outputfolder = outputfolder, configfile = config_path)
}

# Lists to hold file paths
edivisiveFiles <- list()
dnaFiles <- list()
hmmFiles <- list()

# List to hold concordance results
concordanceList <- list()

# Loading file paths
cat("Loading file paths for models...\n")
tryCatch({
  edivisiveFiles[[sample]] <- list.files(paste0(basePath, sample, "/processed/CNVs/MODELS/method-edivisive/"), full = TRUE)
  dnaFiles[[sample]] <- list.files(paste0(basePath, sample, "/processed/CNVs/MODELS/method-dnacopy/"), full = TRUE)
  hmmFiles[[sample]] <- list.files(paste0(basePath, sample, "/processed/CNVs/MODELS/method-HMM/"), full = TRUE)
}, error = function(e) {
  stop("Error while loading file paths:", e$message, ".\n it seems input data is missing, or in another path. Perhaps try running callCNV first?", call.=FALSE)
})

cat("Number of edivisive files:", length(edivisiveFiles), "\n")
cat("Number of dnacopy files:", length(dnaFiles), "\n")
cat("Number of HMM files:", length(hmmFiles), "\n")

# Performing comparisons
cat("Performing comparisons between models...\n")

# Comparison 1: edivisive vs. dnacopy
concordanceList[[paste0(sample, "_edivisive_vs_dnacopy")]] <- compareMethods(edivisiveFiles[[sample]], dnaFiles[[sample]])

# Comparison 2: edivisive vs. HMM
concordanceList[[paste0(sample, "_edivisive_vs_HMM")]] <- compareMethods(edivisiveFiles[[sample]], hmmFiles[[sample]])

# Comparison 3: dnacopy vs. HMM
concordanceList[[paste0(sample, "_dnacopy_vs_HMM")]] <- compareMethods(dnaFiles[[sample]], hmmFiles[[sample]])

if(should_run("plots")) {
  cat("Generating concordance plots...\n")
  pdf(file = paste0(outputfolder_plots, sample, "_Concordance.pdf"))
  par(mfrow = c(1,3))
  Concordance_plot(concordanceList[[paste0(sample, "_edivisive_vs_dnacopy")]], sample, "Edivisive vs dnacopy")
  Concordance_plot(concordanceList[[paste0(sample, "_edivisive_vs_HMM")]], sample, "Edivisive vs HMM")
  Concordance_plot(concordanceList[[paste0(sample, "_dnacopy_vs_HMM")]], sample, "dnacopy vs HMM")
  dev.off()
  cat("Concordance plot saved to:", paste0(outputfolder_plots, sample, "_Concordance.pdf"), "\n")
}


# Select best model based on user input
cat("Selecting the best model based on user preference...\n")
bestModelFiles <- switch(model,
edivisive = edivisiveFiles,
dnacopy = dnaFiles,
HMM = hmmFiles,
stop("Model not recognized", call.=FALSE))


# Quality control
if (should_run("qc")) {
    cat("Performing quality control...\n")

    # Print QC parameters
    cat("Quality control parameters:\n")
    cat("Spikiness threshold:", qc_params$spikiness, "\n")
    cat("Number of segments threshold:", qc_params$num.segments, "\n")

    # Curating the list of files based on concordance
    cat("Curating files based on concordance threshold... (", qc_params$concordanceCutoff,")\n")

    curatedList <- list()
    curatedList[[sample]] <- bestModelFiles[[sample]][concordanceList[[sample]]$concordance > qc_params$concordanceCutoff]
    print(length(curatedList[[sample]]))
    cat("Number of files after applying concordance cutoff:", length(curatedList), "\n")
    cat("Proportion of files retained after curation:", length(curatedList[[sample]]) / length(bestModelFiles[[sample]]), "\n")
    
    qcList <- list()
    cat("Performing quality control on curated files...\n")
    qc <- clusterByQuality(curatedList[[sample]], measures=c('spikiness','num.segments','entropy','bhattacharyya','sos', 'total.read.count'))
    print(qc$parameters)
    data <- as.data.frame(qc$Mclust$data)
    filterData <- data%>%filter(spikiness < qc_params$spikiness & bhattacharyya > qc_params$bhattacharyya)
    qcList[[sample]] <- row.names(filterData)
    # Print post-QC stats
    cat("Number of files after quality control:", length(qcList[[sample]]), "\n")
    cat("Proportion of files retained after quality control:", length(qcList[[sample]]) / length(curatedList[[sample]]), "\n")
} else {
    qcList <- list()

    qcList[[sample]] <- bestModelFiles[[sample]]
}

# Generate heatmap
if (should_run("plots")) {
    cat("Generating heatmap...\n")
    label <- FileLabel(qcList[[sample]], patterns)
    print(factor(label))
    heatmapGenomewide(qcList[[sample]], classes = factor(label), classes.color = c('red','blue'), file = paste0(outputfolder_plots, model,"_binsize_5Mb_genomewide_0.90dnacopy.pdf"))
}


# Generate karyotype measures
if (should_run("plots")) {
    cat("Generating karyotype measures...\n")
    files.tripolar <- qcList[[sample]][str_detect(qcList[[sample]], "_Tripolar")]
    files.bipolar <- qcList[[sample]][!str_detect(qcList[[sample]], "_Bipolar")]
    files.control <- qcList[[sample]][str_detect(qcList[[sample]], "_Control")]

    plot_file_path <- paste0(outputfolder_plots, "Heterogeneity_Plot_", sample, ".pdf")
    pdf(file = plot_file_path)
    plotHeterogeneity(hmms.list = list(Tripolar = files.tripolar, Bipolar = files.bipolar, Control = files.control))
    dev.off()
    cat("Heterogeneity plot saved to:", plot_file_path, "\n")
}


# Save the curated list of files
if(should_run("save")) {
    out_name <- paste0(outputfolder, sample)
    cat("Saving curated list of files to: ", out_name,"\n")
    exportCNVs(qcList[[sample]], out_name)
}
cat("Done!\n")
