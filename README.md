# README
'---
(4th time)
title: "Acute Myeloid Leukemia Heatmap"
author: "CCDL for ALSF - Adapted for this repository by Candace Savonen"
date: "October 2021"
output:
html_notebook:
toc: true
toc_float: true
number_sections: true
This analysis has been adapted from this refine.bio-examples notebook

Purpose of the analysis
In this analysis, we will use this acute myeloid leukemia sample dataset from Shih et al., 2017 and pre-processed by refinebio.

The data that we downloaded from refine.bio for this analysis has 19 samples (obtained from 19 acute myeloid leukemia (AML) model mice), containing RNA-sequencing results for types of AML under controlled treatment conditions.

This dataset can be downloaded from this page on refine.bio.
We will download it already processed and quantile normalized.

Set up analysis folders
python

Copy
# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
applescript

Copy
# Define the file path to the data directory
data_dir <- file.path("data", "SRP070849")

# Declare the file path to the gene expression matrix file
data_file <- file.path(data_dir, "SRP070849.tsv")

# Declare the file path to the metadata file
# inside the directory saved as `data_dir`
metadata_file <- file.path(data_dir, "metadata_SRP070849.tsv")
Clustering Heatmap - RNA-seq
Install libraries
We will use []pheatmap by Slowikowski et al., 2017](https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap) for clustering and creating a heatmap.


Copy
if (!("pheatmap" %in% installed.packages())) {
  # Install pheatmap
  install.packages("pheatmap", update = FALSE)
}
Attach the pheatmap and magrittr libraries:

apache

Copy
# Attach the `pheatmap` library
library(pheatmap)

# We will need this so we can use the pipe: %>%
library(magrittr)

# Set the seed so our results are reproducible:
set.seed(12345)
Import and set up data
This chunk of code will read in both TSV files and add them as data frames to your environment.

pgsql

Copy
# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)

# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  # Here we are going to store the gene IDs as row names so that
  # we can have only numeric values to perform calculations on later
  tibble::column_to_rownames("Gene")
Let's take a look at the metadata object that we read into the R environment.


Copy
head(metadata)
Now let's ensure that the metadata and data are in the same sample order.

nsis

Copy
# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)
Now we are going to use the pheatmap package to look at how are samples and genes are clustering.

Choose genes of interest
For this example, we will sort genes by variance and select genes in the upper quartile, but there are many alternative criterion by which you may want to sort your genes, <i>e.g.</i> fold change, t-statistic, membership in a particular gene ontology, so on.

r

Copy
# Calculate the variance for each gene
variances <- apply(expression_df, 1, var)

# Determine the upper quartile variance cutoff value
upper_var <- quantile(variances, 0.75)

# Filter the data choosing only genes whose variances are in the upper quartile
df_by_var <- data.frame(expression_df) %>%
  dplyr::filter(variances > upper_var)
Let's save these to our results folder as a TSV file.


Copy
readr::write_tsv(df_by_var, file.path(results_dir, "top_90_var_genes.tsv"))
Prepare metadata for annotation
From the accompanying paper, we know that the mice with IDH2 mutant AML were treated with vehicle or AG-221 (the first small molecule in-vivo inhibitor of IDH2 to enter clinical trials) and the mice with TET2 mutant AML were treated with vehicle or 5-Azacytidine (Decitabine, hypomethylating agent) Shih et al., 2017.

php

Copy
# Let's prepare the annotation for the uncollapsed `DESeqData` set object
# which will be used to annotate the heatmap
annotation_df <- metadata %>%
  # Create a variable to store the cancer type information
  dplyr::mutate(
    mutation = dplyr::case_when(
      startsWith(refinebio_title, "TET2") ~ "TET2",
      startsWith(refinebio_title, "IDH2") ~ "IDH2",
      startsWith(refinebio_title, "WT") ~ "WT",
      # If none of the above criteria are satisfied,
      # we mark the `mutation` variable as "unknown"
      TRUE ~ "unknown"
    )
  ) %>%
  # select only the columns we need for annotation
  dplyr::select(
    refinebio_accession_code,
    mutation,
    refinebio_treatment
  ) %>%
  # The `pheatmap()` function requires that the row names of our annotation
  # data frame match the column names of our `DESeaDataSet` object
  tibble::column_to_rownames("refinebio_accession_code")
Create annotated heatmap
ini

Copy
# Create and store the annotated heatmap object
heatmap_annotated <-
  pheatmap(
    df_by_var,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    annotation_col = annotation_df, # Specify our annotation here
    main = "Annotated Heatmap",
    colorRampPalette(c(
      "deepskyblue",
      "black",
      "yellow"
    ))(25
    ),
    scale = "row" # Scale values in the direction of genes (rows)
  )
Save annotated heatmap as a PNG
You can switch this to save to a JPEG or TIFF by changing the function and file name within the function to the respective file suffix.

cmake

Copy
# Open a PNG file
png(file.path(
  plots_dir,
  "aml_heatmap.png" # Replace with a relevant file name
))

# Print the heatmap
heatmap_annotated

# Close the PNG file:
dev.off()
Session info

Copy
# Print session info
sessioninfo::session_info()
'

Bot image for Assistant
Assistant

Poe
README for Acute Myeloid Leukemia Heatmap Script
Overview
This script, 01-heatmap.Rmd, is designed to create a clustered heatmap visualization using RNA-seq data from a study on Acute Myeloid Leukemia (AML). The analysis utilizes a dataset from Shih et al., 2017 that is available on refine.bio.

Purpose
The primary objective of this analysis is to explore transcriptional patterns in AML samples by generating a heatmap that clusters genes and samples based on their expression profiles. This will help visualize the differences and similarities in gene expression among various treatment conditions.

Prerequisites
Before running this script, ensure you have the following:

R installed on your machine.
RStudio or any compatible IDE for R.
Required R packages: pheatmap, magrittr, dplyr, readr, tibble, and sessioninfo.
Installation of Required Packages
Use the following command to install any missing packages:

r

Copy
if (!("pheatmap" %in% installed.packages())) {
  install.packages("pheatmap", update = FALSE)
}
File Structure
The script expects a specific directory structure:


Copy
project_root/
├── data/
│   └── SRP070849/
│       ├── SRP070849.tsv
│       └── metadata_SRP070849.tsv
├── plots/
├── results/
└── 01-heatmap.Rmd
data/: Contains input data files.
plots/: Will store the generated heatmap images.
results/: Will contain processed data files.
Workflow
Setup Analysis Folders: The script will create necessary folders for data, plots, and results if they don't already exist.
Load Libraries: Required libraries are loaded, and a reproducible seed is set.
Import Data: The script reads the gene expression and metadata files into R.
Data Processing:
The script filters genes based on variance to focus on the most variable genes, which are likely to be biologically relevant.
It prepares metadata for annotation based on the treatment conditions.
Generate Heatmap: A clustered heatmap is created using the pheatmap function, with annotations for treatment conditions added to the plot.
Save Heatmap: The heatmap is saved as a PNG file in the plots directory.
Session Information: The script prints session information at the end for reproducibility.
Usage
To run the analysis:

Open the 01-heatmap.Rmd file in RStudio.
Execute each code chunk sequentially by clicking the "Run" button or using the keyboard shortcut (Ctrl + Shift + Enter).
Check the plots/ directory for the saved heatmap.
References
Dataset: refine.bio SRP070849
Publication: Shih et al., 2017
pheatmap Documentation: pheatmap RDocumentation
Contact
For questions or feedback, please reach out to Candace Savonen at the Childhood Cancer Data Lab (CCDL) for ALSF.

This README provides a comprehensive guide to using the 01-heatmap.Rmd script for generating a heatmap based on AML RNA-seq data. Follow the instructions carefully to ensure successful execution of the analysis.