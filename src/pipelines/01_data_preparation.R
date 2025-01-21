# @ Rik Verdonck 20250121

# Purpose: Handles raw data and metadata import, initial checks, and setup.

# Load raw sequence data (FASTQ files) and metadata.
# Perform quality control checks on metadata (metadata.csv).
# Filter and clean metadata for downstream analyses.
# Save processed metadata and link it to sequence data.


# Load libraries
library(tidyverse)
library(data.table)
library(mia)


# 1. Import metadata
metadata <- fread("./data/raw/metadata.csv")

# 2. Add a column to flag mock samples
metadata <- metadata %>%
  mutate(mock = ifelse(grepl("mock", DNA_extract, ignore.case = TRUE), TRUE, FALSE))

# 3. Check for the presence of all libraries in the raw data folders
# Paths to the raw data directories
path_16S <- "./data/raw/16S_reads"
path_ITS <- "./data/raw/ITS_reads"

# Get list of FASTQ files in each directory
raw_16S_files <- list.files(path_16S, pattern = "\\.fq.gz$", full.names = FALSE)
raw_ITS_files <- list.files(path_ITS, pattern = "\\.fq.gz$", full.names = FALSE)

# Extract sample IDs from file names (assumes SampleID matches filenames without extensions)
raw_16S_sample_ids <- sub("\\.fq.gz$", "", raw_16S_files)
raw_ITS_sample_ids <- sub("\\.fq.gz$", "", raw_ITS_files)

# Identify missing libraries for 16S and ITS
missing_16S <- metadata$SampleID[!metadata$SampleID %in% raw_16S_sample_ids]
missing_ITS <- metadata$SampleID[!metadata$SampleID %in% raw_ITS_sample_ids]

# Warnings for missing libraries
if (length(missing_16S) > 0) {
  warning("The following 16S libraries are missing from 'data/raw/16S_reads/':")
  print(missing_16S)
} else {
  print("All 16S libraries are present in 'data/raw/16S_reads/'.")
}

if (length(missing_ITS) > 0) {
  warning("The following ITS libraries are missing from 'data/raw/ITS_reads/':")
  print(missing_ITS)
} else {
  print("All ITS libraries are present in 'data/raw/ITS_reads/'.")
}

# 4. Save the cleaned metadata
fwrite(metadata, "./data/intermediate/cleaned_metadata.csv")

# Print confirmation
print("Cleaned metadata with 'mock' column saved as 'cleaned_metadata.csv'")
