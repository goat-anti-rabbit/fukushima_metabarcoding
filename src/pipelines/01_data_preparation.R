# @ Rik Verdonck 20260313

# Purpose: Handles raw data and metadata import, initial checks, and setup.

# Load raw sequence data (FASTQ files) and metadata.
# Perform quality control checks on metadata (metadata.csv).
# Filter and clean metadata for downstream analyses.
# Save processed metadata and link it to sequence data.
# We already did some concatenations and kicked out some samples before getting started.
# See the bash script concatenate_split_lanes.sh in the src/bash_scripts directory for details
# This was based on previous pilot analyses. 
# See a very quick analysis at the end of this script to justify outlier removal



# Load libraries
library("tidyverse")
library("data.table")
library("ShortRead")

# 0. Set working directory and read in necessary functions. 
source("src/functions/functions.R")

# 1. Import metadata
metadata <- fread("./data/metadata.csv")

# 2. Add a column to flag mock samples
metadata <- metadata %>%
  mutate(mock = ifelse(grepl("mock", DNA_extract, ignore.case = TRUE), TRUE, FALSE))


# 3. Check for the presence of all libraries in the raw data folders
# Paths to the raw data directories
path_16S <- "./data/reads_16S/raw"
path_ITS <- "./data/reads_ITS/raw"

# Get list of FASTQ files in each directory
raw_16S_files <- list.files(path_16S, pattern = "\\.fq.gz$", full.names = FALSE)
raw_ITS_files <- list.files(path_ITS, pattern = "\\.fq.gz$", full.names = FALSE)

# Extract sample IDs from file names (assumes SampleID matches filenames without extensions)
raw_16S_sample_ids <- sub("\\.fq.gz$", "", raw_16S_files)
raw_ITS_sample_ids <- sub("\\.fq.gz$", "", raw_ITS_files)

# Which fq files originate from concatenated fq files?
# Get list of FASTQ files in each directory
raw_16S_concat <- list.files("./data/reads_16S/unconcatenated", pattern = "\\.fq.gz$", full.names = FALSE)
raw_ITS_concat <- list.files("./data/reads_16S/unconcatenated", pattern = "\\.fq.gz$", full.names = FALSE)



# Identify missing libraries for 16S and ITS
missing_from_metadata_16S <- setdiff(raw_16S_files,c(metadata$lib_16S_forward,metadata$lib_16S_reverse))
missing_from_metadata_ITS <- setdiff(raw_ITS_files,c(metadata$lib_ITS_forward,metadata$lib_ITS_reverse))

# Identify missing libraries for 16S and ITS
missing_from_files_16S <- setdiff(c(metadata$lib_16S_forward,metadata$lib_16S_reverse),raw_16S_files)
missing_from_files_ITS <- setdiff(c(metadata$lib_ITS_forward,metadata$lib_ITS_reverse),raw_ITS_files)

# We define a metadata frame for 16S
# We filter the metadata to include only rows where either the forward or reverse library matches 
# a file in the raw 16S files.
# We also remove ITS-specific information
# We add a column that flags whether the library is from a concatenated file or not, 
# based on the presence of the library in the list of concatenated files.
# We also add a column with read counts

meta_16S <- metadata %>%
  filter(lib_16S_forward %in% raw_16S_files | lib_16S_reverse %in% raw_16S_files) %>%
  filter(soil_sample != "A5T" & soil_sample != "T3T") %>% # Remove outliers A5T and T3T)
  mutate(concatenated = ifelse(lib_16S_forward %in% raw_16S_concat, TRUE, FALSE)) %>%
  select(-contains("ITS")) %>%
  dplyr::rename(lib_forward = lib_16S_forward, lib_reverse = lib_16S_reverse) %>%
  dplyr::rename(primer_fw = primer_16S_forward, primer_rv = primer_16S_reverse) %>%
  rowwise() %>%
  mutate(Nreads = ShortRead::countFastq(path_16S, lib_reverse)$records) %>%
  mutate(Nnucl  = ShortRead::countFastq(path_16S, lib_reverse)$nucleotides) %>%
  ungroup()


### Now we do the same thing for ITS:

meta_ITS <- metadata %>%
  filter(lib_ITS_forward %in% raw_ITS_files | lib_ITS_reverse %in% raw_ITS_files) %>%
  mutate(concatenated = ifelse(lib_ITS_forward %in% raw_ITS_concat, TRUE, FALSE)) %>%
  select(-contains("16S")) %>%
  dplyr::rename(lib_forward = lib_ITS_forward, lib_reverse = lib_ITS_reverse) %>%
  dplyr::rename(primer_fw = primer_ITS_forward, primer_rv = primer_ITS_reverse) %>%
  rowwise() %>%
  mutate(Nreads = ShortRead::countFastq(path_ITS, lib_reverse)$records) %>%
  mutate(Nnucl  = ShortRead::countFastq(path_ITS, lib_reverse)$nucleotides) %>%
  ungroup()


# 4. Save the cleaned metadata
fwrite(meta_16S, "./data/metadata_16S.csv",na="NA")
fwrite(meta_ITS, "./data/metadata_ITS.csv",na="NA")


# Print confirmation
print("Cleaned metadata, added columns for 'mock', and read counts for both primer pairs.")



### Removal of outliers: 
# I find little evidence for outliers at the level of ITS, but for 16S there are two samples 
# that are clearly outliers: A5T and T3T 
# This becomes clear with a quick pairs plot presenting a sample of alpha diversity measures, 
# where A5T is colored red, and T3T orange. A5T and T3T must be dominated by a couple of 
# very abundant OTU’s, while otherwise being composed of many rare taxa.
# We will also find out that these samples have an extreme influence on ordination 
# of beta diversity measures, and I think it would be best to exclude them altogether, 
# also because they don’t have particularly interesting values for 137Cs.
# And coincidentally, we don't have libraries for these samples for ITS, so perhaps the extracts were really poor. 
# So I moved these library files to the "excluded" subdirectory in the data directory
# And removed the entries from the metadata as well.
# You can find the entire rationale in earlier versions of this workflow. 
# For reasons of wanting to keep this concise, it's not all explained here. 







