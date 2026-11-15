# @ Rik Verdonck 20261026
# Purpose: collect the read files that go into denoising: one pair of fastq files per library (DNA extract x marker).
#
# Ten soil samples were extracted twice (technical replicates, e.g. A2T and A2T1), and every DNA extract was sequenced
# as a separate library. 
# Every library is processed separately, the replicate extracts included.
# The main analyses use one extract per sample (column primary in data/metadata.csv): the original extract
# (A2T, not A2T1; F3T only has F3T1), so that all samples have the same sampling effort.
# Libraries of the excluded samples (column excluded) are left out. The files are linked, not copied.
#
# Input:  data/metadata.csv, data/reads_<marker>/delivered/ (the files named in the lib_<marker>_... columns)
# Output: data/reads_<marker>/input/

source("src/config.R")

meta <- data.table::fread("data/metadata.csv", na.strings = "", data.table = FALSE)
meta <- meta[!meta$excluded, ]

for (marker in c("16S", "ITS")) {
  path.delivered <- PATH$reads(marker, "delivered")
  path.input     <- PATH$reads(marker, "input")
  dir.create(path.input, recursive = TRUE, showWarnings = FALSE)

  # Start from an empty folder, so that no file from an earlier run remains
  unlink(list.files(path.input, pattern = "\\.fq\\.gz$", full.names = TRUE))

  # The read files of this marker (empty when an extract has no library for this marker).
  files <- c(meta[[paste0("lib_", marker, "_forward")]], meta[[paste0("lib_", marker, "_reverse")]])
  files <- files[!is.na(files)]

  for (f in files) {
    file.symlink(normalizePath(file.path(path.delivered, f)), file.path(path.input, f))
  }

  cat(marker, ":", length(list.files(path.input)), "read files linked\n")
}

# Nothing of this script is needed further on (the next scripts read data/metadata.csv themselves)
rm(meta, marker, path.delivered, path.input, files, f)
