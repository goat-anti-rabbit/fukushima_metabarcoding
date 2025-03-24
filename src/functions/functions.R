FUNS <- list()

# Count the reads in a fq.gz file based on filename
FUNS$count_reads <- function(file_path) {
  if (file.exists(file_path)) {
    fq <- ShortRead::FastqStreamer(file_path, n = 1e6)
    count <- 0
    while (length(reads <- ShortRead::yield(fq)) > 0) {
      count <- count + length(reads)
    }
    close(fq)
    return(count)
  } else {
    return(NA)
  }
}

FUNS$get.sample.name      <- function(fname) strsplit(basename(fname), "_")[[1]][1]

# Helper function to save paired quality profiles in one PDF
FUNS$save_paired_quality_profiles <- function(forward_files, reverse_files, library_name) {
  # Ensure the lists of forward and reverse files are the same length
  if (length(forward_files) != length(reverse_files)) {
    stop("Mismatch in the number of forward and reverse files.")
  }
  
  # Create a PDF to save plots
  output_filename <- file.path(output_dir, paste0(library_name, "_paired_quality_profiles.pdf"))
  pdf(output_filename, width = 8, height = 6)  # Set dimensions
  
  for (i in seq_along(forward_files)) {
    # Extract sample name from filename
    sample_name <- basename(forward_files[i])
    sample_name <- sub("_R[12]_001.fastq", "", sample_name)
    
    # Generate plots for forward and reverse reads
    p1 <- plotQualityProfile(forward_files[i]) + ggtitle(paste0(library_name, " - ", sample_name, " - Forward"))
    p2 <- plotQualityProfile(reverse_files[i]) + ggtitle(paste0(library_name, " - ", sample_name, " - Reverse"))
    
    # Arrange both plots on a single page
    grid.arrange(p1, p2, ncol = 1)
  }
  
  dev.off()
  cat("Paired quality profiles saved to", output_filename, "\n")
}

# Takes a primer, and returns it alongside its reverse, reverse complement and complement
FUNS$allOrients <- function(primer) {
  require("Biostrings")
  dna <- DNAString(primer)  
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  
}

# Takes a single (primer) sequence and counts the number of occurences in fn (filepath to read library in fq.gz)
FUNS$primerHits <- function(primer, fn) {
  require("Biostrings")
  require("ShortRead")
  # Counts number of reads in which the primer is found
  nhits <- Biostrings::vcountPattern(primer, ShortRead::sread(ShortRead::readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Takes a dataset (XITS or X16S), the name the reads fw, the name of the reads rev
# This function generates a table in which we count the number of reads that contain a primer in any direction. 
# The column names are a bit of a mess:
# Upper case FW or RV: whether we are looking for the forward or the reverse primers
# Lower case fw or rv: whether we are looking into the forward or reverse reads (you can for example find forward primers in reverse reads in case of readthrough)
# The direction of the primer. Whether it's oriented the way it should, reversed, complemented or reverse complemented. 

FUNS$generate_primer_summary <- function(dataset, fnFs_col, fnRs_col, multithread = 4) 
{
  require("doParallel")
  cl <- makeCluster(multithread)
  registerDoParallel(cl)
  # Perform primer hits calculation in parallel
  primer_summary <- foreach(i = 1:length(dataset[[fnFs_col]]), .combine = rbind, .export=c("FUNS")) %dopar% {
    c( FW.fw = sapply(dataset$FWD.orients, FUNS$primerHits, fn = dataset[[fnFs_col]][[i]]),
       FW.rv = sapply(dataset$FWD.orients, FUNS$primerHits, fn = dataset[[fnRs_col]][[i]]),
       RV.fw = sapply(dataset$REV.orients, FUNS$primerHits, fn = dataset[[fnFs_col]][[i]]),
       RV.rv = sapply(dataset$REV.orients, FUNS$primerHits, fn = dataset[[fnRs_col]][[i]])
    )
  }
  
  rownames(primer_summary) <- dataset$sample.names
  stopCluster(cl)
  return(primer_summary)
}

# Learn DADA2 error models from a balanced subset of all filtered libraries.
# Instead of letting DADA2 use reads from only the largest few files (which
# happens when nbases is hit early), this function samples a fixed number of
# reads from each input FASTQ file, combines them, and estimates error rates.
#
# Args:
#   filt_files           Character vector with paths to filtered FASTQ files
#   n_reads_per_sample   Integer: how many reads to sample per file (default = 10000)
#   nbases               Integer: max number of bases to use in training (default = 1e8)
#   multithread          Logical: use multiple cores (default = TRUE)
#   verbose              Logical: print progress messages (default = TRUE)
#
# Returns:
#   An error model object as returned by dada2::learnErrors()
#
# Example:
#   errF <- learnErrorsBalanced(filtFs, n_reads_per_sample = 10000, multithread = TRUE)
# ------------------------------------------------------------------------------


FUNS$learnErrorsBalanced <- function(filt_files, 
                                     n_reads_per_sample = 10000, 
                                     nbases = 1e8, 
                                     multithread = TRUE, 
                                     verbose = TRUE,
                                     MAX_CONSIST = 10) {
  
  library(ShortRead)
  library(dada2)
  
  tmp <- tempfile(fileext = ".fastq.gz")
  
  for (f in filt_files) {
    fq <- readFastq(f)
    n <- min(n_reads_per_sample, length(fq))
    fq_sub <- fq[1:n]
    writeFastq(fq_sub, tmp, mode = "a", compress = TRUE)
  }
  
  err <- learnErrors(tmp, multithread = multithread, nbases = nbases, 
                     verbose = verbose, MAX_CONSIST = MAX_CONSIST)
  return(err)
}



