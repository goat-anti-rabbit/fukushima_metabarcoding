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


FUNS$plot_asv_length_abundance <- function(seqtab, META, mock_column = "mock", main = "ASV abundance vs. length") {
  require(scales)
  
  asv_seqs <- getSequences(seqtab)
  asv_lengths <- nchar(asv_seqs)
  asv_abundance <- colSums(seqtab)
  asv_prevalence <- colSums(seqtab > 0)
  
  # Identify which samples are mock
  mock_samples <- which(META[[mock_column]])
  mock_asvs <- colSums(seqtab[mock_samples, , drop = FALSE]) > 0
  
  # Assign colors
  col_vec <- rep(alpha("blue", 0.4), length(asv_seqs))          # default: in >1 sample
  col_vec[asv_prevalence == 1] <- alpha("red", 0.4)             # found in 1 sample only
  col_vec[mock_asvs] <- alpha("green", 0.6)                     # in mock samples
  
  # Plot
  plot(asv_lengths, asv_abundance,
       xlab = "ASV length (bp)",
       ylab = "Total abundance",
       main = main,
       log = "y",
       col = col_vec,
       pch = 20)
  
  legend("topleft", legend = c("Mock ASVs", "Single-sample ASVs", "Multi-sample ASVs"),
         col = c(alpha("green", 0.6), alpha("red", 0.4), alpha("blue", 0.4)),
         pch = 20, bty = "n")
}


# Helper: aggregate counts by soil_sample
FUNS$aggregate_counts_by_sample <- function(seqtab, meta) {
  meta$soil_sample <- factor(meta$soil_sample, levels = unique(meta$soil_sample))
  agg_counts <- rowsum(seqtab, group = meta$soil_sample)
  return(agg_counts)
}

# Helper: split TSE into mock and main + remove mock-only ASVs
FUNS$split_TSE_remove_mockASVs <- function(counts, meta) {
  is_mock <- meta$mock
  mock_counts <- counts[is_mock, , drop = FALSE]
  main_counts <- counts[!is_mock, , drop = FALSE]
  
  # Drop ASVs only found in mocks
  main_counts <- main_counts[, colSums(main_counts)>0, drop = FALSE]
  mock_counts <- mock_counts[, colSums(mock_counts)>0, drop = FALSE]  
  
  list(
    main = list(counts = main_counts, meta = meta[!is_mock, ]),
    mock = list(counts = mock_counts, meta = meta[is_mock, ])
  )
}

# ---- Create TSEs ----
FUNS$make_TSE <- function(counts, meta, prefix) {
  tse <- TreeSummarizedExperiment(
    assays = list(counts = t(counts)),
    colData = DataFrame(meta),
    rowData = DataFrame(sequence = colnames(counts))
  )
  # Rename ASVs and store sequences
  asv_seqs <- rownames(tse)
  asv_ids <- paste0(prefix, "_asv", seq_along(asv_seqs))
  rownames(tse) <- asv_ids
  rowData(tse)$sequence <- asv_seqs
  return(tse)
}




FUNS$getN     <- function(x) sum(getUniques(x))


FUNS$plot_read_retention <- function(track_df, pdf_file = NULL, title_prefix = "") {
  library(data.table)
  library(ggplot2)
  
  # Convert to data.table if needed
  track_dt <- as.data.table(track_df, keep.rownames = "sample")
  
  # Melt to long format
  track_long <- melt(track_dt, id.vars = "sample", variable.name = "step", value.name = "reads")
  
  # Normalize to percent of input
  track_long[, input_reads := reads[step == "input"], by = sample]
  track_long[, pct := 100 * reads / input_reads]
  
  # Plot 1: Absolute read counts with labels
  p1 <- ggplot(track_long, aes(x = step, y = reads, group = sample)) +
    geom_line(alpha = 0.5, color = "steelblue") +
    geom_point(size = 0.8, color = "steelblue") +
    geom_text(
      data = track_long[step == "nonchim"],
      aes(label = sample),
      hjust = 0, nudge_x = 0.1, size = 2
    ) +
    labs(title = paste(title_prefix, "Read Retention Across Pipeline (Absolute)"),
         x = "Pipeline Step", y = "Read Count") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(5.5, 60, 5.5, 5.5)
    )
  
  # Plot 2: Percentage of input
  p2 <- ggplot(track_long, aes(x = step, y = pct, group = sample)) +
    geom_line(alpha = 0.5, color = "steelblue") +
    geom_point(size = 0.8, color = "steelblue") +
    scale_y_continuous(limits = c(0, 110)) +
    labs(title = paste(title_prefix, "Read Retention Across Pipeline (% of Input)"),
         x = "Pipeline Step", y = "% of Input Reads") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save to PDF if filename provided
  if (!is.null(pdf_file)) {
    pdf(pdf_file, width = 12, height = 8)
    print(p1)
    print(p2)
    dev.off()
  }
  
  return(list(abs = p1, pct = p2))
}


# Helper to calculate % identification at taxonomic ranks
summarize_taxonomy <- function(tse, ranks = c("order", "family", "genus", "species")) {
  tax <- as.data.frame(rowData(tse))
  
  if (!all(ranks %in% colnames(tax))) {
    stop("Not all requested ranks are present in rowData.")
  }
  
  n_total <- nrow(tax)
  summary <- sapply(ranks, function(rank) {
    sum(!is.na(tax[[rank]]) & tax[[rank]] != "") / n_total * 100
  })
  return(round(summary, 1))
}

# Helper to calculate % identification at taxonomic ranks
FUNS$summarize_taxonomy <- function(tse, ranks = c("Order", "Family", "Genus", "Species")) {
  tax <- as.data.frame(rowData(tse))
  
  if (!all(ranks %in% colnames(tax))) {
    stop("Not all requested ranks are present in rowData.")
  }
  
  n_total <- nrow(tax)
  summary <- sapply(ranks, function(rank) {
    sum(!is.na(tax[[rank]]) & tax[[rank]] != "") / n_total * 100
  })
  return(round(summary, 1))
}


### Function that compares taxonomy assignments in two similar TSE objects in which 
### taxonomy assignment has been completed with different settings
### Selects n sequences at tax_level that have an assignment in one, but not in the other TSE object
### Outputs them to a fasta file such that they can be investigated with e.g. BLASTn


FUNS$sample_taxonomic_gain_ASVs <- function(strict_TSE, lenient_TSE, tax_level = "Genus", n = 20, outfile = "sampled_gain.fasta") {
  stopifnot(tax_level %in% colnames(rowData(strict_TSE)))
  
  strict_tax  <- rowData(strict_TSE)[[tax_level]]
  lenient_tax <- rowData(lenient_TSE)[[tax_level]]
  
  # Find ASVs with NA in strict but not in lenient taxonomy
  gained_tax <- is.na(strict_tax) & !is.na(lenient_tax)
  gained_ids <- rownames(strict_TSE)[gained_tax]
  gained_seqs <- rowData(lenient_TSE)$sequence[gained_tax]
  
  # Spread across taxonomic space using kmer distance clustering
  kmers <- Biostrings::oligonucleotideFrequency(Biostrings::DNAStringSet(gained_seqs), width = 5)
  dists <- dist(kmers)
  clusters <- cutree(hclust(dists), k = min(n, length(gained_seqs)))
  #print(table(clusters))
  
  # Sample one ASV per cluster (or fewer if n > total)
  sampled_indices   <- unlist(lapply(unique(clusters), function(cl) {
    cluster_indices <- which(clusters == cl)
    sample(cluster_indices, 1)
  }))
  
  sampled_seqs <- rowData(lenient_TSE)$sequence[sampled_indices]
  
  # Make informative headers
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  header_strings <- vapply(sampled_indices, function(id) {
    taxonomy <- unlist(rowData(lenient_TSE)[id, tax_levels])
    taxonomy <- taxonomy[!is.na(taxonomy)]
    paste(unlist(rownames(rowData(lenient_TSE))[id]), paste(taxonomy, collapse = "; "), sep = " | ")
  }, character(1))
  
  names(sampled_seqs) <- header_strings
  
  # Write to FASTA
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(sampled_seqs),
    filepath = outfile
  )
  
  message(length(sampled_indices), " sequences written to ", outfile)
  invisible(sampled_indices)
}

FUNS$sample_taxonomic_gain_ASVs <- function(strict_TSE, lenient_TSE, tax_level = "Genus", n = 20, outfile = "sampled_gain.fasta") {
  # Safety checks
  stopifnot(tax_level %in% colnames(rowData(strict_TSE)))
  stopifnot(tax_level %in% colnames(rowData(lenient_TSE)))
  
  # Define hierarchy
  tax_hierarchy <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  target_idx <- match(tax_level, tax_hierarchy)
  if (is.na(target_idx)) stop("Invalid taxonomic level.")
  
  higher_levels <- tax_hierarchy[1:(target_idx - 1)]
  
  # Extract taxonomy columns
  strict_tax_all  <- rowData(strict_TSE)[, c(higher_levels, tax_level), drop = FALSE]
  lenient_tax_all <- rowData(lenient_TSE)[, c(higher_levels, tax_level), drop = FALSE]
  
  # Keep ASVs that match higher levels, and gain assignment at the target level
  keep <- apply(strict_tax_all[, higher_levels, drop=FALSE], 1, function(x) all(!is.na(x))) &
    apply(lenient_tax_all[, higher_levels, drop=FALSE], 1, function(x) all(!is.na(x))) &
    is.na(strict_tax_all[[tax_level]]) &
    !is.na(lenient_tax_all[[tax_level]])
  
  if (sum(keep) == 0) stop("No ASVs meet the gain-at-this-level criterion.")
  
  gained_ids <- rownames(strict_TSE)[keep]
  gained_seqs <- rowData(lenient_TSE)$sequence[keep]
  
  # Spread across taxonomic space using k-mer distance clustering
  kmers <- Biostrings::oligonucleotideFrequency(Biostrings::DNAStringSet(gained_seqs), width = 5)
  dists <- dist(kmers)
  clusters <- cutree(hclust(dists), k = min(n, length(gained_seqs)))
  
  sampled_indices <- unlist(lapply(unique(clusters), function(cl) {
    cluster_indices <- which(clusters == cl)
    sample(cluster_indices, 1)
  }))
  
  sampled_seqs <- rowData(lenient_TSE)$sequence[keep][sampled_indices]
  
  # Build informative headers
  tax_levels <- tax_hierarchy
  header_strings <- vapply(sampled_indices, function(idx) {
    tax_row <- rowData(lenient_TSE)[keep, , drop=FALSE][idx, tax_levels]
    tax_row_vec <- as.vector(as.matrix(tax_row))
    tax_row_vec <- tax_row_vec[!is.na(tax_row_vec)]
    asv_id <- rownames(strict_TSE)[keep][idx]
    paste0(asv_id, " | ", paste(tax_row_vec, collapse = "; "))
  }, character(1))
  
  names(sampled_seqs) <- header_strings
  
  # Write FASTA
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(sampled_seqs),
    filepath = outfile
  )
  
  message(length(sampled_seqs), " sequences written to ", outfile)
  invisible(sampled_indices)
}

