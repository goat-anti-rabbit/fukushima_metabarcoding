# @ Rik Verdonck 20261026
# Purpose: functions used by the scripts in src/. Load with source("src/functions/functions.R").
# All functions are stored in the list FUNS.

FUNS <- list()


# Number of reads in a dada, derep or mergePairs object (as in the DADA2 tutorial)
FUNS$getN <- function(x) sum(dada2::getUniques(x))


# Runs an external program (cutadapt, vsearch) and stops with an error when it fails, 
# so that a failed call can never go unnoticed. 
# The files in `outputs` are deleted before the call (so that an old file cannot be taken for
# a new one) and must exist afterwards. stdout and stderr are passed on to system2 (a log file, "" for the console,
# or FALSE to discard).
FUNS$run_tool <- function(command, args, outputs = character(0), stdout = FALSE, stderr = FALSE) {
  unlink(outputs)
  status <- system2(command, shQuote(args), stdout = stdout, stderr = stderr)
  if (!identical(as.integer(status), 0L)) {
    stop(command, " failed with exit status ", status, if (is.character(stdout) && nzchar(stdout)) paste0("; see ", stdout))
  }
  missing <- outputs[!file.exists(outputs)]
  if (length(missing) > 0) stop(command, " did not write ", paste(missing, collapse = ", "))
  invisible(TRUE)
}


# Read totals in the colData of a TSE. "reads" is the current number of reads of each library (the column sums of
# the counts) and must be updated whenever ASVs are removed. "reads_denoised" is the number of reads after
# denoising, the last step of the read tracking (this column was called "final" in objects made before 4 Nov 2026).
FUNS$update_reads <- function(tse) {
  cd <- SummarizedExperiment::colData(tse)
  if ("final" %in% names(cd)) {
    if ("reads_denoised" %in% names(cd)) stop("colData has both 'final' and 'reads_denoised'.")
    names(cd)[names(cd) == "final"] <- "reads_denoised"
  }
  cd$reads <- colSums(SummarizedExperiment::assay(tse, "counts"))
  SummarizedExperiment::colData(tse) <- cd
  tse
}


# Turns a DADA2 ASV table (rows = libraries, columns = ASV sequences) into a TreeSummarizedExperiment:
#   assay "counts": rows = ASVs, columns = libraries
#   colData:        the metadata of each library (one row of data/metadata.csv), plus the read tracking columns
#                   ("final" becomes "reads_denoised"), plus "reads" (see FUNS$update_reads)
#   rowData:        the length of each ASV; the sequences themselves are stored in referenceSeq
# The ASVs get short names (id_prefix + number, most abundant first), e.g. 16S_ASV00001.
FUNS$make_tse <- function(asv_table, meta, track, id_prefix = "ASV") {
  if (!identical(rownames(asv_table), meta$DNA_extract) || !identical(meta$DNA_extract, track$DNA_extract)) {
    stop("ASV table, metadata and read tracking must have the same libraries in the same order.")
  }

  sequences <- colnames(asv_table)
  ids       <- sprintf("%s%05d", id_prefix, seq_along(sequences))

  counts <- t(asv_table)
  rownames(counts) <- ids

  coldata <- cbind(meta, track[, setdiff(names(track), names(meta)), drop = FALSE])
  rownames(coldata) <- meta$DNA_extract

  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays       = list(counts = counts),
    colData      = S4Vectors::DataFrame(coldata),
    rowData      = S4Vectors::DataFrame(length = nchar(sequences), row.names = ids),
    referenceSeq = Biostrings::DNAStringSet(setNames(sequences, ids))
  )
  FUNS$update_reads(tse)
}


# Compares the ASVs of a mock community with the reference sequences of its members (vsearch, global alignment,
# both strands), in both directions. Returns a list of four tables:
#   hits:       ASV -> reference: for every ASV all equally good references (ties are all kept), one row each
#   asvs:       one row per ASV: its reads per library and its best reference(s), ties separated by " | "
#   ref_hits:   reference -> ASV: for every reference all equally good ASVs, one row each
#   references: one row per reference sequence, with its best ASV(s), ties separated by " | "
# hits and ref_hits give the alignment details and match class of every single hit. Tied hits have the same identity
# but can differ in coverage; in the one-row summaries (asvs, references), asv_coverage and match are those of the
# most complete alignment among the ties. Decide per hit (from hits or ref_hits) whether something is recovered.
# identity: % identity over the aligned part, as reported by vsearch (terminal gaps not counted);
# asv_coverage, reference_coverage: % of the ASV / of the reference covered by the alignment. The Sanger references of
# the fungal isolates can be shorter than the amplicon, so a high identity can concern only part of the ASV.
# strand: "+" or "-", the strand of the query sequence that matched. See FUNS$match_class for the classes.
FUNS$mock_check <- function(tse_mock, reference_fasta) {
  asv_fasta <- tempfile(fileext = ".fasta")
  ref_fasta <- tempfile(fileext = ".fasta")
  on.exit(unlink(c(asv_fasta, ref_fasta)), add = TRUE)
  Biostrings::writeXStringSet(TreeSummarizedExperiment::referenceSeq(tse_mock), asv_fasta)
  ref <- Biostrings::readDNAStringSet(reference_fasta)
  names(ref) <- sub("\\s.*$", "", trimws(names(ref)))   # vsearch also cuts names at the first space
  Biostrings::writeXStringSet(ref, ref_fasta)

  # All best hits (equal identity) of each query sequence in the database; queries without a hit of at least 70 % are absent
  best_hits <- function(query_fasta, db_fasta) {
    hits_file <- tempfile(fileext = ".tsv")
    on.exit(unlink(hits_file))
    FUNS$run_tool("vsearch", c("--usearch_global", query_fasta, "--db", db_fasta, "--id", 0.7, "--strand", "both",
                               "--maxaccepts", 0, "--maxrejects", 0, "--top_hits_only", "--userout", hits_file,
                               "--userfields", "query+target+id+alnlen+mism+gaps+qcov+tcov+qstrand", "--quiet"),
                  outputs = hits_file)
    columns <- c("query", "target", "identity", "alignment_length", "mismatches", "gaps", "query_coverage",
                 "target_coverage", "strand")
    if (file.size(hits_file) == 0) return(setNames(data.frame(matrix(nrow = 0, ncol = 9)), columns))
    read.table(hits_file, sep = "\t", col.names = columns, quote = "", comment.char = "", colClasses = c(strand = "character"))
  }
  details <- c("identity", "alignment_length", "mismatches", "gaps")

  # ASV -> reference
  h <- best_hits(asv_fasta, ref_fasta)
  hits <- data.frame(ASV = h$query, reference = h$target, h[, details], asv_coverage = h$query_coverage,
                     reference_coverage = h$target_coverage, strand = h$strand)
  hits$match <- FUNS$match_class(hits$identity, hits$mismatches, hits$gaps, hits$asv_coverage)
  hits <- hits[order(hits$ASV, -hits$asv_coverage), ]           # per ASV, the most complete alignment first

  # reference -> ASV (shows the closest ASV also for references that are nobody's best reference)
  h <- best_hits(ref_fasta, asv_fasta)
  ref_hits <- data.frame(reference = h$query, ASV = h$target, h[, details], asv_coverage = h$target_coverage,
                         reference_coverage = h$query_coverage, strand = h$strand)
  ref_hits$match <- FUNS$match_class(ref_hits$identity, ref_hits$mismatches, ref_hits$gaps, ref_hits$asv_coverage)
  ref_hits <- ref_hits[order(ref_hits$reference, -ref_hits$asv_coverage), ]

  # One row per query sequence in `ids`: its best hits joined by " | " (or "none"), their number, and identity,
  # asv_coverage and match of the first (= most complete) of them
  summarise <- function(h, query, target, ids) {
    first <- h[match(ids, h[[query]]), ]
    best  <- vapply(ids, function(x) paste(h[[target]][h[[query]] == x], collapse = " | "), "", USE.NAMES = FALSE)
    data.frame(best = ifelse(best == "", "none", best),
               n_best = vapply(ids, function(x) sum(h[[query]] == x), 0L, USE.NAMES = FALSE),
               identity = first$identity, asv_coverage = first$asv_coverage,
               match = ifelse(is.na(first$identity), "no hit", first$match))
  }

  counts <- SummarizedExperiment::assay(tse_mock, "counts")
  s <- summarise(hits, "ASV", "reference", rownames(counts))
  asvs <- data.frame(ASV = rownames(counts), counts, reference = s$best, s[, -1], check.names = FALSE, row.names = NULL)
  s <- summarise(ref_hits, "reference", "ASV", names(ref))
  references <- data.frame(reference = names(ref), ASV = s$best, s[, -1])

  list(hits = hits, asvs = asvs, ref_hits = ref_hits, references = references)
}


# Match classes for FUNS$mock_check, by explicit comparisons:
#   "identical"        the whole ASV is aligned, without mismatches or gaps
#   "99 to < 100 %", "97 to < 99 %"   the whole ASV is aligned, with this identity
#   "partial overlap"  at least 97 % identity, but only over part of the ASV
#   "< 97 %"           less than 97 % identity
#   "no hit"           nothing with at least 70 % identity
FUNS$match_class <- function(identity, mismatches, gaps, asv_coverage) {
  dplyr::case_when(
    is.na(identity)         ~ "no hit",
    identity < 97           ~ "< 97 %",
    asv_coverage < 100      ~ "partial overlap",
    mismatches + gaps == 0  ~ "identical",
    identity >= 99          ~ "99 to < 100 %",
    identity >= 97          ~ "97 to < 99 %",
    TRUE                    ~ "< 97 %"
  )
}


# For every soil sample that was extracted twice: how similar are its two extracts?
# tse: object with all extracts (columns), with colData columns soil_sample, primary and sampling_site.
# Per sample: ASVs in each extract, the share of ASVs (and of reads) shared by both extracts, the Bray-Curtis
# dissimilarity between the two extracts, and for comparison the median Bray-Curtis dissimilarity between the
# original extract and the other samples of the same site.
FUNS$compare_extracts <- function(tse) {
  counts <- SummarizedExperiment::assay(tse, "counts")
  cd     <- as.data.frame(SummarizedExperiment::colData(tse))
  bray   <- as.matrix(vegan::vegdist(t(counts) / colSums(counts), "bray"))

  samples_with_two <- names(which(table(cd$soil_sample) == 2))
  result <- NULL
  for (s in samples_with_two) {
    original  <- cd$DNA_extract[cd$soil_sample == s & cd$primary]
    replicate <- cd$DNA_extract[cd$soil_sample == s & !cd$primary]
    site      <- cd$sampling_site[cd$DNA_extract == original]
    others    <- cd$DNA_extract[cd$sampling_site == site & cd$primary & cd$soil_sample != s]
    in_original  <- counts[, original] > 0
    in_replicate <- counts[, replicate] > 0
    shared       <- in_original & in_replicate
    result <- rbind(result, data.frame(
      soil_sample                     = s,
      site                            = site,
      asvs_original                   = sum(in_original),
      asvs_replicate                  = sum(in_replicate),
      shared_asvs_pct                 = round(100 * sum(shared) / sum(in_original | in_replicate), 1),
      reads_in_shared_asvs_pct        = round(100 * sum(counts[shared, c(original, replicate)]) / sum(counts[, c(original, replicate)]), 1),
      bray_between_extracts           = round(bray[original, replicate], 3),
      bray_to_other_samples_same_site = round(median(bray[original, others]), 3)))
  }
  result
}


# Compares several versions of a count table of the same soil samples (rows = samples in the same order, columns = ASVs).
# All tables are rarefied to the same depth for richness. Per table: the site effect on Bray-Curtis dissimilarities
# (PERMANOVA) and the agreement of its Bray-Curtis dissimilarities with those of the first table (Mantel).
FUNS$compare_tables <- function(tables, site) {
  depth    <- min(sapply(tables, function(x) min(rowSums(x))))
  richness <- data.frame(site = site, sapply(tables, function(x) vegan::rarefy(x, depth)))
  bray     <- lapply(tables, function(x) vegan::vegdist(x / rowSums(x), "bray"))
  stats <- NULL
  for (n in names(tables)) {
    permanova <- vegan::adonis2(bray[[n]] ~ site, permutations = 999)
    stats <- rbind(stats, data.frame(
      table               = n,
      site_R2             = round(permanova$R2[1], 3),
      site_p              = permanova$`Pr(>F)`[1],
      mantel_r_with_first = round(vegan::mantel(bray[[1]], bray[[n]], permutations = 999)$statistic, 3)))
  }
  list(depth = depth, richness = richness, stats = stats)
}


# Saves the quality profiles of forward and reverse reads in one pdf, one page per library (forward above, reverse below)
FUNS$save_paired_quality_profiles <- function(forward_files, reverse_files, library_names, output_file) {
  if (length(forward_files) != length(reverse_files)) {
    stop("Mismatch in the number of forward and reverse files.")
  }

  pdf(output_file, width = 8, height = 6)
  for (i in seq_along(forward_files)) {
    p1 <- dada2::plotQualityProfile(forward_files[i]) + ggplot2::ggtitle(paste(library_names[i], "- forward"))
    p2 <- dada2::plotQualityProfile(reverse_files[i]) + ggplot2::ggtitle(paste(library_names[i], "- reverse"))
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  }
  dev.off()
}


# Plots the share of read pairs left after each step, one line per library (mocks in red)
# track: data frame with a column per step (read pairs left), the first step being "input"
FUNS$plot_read_tracking <- function(track, steps, is_mock, main, output_file) {
  kept <- 100 * track[, steps] / track$input
  pdf(output_file, width = 8, height = 6)
  par(mar = c(9, 5, 3, 1))
  matplot(t(kept), type = "l", lty = 1, col = ifelse(is_mock, "red", "grey40"),
          xaxt = "n", ylim = c(0, 100), xlab = "", ylab = "% of read pairs left", main = main)
  axis(1, at = seq_along(steps), labels = steps, las = 2)
  legend("bottomleft", legend = c("soil libraries", "mocks"), col = c("grey40", "red"), lty = 1, bty = "n")
  dev.off()
}


# Plots the length distribution of the ASVs: number of ASVs and number of reads per length
FUNS$plot_asv_lengths <- function(asv_table, main, output_file) {
  len <- nchar(colnames(asv_table))
  pdf(output_file, width = 8, height = 8)
  par(mfrow = c(2, 1))
  barplot(table(len), main = paste(main, "- number of ASVs per length"), xlab = "length (nt)", ylab = "ASVs")
  barplot(tapply(colSums(asv_table), len, sum), main = paste(main, "- number of reads per length"), xlab = "length (nt)", ylab = "reads")
  dev.off()
}


# Learn DADA2 error models from a balanced subset of all filtered libraries.
# Instead of letting DADA2 use reads from only the first few files (which
# happens when nbases is hit early), this function takes a fixed number of
# reads (the first ones) from each input FASTQ file, combines them, and estimates error rates.
#
# Args:
#   filt_files           Character vector with paths to filtered FASTQ files
#   n_reads_per_sample   Integer: how many reads to take per file (default = 10000)
#   nbases               Integer: max number of bases to use in training (default = 1e8)
#   multithread          Logical: use multiple cores (default = TRUE)
#   verbose              Logical: print progress messages (default = TRUE)
#   MAX_CONSIST          Integer: maximum number of rounds of error learning (default = 10)
#
# Returns:
#   An error model object as returned by dada2::learnErrors()
#
# Example:
#   errF <- FUNS$learnErrorsBalanced(filtFs, n_reads_per_sample = 10000, multithread = TRUE)
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


# How much the error rates of a DADA2 error model still changed in each round ("selfConsist step") of learnErrors.
# DADA2 stops early only when a new error matrix is identical to an earlier one; otherwise it stops after MAX_CONSIST
# rounds with "Self-consistency loop terminated before convergence".
# err: output of dada2::learnErrors. err$err_in holds the input matrix of each round (the first one is the initial,
# maximal estimate) and err$err_out the output of the last round, so round k changes err_in[[k]] into err_in[[k + 1]]
# (or into err_out, for the last round).
# Per round, for the substitution rates only (A2C, A2G, ...; the rates of correct base calls are left out):
#   max_abs_change:  the largest absolute change of any substitution rate
#   n_rates_above:   the number of rates above 0.001 before and after the round
#   max_rel_change, median_rel_change: largest and median relative change of these rates
FUNS$error_convergence <- function(err) {
  matrices <- c(err$err_in, list(err$err_out))
  substitution <- !(rownames(err$err_out) %in% c("A2A", "C2C", "G2G", "T2T"))
  result <- NULL
  for (k in seq_along(err$err_in)) {
    before <- matrices[[k]][substitution, ]
    after  <- matrices[[k + 1]][substitution, ]
    above  <- before > 0.001 & after > 0.001
    relative <- abs(after - before)[above] / before[above]
    result <- rbind(result, data.frame(
      round             = k,
      max_abs_change    = signif(max(abs(after - before)), 3),
      n_rates_above     = sum(above),
      max_rel_change    = if (any(above)) signif(max(relative), 3) else NA,      # NA when no rate is above 0.001
      median_rel_change = if (any(above)) signif(median(relative), 3) else NA))
  }
  result
}


# Assigns taxonomy with dada2::assignTaxonomy (naive Bayesian classifier, 100 bootstraps) against one reference file.
# Returns a data frame with one row per sequence: a column per rank with the name (e.g. family) and a column per rank
# with its bootstrap confidence, 0-100 (e.g. family_confidence). Names with a confidence below min_boot are NA.
# The rank prefixes of the reference files (d__, p__, ... in GG2; k__, p__, ... in UNITE) are removed.
# All ASVs are in the orientation of the reference sequences (forward primer first), so tryRC is not needed.
# Stops when the reference has more different taxonomies than assignTaxonomy can handle (see FUNS$unite_to_genus).
FUNS$assign_taxonomy <- function(sequences, reference, ranks, min_boot = 50, multithread = 20) {
  n_taxonomies <- FUNS$count_reference_taxonomies(reference)
  cat(reference, ":", n_taxonomies, "different taxonomies\n")
  if (n_taxonomies > 65536) {
    stop(reference, " has ", n_taxonomies, " different taxonomies; dada2::assignTaxonomy 1.40.0 is only correct up to 65,536.")
  }
  tax <- dada2::assignTaxonomy(sequences, reference, taxLevels = ranks, minBoot = min_boot,
                               outputBootstraps = TRUE, multithread = multithread, verbose = TRUE)
  taxon_names <- sub("^[a-z]__", "", tax$tax)
  taxon_names[taxon_names == ""] <- NA
  colnames(taxon_names) <- tolower(ranks)
  confidence <- tax$boot
  colnames(confidence) <- paste0(tolower(ranks), "_confidence")
  data.frame(taxon_names, confidence, row.names = NULL)
}


# Rule of thumb for a "real" taxon name (used for the assignment rates only). Not counted as real: names with
# Incertae sedis, unidentified, uncultured, metagenome or "_sp" (species unknown), and the codes given to groups without
# a Latin name, which contain digits (e.g. UBA2421 and Palsa-1150 in GG2, KD4-96 and Subgroup_2 in SILVA).
# The numbers that GG2 appends to some real names (e.g. Bradyrhizobium_503403) are ignored.
FUNS$is_real_name <- function(x) {
  !is.na(x) &
    !grepl("incertae|unidentified|uncultured|metagenome|_sp$", x, ignore.case = TRUE) &
    !grepl("[0-9]", sub("_[0-9]{4,}$", "", x))
}


# Share of the ASVs and of the reads that get a name at each rank, for a series of confidence cutoffs.
# tax: data frame from FUNS$assign_taxonomy; reads: total number of reads of each ASV (same order).
# "named": any name; "real_name": only names that pass FUNS$is_real_name.
FUNS$taxonomy_rates <- function(tax, reads, ranks, cutoffs = c(50, 60, 70, 80, 90)) {
  result <- NULL
  for (r in ranks) {
    for (cutoff in cutoffs) {
      named <- !is.na(tax[[r]]) & tax[[paste0(r, "_confidence")]] >= cutoff
      real  <- named & FUNS$is_real_name(tax[[r]])
      result <- rbind(result, data.frame(
        rank                = r,
        cutoff              = cutoff,
        asvs_named_pct      = round(100 * mean(named), 1),
        reads_named_pct     = round(100 * sum(reads[named]) / sum(reads), 1),
        asvs_real_name_pct  = round(100 * mean(real), 1),
        reads_real_name_pct = round(100 * sum(reads[real]) / sum(reads), 1)))
    }
  }
  result
}


# Plots the rates of FUNS$taxonomy_rates against the confidence cutoff, one line per rank
# (solid: any name; dashed: real names only); left the share of ASVs, right the share of reads.
FUNS$plot_taxonomy_rates <- function(rates, main) {
  ranks <- unique(rates$rank)
  par(mfrow = c(1, 2))
  for (what in c("asvs", "reads")) {
    plot(NA, xlim = range(rates$cutoff), ylim = c(0, 100), xlab = "confidence cutoff",
         ylab = paste("% of", what, "with a name"), main = main)
    for (i in seq_along(ranks)) {
      x <- rates[rates$rank == ranks[i], ]
      lines(x$cutoff, x[[paste0(what, "_named_pct")]],     col = i, lwd = 2)
      lines(x$cutoff, x[[paste0(what, "_real_name_pct")]], col = i, lwd = 2, lty = 2)
    }
    legend("bottomleft", legend = c(ranks, "any name", "real names only"), col = c(seq_along(ranks), 1, 1),
           lty = c(rep(1, length(ranks)), 1, 2), lwd = 2, bty = "n", cex = 0.7)
  }
}


# Agreement between two taxonomies of the same ASVs (tables from FUNS$assign_taxonomy, rows in the same order), per
# rank and confidence cutoff: among the ASVs to which both give a real name (FUNS$is_real_name) with a confidence of
# at least the cutoff, their number (asvs_compared) and the % with the same name (NA when there is nothing to compare).
FUNS$name_agreement <- function(a, b, ranks, cutoffs = c(50, 80)) {
  result <- NULL
  for (r in ranks) {
    for (cutoff in cutoffs) {
      both <- FUNS$is_real_name(a[[r]]) & FUNS$is_real_name(b[[r]]) &
        a[[paste0(r, "_confidence")]] >= cutoff & b[[paste0(r, "_confidence")]] >= cutoff
      result <- rbind(result, data.frame(
        rank          = r,
        cutoff        = cutoff,
        asvs_compared = sum(both),
        same_name_pct = if (any(both)) round(100 * mean(a[[r]][both] == b[[r]][both]), 1) else NA))
    }
  }
  result
}


# Writes a copy of a UNITE reference file with the taxonomy cut after the genus (the species field is removed).
# Why: dada2::assignTaxonomy (1.40.0) gives wrong results when a reference has more than 65,536 different taxonomies.
# In get_best_genus() in src/taxonomy.cpp, the position of taxonomy g in the table of k-mer probabilities is
# computed as g * n_kmers in 32-bit unsigned integers, with n_kmers = 4^8 = 65,536. From the 65,537th taxonomy on
# (g = 65,536, counting from 0) this exceeds 2^32 and wraps around, so the classifier reads the probabilities of
# another taxonomy. (The table itself is filled correctly: that code uses 64-bit integers.) Cut after the genus,
# the UNITE files have far fewer different taxonomies.
FUNS$unite_to_genus <- function(input, output) {
  ref <- Biostrings::readDNAStringSet(input)
  names(ref) <- sub(";s__[^;]*$", "", names(ref))
  Biostrings::writeXStringSet(ref, output)
}


# Number of different taxonomies in a reference file, counted as dada2::assignTaxonomy (1.40.0) counts them:
# sequences shorter than 20 nt are left out, UNITE headers are reduced to their taxonomy with the same substitutions
# as assignTaxonomy makes, and shorter taxonomies are padded to the full depth.
FUNS$count_reference_taxonomies <- function(reference, min_ref_len = 20) {
  ref <- Biostrings::readDNAStringSet(reference)
  tax <- trimws(names(ref)[Biostrings::width(ref) >= min_ref_len])
  if (all(grepl("FU\\|re[pf]s", tax[1:10]))) {   # UNITE format, as detected by assignTaxonomy
    tax <- sapply(strsplit(tax, "\\|"), `[`, 5)
    tax <- gsub("[pcofg]__unidentified;", "_DADA2_UNSPECIFIED;", tax)
    tax <- gsub(";s__(\\w+)_", ";s__", tax)
    tax <- gsub(";s__sp$", ";_DADA2_UNSPECIFIED", tax)
  }
  depth <- lengths(strsplit(tax, ";"))
  tax <- paste0(tax, strrep("_DADA2_UNSPECIFIED;", max(depth) - depth))
  length(unique(tax))
}


# Clusters the ASVs of a TSE into OTUs with vsearch (--cluster_size), at one or more identity thresholds.
# The ASVs are taken from the most to the least abundant (reads over all libraries of the TSE).
# An ASV that is at least `id` identical to an existing centroid joins it. Otherwise it becomes a new centroid.
# So the centroid of an OTU is its most abundant ASV.
# vsearch settings are its defaults. It compares the centroids in order of shared k-mers and stops at the first one
# that is close enough, which is usually the most similar one. (Comparing all centroids placed 18, 5 and 2 of 4,874
# ITS ASVs differently at 97, 98 and 99 %, and is far too slow for 16S.)
# Identity (--iddef 2): matching columns / alignment length, terminal gaps not counted.
# Returns a data frame with one row per ASV (same order as the TSE) and one column per threshold (otu97 for 0.97, ...).
# OTUs are named after their centroid, with "ASV" replaced by "OTU" plus the threshold:
# 16S_OTU97_00012 is the 97 % OTU around 16S_ASV00012.
FUNS$cluster_otus <- function(tse, ids = c(0.97, 0.98, 0.99), threads = 20) {
  if (!all(grepl("ASV", rownames(tse)))) stop("ASV names must contain 'ASV'.")
  fasta <- tempfile(fileext = ".fasta")
  uc    <- tempfile(fileext = ".uc")
  on.exit(unlink(c(fasta, uc)), add = TRUE)
  seqs <- TreeSummarizedExperiment::referenceSeq(tse)[rownames(tse)]
  names(seqs) <- paste0(rownames(tse), ";size=", rowSums(SummarizedExperiment::assay(tse, "counts")))
  Biostrings::writeXStringSet(seqs, fasta)

  result <- data.frame(row.names = rownames(tse))
  for (id in ids) {
    FUNS$run_tool("vsearch", c("--cluster_size", fasta, "--id", id, "--sizein", "--strand", "plus",
                               "--threads", threads, "--uc", uc, "--quiet"),
                  outputs = uc)
    # uc file: one line per ASV, "S" for a centroid, "H" for an ASV that joined the centroid in column 10
    x <- read.table(uc, sep = "\t", quote = "", comment.char = "", colClasses = "character")
    x <- x[x$V1 %in% c("S", "H"), ]
    asv      <- sub(";size=.*$", "", x$V9)
    centroid <- ifelse(x$V1 == "S", asv, sub(";size=.*$", "", x$V10))
    level    <- round(100 * id)
    result[[paste0("otu", level)]] <- sub("ASV", paste0("OTU", level, "_"), centroid[match(rownames(tse), asv)])
  }
  if (anyNA(as.matrix(result))) stop("vsearch did not place every ASV in an OTU.")
  result
}


# Filters a TSE for a downstream analysis. Returns a new TSE. The input object is not changed.
# Step 1 removes ASVs:
#   - ASVs with fewer than min_reads reads in total, summed over all samples
#   - ASVs found (at least 1 read) in fewer than min_samples samples
#   Both are counted over the samples in the TSE as it is given. So remove samples first, then filter.
#   The defaults (2 and 1) remove nothing from the objects of 03a and 03b, which only contain ASVs with 2 reads or more.
# Step 2 removes no ASVs. It replaces taxon names by `label`:
#   - when the confidence of that rank is below min_confidence
#   - when the name or its confidence is missing
#   - when the name itself means "no name". Only SILVA has such names in our objects: "Unknown Family", and
#     "Incertae Sedis" without the name of a higher rank (the same label in unrelated lineages).
#     GG2 and UNITE leave such ranks empty (NA), so they are covered by the rule above.
#   So after this step, "no name" is always `label`, never NA.
#   The confidence is the bootstrap value of dada2::assignTaxonomy (0-100): in how many of 100 tries, each with a
#   random eighth of the ASV's k-mers, the classifier chose the same taxon at that rank.
#   It never rises towards the lower ranks. So when a name is replaced for a low confidence,
#   the names below it are replaced too.
#   Each name is checked against the confidence of its own database:
#     family         with gg2_family_confidence (16S) or unite_family_confidence (ITS)
#     silva_family   with silva_family_confidence
#     unite_s_family with unite_s_family_confidence
#   The confidence columns themselves are kept as they are.
# The read totals in the colData are updated (FUNS$update_reads). The settings are stored in metadata(tse)$filter.
# Example: tse80 <- FUNS$filter_tse(readRDS("results/objects/TSE_16S_taxonomy.rds"), min_samples = 3)
FUNS$filter_tse <- function(tse, min_reads = 2, min_samples = 1, min_confidence = 80, label = "unclassified") {
  counts <- SummarizedExperiment::assay(tse, "counts")
  keep   <- rowSums(counts) >= min_reads & rowSums(counts > 0) >= min_samples
  message("filter_tse: kept ", sum(keep), " of ", length(keep), " ASVs, ",
          round(100 * sum(counts[keep, ]) / sum(counts), 2), " % of the reads")
  tse <- tse[keep, ]

  rd <- SummarizedExperiment::rowData(tse)
  for (conf_col in grep("_confidence$", names(rd), value = TRUE)) {
    # Name column of this confidence column. It is silva_family for silva_family_confidence,
    # but family (without prefix) for gg2_family_confidence and unite_family_confidence.
    name_col <- sub("_confidence$", "", conf_col)
    if (!name_col %in% names(rd)) name_col <- sub("^[^_]+_", "", name_col)
    if (!name_col %in% names(rd)) stop("No name column found for ", conf_col)
    low <- is.na(rd[[name_col]]) | is.na(rd[[conf_col]]) | rd[[conf_col]] < min_confidence |
      rd[[name_col]] %in% c("Unknown Family", "Incertae Sedis")
    rd[[name_col]][low] <- label
  }
  SummarizedExperiment::rowData(tse) <- rd

  tse <- FUNS$update_reads(tse)
  S4Vectors::metadata(tse)$filter <- list(min_reads = min_reads, min_samples = min_samples,
                                          min_confidence = min_confidence, label = label,
                                          asvs_before = length(keep), asvs_after = sum(keep))
  tse
}
