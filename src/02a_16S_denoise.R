# @ Rik Verdonck 20261027
# Purpose: denoising of the 16S reads with DADA2, from the raw reads to a table of ASVs.
#
# The main analysis uses one DNA extract per soil sample, the original one (column primary in data/metadata.csv),
# so that all samples have the same sampling effort. 
# The replicate extracts of 10 samples are left out here, also from the ASV inference: 
# they would add reads (and so ASV resolution) to only some samples and sites.
# They are used in side_analyses/02a_16S_replicates.R.
# The soil libraries are denoised together (pooled), so that a variant that occurs in several samples is recognised
# even when it is rare in each of them. The two mock communities are denoised separately, with the error models
# of the soil libraries, so that they cannot influence the soil ASVs.
#
# Input:  data/metadata.csv, data/reads_16S/input/ (from 01_input_reads.R)
# Output: results/objects/TSE_16S.rds           TreeSummarizedExperiment of the 39 soil samples, passed on to the analyses
#         results/objects/TSE_16S_mock.rds      the same for the 2 mocks
#         results/objects/asv_16S.fasta         soil ASV sequences
#         results/objects/X16S.rds              all working objects of this script (for troubleshooting)
#         results/objects/derep_16S.rds, derep_16S_mock.rds   dereplicated reads (for troubleshooting)
#         results/tables/read_tracking_16S.csv  read pairs left after each step, for every library
#         results/tables/error_model_convergence_16S.csv   change of the error rates in each round of learnErrors
#         results/figures/16S_*.pdf             quality profiles, error models, read tracking, ASV lengths
#         logs/cutadapt/
# In the R session, only X16S (plus FUNS and PATH) is left at the end.
# Runtime: about 5 hours on 20 threads (error models about 2.5 h, denoising about 2 h), about 25 GB of memory.

source("src/config.R")
source("src/functions/functions.R")
library("dada2")


### 1. Libraries and file names
###############################

X16S <- list()

# The 39 soil libraries (the primary extract of each sample) and the 2 mocks, mocks last
X16S$meta <- data.table::fread("data/metadata.csv", na.strings = "", data.table = FALSE)
X16S$meta <- X16S$meta[!X16S$meta$excluded & X16S$meta$primary & !is.na(X16S$meta$lib_16S_forward), ]
X16S$meta <- X16S$meta[order(grepl("^mock", X16S$meta$soil_sample)), ]
stopifnot(sum(!grepl("^mock", X16S$meta$soil_sample)) == 39)

X16S$libraries     <- X16S$meta$DNA_extract
X16S$soil          <- !grepl("^mock", X16S$meta$soil_sample)   # TRUE for soil libraries, FALSE for the mocks
X16S$path.input    <- PATH$reads("16S", "input")
X16S$path.cutadapt <- PATH$reads("16S", "cutadapt")
X16S$path.filtered <- PATH$reads("16S", "filtered")

X16S$fnFs     <- file.path(X16S$path.input,    X16S$meta$lib_16S_forward)
X16S$fnRs     <- file.path(X16S$path.input,    X16S$meta$lib_16S_reverse)
X16S$fnFs.cut <- file.path(X16S$path.cutadapt, X16S$meta$lib_16S_forward)
X16S$fnRs.cut <- file.path(X16S$path.cutadapt, X16S$meta$lib_16S_reverse)
X16S$filtFs   <- file.path(X16S$path.filtered, X16S$meta$lib_16S_forward)
X16S$filtRs   <- file.path(X16S$path.filtered, X16S$meta$lib_16S_reverse)

for (d in c(X16S$path.cutadapt, file.path(PATH$logs, "cutadapt"), PATH$objects, PATH$figures, PATH$tables)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
rm(d)

# Quality profiles of the raw reads (primers included: 17 nt forward, 21 nt reverse).
# These are the basis for the truncation lengths in section 3.
FUNS$save_paired_quality_profiles(X16S$fnFs, X16S$fnRs, X16S$libraries, file.path(PATH$figures, "16S_quality_profiles.pdf"))


### 2. Primer removal (cutadapt)
################################

# Primers 341F and 805R, as written in the lab protocol (5'-3').
# In the reads, the reverse primer appears as its reverse complement, hence the .RC versions.
X16S$FWD    <- "CCTACGGGNGGCWGCAG"
X16S$REV    <- "GACTACHVGGGTATCTAATCC"
X16S$FWD.RC <- dada2:::rc(X16S$FWD)
X16S$REV.RC <- dada2:::rc(X16S$REV)

# Read 1 starts with the forward primer, and ends with the reverse complement of the reverse primer
# if the insert is short enough to read through. Read 2 the other way round.
# Pairs in which a primer is missing from either read are discarded.
FUNS$run_tool("cutadapt", "--version", stdout = "", stderr = "")
for (i in seq_along(X16S$libraries)) {
  logfile <- file.path(PATH$logs, "cutadapt", paste0(X16S$libraries[i], "_16S.log"))
  FUNS$run_tool("cutadapt",
                args = c("-g", X16S$FWD, "-a", X16S$REV.RC,
                         "-G", X16S$REV, "-A", X16S$FWD.RC,
                         "-n", 2, "--overlap", 10, "-e", 0.1, "--discard-untrimmed", "--pair-filter=any", "-j", 20,
                         "-o", X16S$fnFs.cut[i], "-p", X16S$fnRs.cut[i],
                         X16S$fnFs[i], X16S$fnRs[i]),
                outputs = c(X16S$fnFs.cut[i], X16S$fnRs.cut[i]), stdout = logfile, stderr = logfile)
}
rm(i, logfile)


### 3. Quality filtering
########################

# Reads are cut to 250 (forward) and 200 (reverse) nucleotides. The insert is 401-430 nt long,
# so the two reads still overlap by at least 20 nt. A read is discarded when it becomes shorter than that
# (cut at the first base with quality 5 or lower), when it contains an N, or when it has more than
# 1 expected error. Reads from the PhiX control are removed too.
X16S$out <- filterAndTrim(
  fwd = X16S$fnFs.cut,
  filt = X16S$filtFs,
  rev = X16S$fnRs.cut,
  filt.rev = X16S$filtRs,
  truncLen = c(250, 200),
  truncQ = 5,
  maxEE = c(1, 1),
  maxN = 0,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = 20,
  matchIDs = TRUE
)


### 4. Error models
###################

# learnErrors() would only use the first few libraries. Instead, we take the first 12,000 filtered reads
# of every soil library, so that all samples contribute equally (see functions.R). The mocks are not used here.
X16S$errF <- FUNS$learnErrorsBalanced(X16S$filtFs[X16S$soil], n_reads_per_sample = 12000, nbases = 1e8, multithread = 20, MAX_CONSIST = 15)
X16S$errR <- FUNS$learnErrorsBalanced(X16S$filtRs[X16S$soil], n_reads_per_sample = 12000, nbases = 1e8, multithread = 20, MAX_CONSIST = 15)

# How much the error rates still changed in each round of learnErrors (see FUNS$error_convergence)
data.table::fwrite(rbind(data.frame(reads = "forward", FUNS$error_convergence(X16S$errF)),
                         data.frame(reads = "reverse", FUNS$error_convergence(X16S$errR))),
                   file.path(PATH$tables, "error_model_convergence_16S.csv"))

pdf(file.path(PATH$figures, "16S_error_models.pdf"), width = 8, height = 8)
print(plotErrors(X16S$errF, nominalQ = TRUE) + ggplot2::ggtitle("16S forward reads"))
print(plotErrors(X16S$errR, nominalQ = TRUE) + ggplot2::ggtitle("16S reverse reads"))
dev.off()


### 5. Denoising, merging and chimera removal: soil libraries (pooled)
#######################################################################

# OMEGA_A = 1e-100: a sequence only becomes a new ASV with very strong support (default 1e-40).
#   This balances the pooling, which gives rare variants much more statistical support than a single library would.
# OMEGA_C = 1e-40 (the default): a read is assigned to its ASV, unless it clearly does not fit (then it is dropped).
# MIN_ABUNDANCE = 2: an ASV must be seen at least twice.
# At merging, a pair is discarded when its two reads were assigned to ASVs that disagree where they overlap.
# These are mostly PCR chimeras, and reads of rare variants that were not accepted as ASVs
# (side_analyses/02a_16S_lost_reads.R shows that losing them does not change the results).
X16S$derepFs <- derepFastq(X16S$filtFs[X16S$soil])
X16S$derepRs <- derepFastq(X16S$filtRs[X16S$soil])
names(X16S$derepFs) <- X16S$libraries[X16S$soil]
names(X16S$derepRs) <- X16S$libraries[X16S$soil]

X16S$dadaFs  <- dada(X16S$derepFs, err = X16S$errF, multithread = 20, pool = TRUE, OMEGA_A = 1e-100, OMEGA_C = 1e-40, MIN_ABUNDANCE = 2)
X16S$dadaRs  <- dada(X16S$derepRs, err = X16S$errR, multithread = 20, pool = TRUE, OMEGA_A = 1e-100, OMEGA_C = 1e-40, MIN_ABUNDANCE = 2)
X16S$mergers <- mergePairs(X16S$dadaFs, X16S$derepFs, X16S$dadaRs, X16S$derepRs, minOverlap = 12, maxMismatch = 1)
X16S$seqtab  <- makeSequenceTable(X16S$mergers)
X16S$nochim  <- removeBimeraDenovo(X16S$seqtab, method = "consensus", multithread = 20)

# The dereplicated reads take most of the memory and are not needed any more: we write them to disk and free the memory
saveRDS(list(fwd = X16S$derepFs, rev = X16S$derepRs), file.path(PATH$objects, "derep_16S.rds"))
X16S$derepFs <- NULL
X16S$derepRs <- NULL
invisible(gc())


### 6. Denoising, merging and chimera removal: mocks (pooled, same settings)
############################################################################

X16S$derepFs.mock <- derepFastq(X16S$filtFs[!X16S$soil])
X16S$derepRs.mock <- derepFastq(X16S$filtRs[!X16S$soil])
names(X16S$derepFs.mock) <- X16S$libraries[!X16S$soil]
names(X16S$derepRs.mock) <- X16S$libraries[!X16S$soil]

X16S$dadaFs.mock  <- dada(X16S$derepFs.mock, err = X16S$errF, multithread = 20, pool = TRUE, OMEGA_A = 1e-100, OMEGA_C = 1e-40, MIN_ABUNDANCE = 2)
X16S$dadaRs.mock  <- dada(X16S$derepRs.mock, err = X16S$errR, multithread = 20, pool = TRUE, OMEGA_A = 1e-100, OMEGA_C = 1e-40, MIN_ABUNDANCE = 2)
X16S$mergers.mock <- mergePairs(X16S$dadaFs.mock, X16S$derepFs.mock, X16S$dadaRs.mock, X16S$derepRs.mock, minOverlap = 12, maxMismatch = 1)
X16S$seqtab.mock  <- makeSequenceTable(X16S$mergers.mock)
X16S$nochim.mock  <- removeBimeraDenovo(X16S$seqtab.mock, method = "consensus", multithread = 20)

saveRDS(list(fwd = X16S$derepFs.mock, rev = X16S$derepRs.mock), file.path(PATH$objects, "derep_16S_mock.rds"))
X16S$derepFs.mock <- NULL
X16S$derepRs.mock <- NULL


### 7. Removing very rare and short ASVs
########################################

# We keep ASVs with at least 2 reads over all soil libraries, and of at least 351 nt. The V3-V4 insert is
# 401-430 nt, so much shorter sequences are off-target products. The read filter is kept lenient on purpose:
# a stricter one can still be applied in the analyses, while removed reads cannot be brought back.
X16S$asv      <- X16S$nochim[, colSums(X16S$nochim) >= 2 & nchar(colnames(X16S$nochim)) >= 351]
X16S$asv.mock <- X16S$nochim.mock[, colSums(X16S$nochim.mock) >= 2 & nchar(colnames(X16S$nochim.mock)) >= 351]

FUNS$plot_asv_lengths(X16S$asv, "16S soil ASVs", file.path(PATH$figures, "16S_asv_lengths.pdf"))


### 8. Read tracking
####################

# Read pairs left after each step, for every library
X16S$track <- data.frame(
  DNA_extract   = X16S$libraries,
  soil_sample   = X16S$meta$soil_sample,
  primary       = X16S$meta$primary,
  input         = ShortRead::countFastq(X16S$fnFs)$records,
  primers_found = X16S$out[, 1],
  filtered      = X16S$out[, 2],
  denoised_fwd  = c(sapply(X16S$dadaFs, FUNS$getN),  sapply(X16S$dadaFs.mock, FUNS$getN)),
  denoised_rev  = c(sapply(X16S$dadaRs, FUNS$getN),  sapply(X16S$dadaRs.mock, FUNS$getN)),
  merged        = c(sapply(X16S$mergers, FUNS$getN), sapply(X16S$mergers.mock, FUNS$getN)),
  no_chimera    = c(rowSums(X16S$nochim),            rowSums(X16S$nochim.mock)),
  final         = c(rowSums(X16S$asv),               rowSums(X16S$asv.mock))
)
print(X16S$track)
data.table::fwrite(X16S$track, file.path(PATH$tables, "read_tracking_16S.csv"))

FUNS$plot_read_tracking(X16S$track, steps = c("input", "primers_found", "filtered", "denoised_fwd", "denoised_rev", "merged", "no_chimera", "final"),
                        is_mock = !X16S$soil, main = "16S: read pairs left after each step", output_file = file.path(PATH$figures, "16S_read_tracking.pdf"))


### 9. Output: TreeSummarizedExperiment objects
###############################################

# One object for the soil libraries and one for the mocks: ASVs as rows (16S_ASV00001 = most abundant),
# libraries as columns, the metadata and read tracking of each library in colData, the sequences in referenceSeq.
X16S$tse      <- FUNS$make_tse(X16S$asv,      X16S$meta[X16S$soil, ],  X16S$track[X16S$soil, ],  id_prefix = "16S_ASV")
X16S$tse.mock <- FUNS$make_tse(X16S$asv.mock, X16S$meta[!X16S$soil, ], X16S$track[!X16S$soil, ], id_prefix = "16S_mockASV")
X16S$tse
cat(nrow(X16S$tse), "soil ASVs and", nrow(X16S$tse.mock), "mock ASVs\n")

saveRDS(X16S$tse,      file.path(PATH$objects, "TSE_16S.rds"))
saveRDS(X16S$tse.mock, file.path(PATH$objects, "TSE_16S_mock.rds"))
Biostrings::writeXStringSet(TreeSummarizedExperiment::referenceSeq(X16S$tse), file.path(PATH$objects, "asv_16S.fasta"))
saveRDS(X16S, file.path(PATH$objects, "X16S.rds"))
