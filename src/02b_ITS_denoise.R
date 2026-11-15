# @ Rik Verdonck 20261027
# Purpose: denoising of the ITS reads with DADA2, from the raw reads to a table of ASVs.
#
# Same structure as 02a_16S_denoise.R. Differences with 16S:
#   * No truncation to a fixed length. The ITS2 insert varies in length (95 % of the inserts are 238-363 nt) and many
#     inserts are shorter than the reads, so their reads end where the insert ends. A fixed truncation length
#     would discard the reads of short inserts (a bias against taxa with a short ITS2), and the share of such reads
#     differs strongly between samples. 
#     Instead, the last 50 nt of the reverse reads are cut off (trimRight), where their quality drops.
#     The reads still overlap for inserts up to about 500 nt. Reverse reads may have up to 5 expected
#     errors (maxEE). These settings were compared with alternatives in side_analyses/02b_ITS_filter_settings.R.
#   * Default OMEGA_A and OMEGA_C. The strict OMEGA_A of 16S balances the pooling of the very diverse bacterial communities.
#     The alpha diversity of fungal communities is extreme, so we did not need this extra strict OMEGA_A.
#   * No length filter for the ASVs.
#
# As for 16S, the main analysis uses one DNA extract per soil sample, the original one (column primary in data/metadata.csv).
# The replicate extracts of 9 samples are left out, also from the ASV inference
# (they are used in side_analyses/02b_ITS_replicates.R). 
# The soil libraries are denoised together (pooled). 
# The two mock communities are denoised separately, with the error models of the soil libraries.
#
# Input:  data/metadata.csv, data/reads_ITS/input/ (from 01_input_reads.R)
# Output: results/objects/TSE_ITS.rds           TreeSummarizedExperiment of the 35 soil samples, passed on to the analyses
#         results/objects/TSE_ITS_mock.rds      the same for the 2 mocks
#         results/objects/asv_ITS.fasta         soil ASV sequences
#         results/objects/XITS.rds              all working objects of this script (for troubleshooting)
#         results/objects/derep_ITS.rds, derep_ITS_mock.rds   dereplicated reads (for troubleshooting)
#         results/tables/read_tracking_ITS.csv  read pairs left after each step, for every library
#         results/tables/error_model_convergence_ITS.csv   change of the error rates in each round of learnErrors
#         results/figures/ITS_*.pdf             quality profiles, error models, read tracking, ASV lengths
#         logs/cutadapt/
# In the R session, only XITS (plus FUNS and PATH) is left at the end.
# Runtime: about 20 minutes on 20 threads.

source("src/config.R")
source("src/functions/functions.R")
library("dada2")


### 1. Libraries and file names
###############################

XITS <- list()

# The 35 soil libraries (the primary extract of each sample) and the 2 mocks, mocks last
XITS$meta <- data.table::fread("data/metadata.csv", na.strings = "", data.table = FALSE)
XITS$meta <- XITS$meta[!XITS$meta$excluded & XITS$meta$primary & !is.na(XITS$meta$lib_ITS_forward), ]
XITS$meta <- XITS$meta[order(grepl("^mock", XITS$meta$soil_sample)), ]
stopifnot(sum(!grepl("^mock", XITS$meta$soil_sample)) == 35)

XITS$libraries     <- XITS$meta$DNA_extract
XITS$soil          <- !grepl("^mock", XITS$meta$soil_sample)   # TRUE for soil libraries, FALSE for the mocks
XITS$path.input    <- PATH$reads("ITS", "input")
XITS$path.cutadapt <- PATH$reads("ITS", "cutadapt")
XITS$path.filtered <- PATH$reads("ITS", "filtered")

XITS$fnFs     <- file.path(XITS$path.input,    XITS$meta$lib_ITS_forward)
XITS$fnRs     <- file.path(XITS$path.input,    XITS$meta$lib_ITS_reverse)
XITS$fnFs.cut <- file.path(XITS$path.cutadapt, XITS$meta$lib_ITS_forward)
XITS$fnRs.cut <- file.path(XITS$path.cutadapt, XITS$meta$lib_ITS_reverse)
XITS$filtFs   <- file.path(XITS$path.filtered, XITS$meta$lib_ITS_forward)
XITS$filtRs   <- file.path(XITS$path.filtered, XITS$meta$lib_ITS_reverse)

for (d in c(XITS$path.cutadapt, file.path(PATH$logs, "cutadapt"), PATH$objects, PATH$figures, PATH$tables)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
rm(d)

# Quality profiles of the raw reads (primers included)
FUNS$save_paired_quality_profiles(XITS$fnFs, XITS$fnRs, XITS$libraries, file.path(PATH$figures, "ITS_quality_profiles.pdf"))


### 2. Primer removal (cutadapt)
################################

# Primers ITS86F and ITS4 (5'-3'). The forward primer is written with two degenerate positions (R), so that it matches
# both ITS86F (GTGAATCATCGAATCTTTGAA) and its degenerate version.
# Many ITS2 inserts are shorter than the reads, so the reads often run through into the other primer.
# Cutadapt removes it at the 3' end (-a, -A), which leaves reads of variable length.
XITS$FWD    <- "GTGARTCATCGARTCTTTGAA"
XITS$REV    <- "TCCTCCGCTTATTGATATGC"
XITS$FWD.RC <- dada2:::rc(XITS$FWD)
XITS$REV.RC <- dada2:::rc(XITS$REV)

FUNS$run_tool("cutadapt", "--version", stdout = "", stderr = "")
for (i in seq_along(XITS$libraries)) {
  logfile <- file.path(PATH$logs, "cutadapt", paste0(XITS$libraries[i], "_ITS.log"))
  FUNS$run_tool("cutadapt",
                args = c("-g", XITS$FWD, "-a", XITS$REV.RC,
                         "-G", XITS$REV, "-A", XITS$FWD.RC,
                         "-n", 2, "--overlap", 10, "-e", 0.1, "--discard-untrimmed", "--pair-filter=any", "-j", 20,
                         "-o", XITS$fnFs.cut[i], "-p", XITS$fnRs.cut[i],
                         XITS$fnFs[i], XITS$fnRs[i]),
                outputs = c(XITS$fnFs.cut[i], XITS$fnRs.cut[i]), stdout = logfile, stderr = logfile)
}
rm(i, logfile)


### 3. Quality filtering
########################

# No truncation to a fixed length (see top). The last 50 nt of the reverse reads are cut off, where their quality drops.
# Reads are also cut at the first base with quality 2 or lower, and discarded when they are shorter than 50 nt,
# contain an N, or have more than 2 (forward) or 5 (reverse) expected errors.
XITS$out <- filterAndTrim(
  fwd = XITS$fnFs.cut,
  filt = XITS$filtFs,
  rev = XITS$fnRs.cut,
  filt.rev = XITS$filtRs,
  trimRight = c(0, 50),
  truncQ = 2,
  minLen = 50,
  maxEE = c(2, 5),
  maxN = 0,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = 20,
  matchIDs = TRUE
)


### 4. Error models
###################

# As for 16S: the first 12,000 filtered reads of every soil library (see functions.R). The mocks are not used here.
XITS$errF <- FUNS$learnErrorsBalanced(XITS$filtFs[XITS$soil], n_reads_per_sample = 12000, nbases = 1e8, multithread = 20, MAX_CONSIST = 15)
XITS$errR <- FUNS$learnErrorsBalanced(XITS$filtRs[XITS$soil], n_reads_per_sample = 12000, nbases = 1e8, multithread = 20, MAX_CONSIST = 15)

# How much the error rates still changed in each round of learnErrors (see FUNS$error_convergence)
data.table::fwrite(rbind(data.frame(reads = "forward", FUNS$error_convergence(XITS$errF)),
                         data.frame(reads = "reverse", FUNS$error_convergence(XITS$errR))),
                   file.path(PATH$tables, "error_model_convergence_ITS.csv"))

pdf(file.path(PATH$figures, "ITS_error_models.pdf"), width = 8, height = 8)
print(plotErrors(XITS$errF, nominalQ = TRUE) + ggplot2::ggtitle("ITS forward reads"))
print(plotErrors(XITS$errR, nominalQ = TRUE) + ggplot2::ggtitle("ITS reverse reads"))
dev.off()


### 5. Denoising, merging and chimera removal: soil libraries (pooled)
#######################################################################

# Default settings of dada() (OMEGA_A = 1e-40, OMEGA_C = 1e-40). Merging and chimera removal as for 16S.
XITS$derepFs <- derepFastq(XITS$filtFs[XITS$soil])
XITS$derepRs <- derepFastq(XITS$filtRs[XITS$soil])
names(XITS$derepFs) <- XITS$libraries[XITS$soil]
names(XITS$derepRs) <- XITS$libraries[XITS$soil]

XITS$dadaFs  <- dada(XITS$derepFs, err = XITS$errF, multithread = 20, pool = TRUE)
XITS$dadaRs  <- dada(XITS$derepRs, err = XITS$errR, multithread = 20, pool = TRUE)
XITS$mergers <- mergePairs(XITS$dadaFs, XITS$derepFs, XITS$dadaRs, XITS$derepRs, minOverlap = 12, maxMismatch = 1)
XITS$seqtab  <- makeSequenceTable(XITS$mergers)
XITS$nochim  <- removeBimeraDenovo(XITS$seqtab, method = "consensus", multithread = 20)

# The dereplicated reads take most of the memory and are not needed any more: we write them to disk and free the memory
saveRDS(list(fwd = XITS$derepFs, rev = XITS$derepRs), file.path(PATH$objects, "derep_ITS.rds"))
XITS$derepFs <- NULL
XITS$derepRs <- NULL
invisible(gc())


### 6. Denoising, merging and chimera removal: mocks (pooled, same settings)
############################################################################

XITS$derepFs.mock <- derepFastq(XITS$filtFs[!XITS$soil])
XITS$derepRs.mock <- derepFastq(XITS$filtRs[!XITS$soil])
names(XITS$derepFs.mock) <- XITS$libraries[!XITS$soil]
names(XITS$derepRs.mock) <- XITS$libraries[!XITS$soil]

XITS$dadaFs.mock  <- dada(XITS$derepFs.mock, err = XITS$errF, multithread = 20, pool = TRUE)
XITS$dadaRs.mock  <- dada(XITS$derepRs.mock, err = XITS$errR, multithread = 20, pool = TRUE)
XITS$mergers.mock <- mergePairs(XITS$dadaFs.mock, XITS$derepFs.mock, XITS$dadaRs.mock, XITS$derepRs.mock, minOverlap = 12, maxMismatch = 1)
XITS$seqtab.mock  <- makeSequenceTable(XITS$mergers.mock)
XITS$nochim.mock  <- removeBimeraDenovo(XITS$seqtab.mock, method = "consensus", multithread = 20)

saveRDS(list(fwd = XITS$derepFs.mock, rev = XITS$derepRs.mock), file.path(PATH$objects, "derep_ITS_mock.rds"))
XITS$derepFs.mock <- NULL
XITS$derepRs.mock <- NULL


### 7. Removing very rare ASVs
##############################

# We keep ASVs with at least 2 reads over all soil libraries. The filter is kept lenient on purpose:
# a stricter one can still be applied in the analyses, while removed reads cannot be brought back.
XITS$asv      <- XITS$nochim[, colSums(XITS$nochim) >= 2]
XITS$asv.mock <- XITS$nochim.mock[, colSums(XITS$nochim.mock) >= 2]

FUNS$plot_asv_lengths(XITS$asv, "ITS soil ASVs", file.path(PATH$figures, "ITS_asv_lengths.pdf"))


### 8. Read tracking
####################

XITS$track <- data.frame(
  DNA_extract   = XITS$libraries,
  soil_sample   = XITS$meta$soil_sample,
  primary       = XITS$meta$primary,
  input         = ShortRead::countFastq(XITS$fnFs)$records,
  primers_found = XITS$out[, 1],
  filtered      = XITS$out[, 2],
  denoised_fwd  = c(sapply(XITS$dadaFs, FUNS$getN),  sapply(XITS$dadaFs.mock, FUNS$getN)),
  denoised_rev  = c(sapply(XITS$dadaRs, FUNS$getN),  sapply(XITS$dadaRs.mock, FUNS$getN)),
  merged        = c(sapply(XITS$mergers, FUNS$getN), sapply(XITS$mergers.mock, FUNS$getN)),
  no_chimera    = c(rowSums(XITS$nochim),            rowSums(XITS$nochim.mock)),
  final         = c(rowSums(XITS$asv),               rowSums(XITS$asv.mock))
)
print(XITS$track)
data.table::fwrite(XITS$track, file.path(PATH$tables, "read_tracking_ITS.csv"))

FUNS$plot_read_tracking(XITS$track, steps = c("input", "primers_found", "filtered", "denoised_fwd", "denoised_rev", "merged", "no_chimera", "final"),
                        is_mock = !XITS$soil, main = "ITS: read pairs left after each step", output_file = file.path(PATH$figures, "ITS_read_tracking.pdf"))


### 9. Output: TreeSummarizedExperiment objects
###############################################

XITS$tse      <- FUNS$make_tse(XITS$asv,      XITS$meta[XITS$soil, ],  XITS$track[XITS$soil, ],  id_prefix = "ITS_ASV")
XITS$tse.mock <- FUNS$make_tse(XITS$asv.mock, XITS$meta[!XITS$soil, ], XITS$track[!XITS$soil, ], id_prefix = "ITS_mockASV")
XITS$tse
cat(nrow(XITS$tse), "soil ASVs and", nrow(XITS$tse.mock), "mock ASVs\n")

saveRDS(XITS$tse,      file.path(PATH$objects, "TSE_ITS.rds"))
saveRDS(XITS$tse.mock, file.path(PATH$objects, "TSE_ITS_mock.rds"))
Biostrings::writeXStringSet(TreeSummarizedExperiment::referenceSeq(XITS$tse), file.path(PATH$objects, "asv_ITS.fasta"))
saveRDS(XITS, file.path(PATH$objects, "XITS.rds"))
