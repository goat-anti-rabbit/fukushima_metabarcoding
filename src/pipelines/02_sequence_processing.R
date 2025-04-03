# @ Rik Verdonck 20250121
# Purpose: Preprocess raw sequencing data (DADA2 pipeline).

# Quality filtering of reads (trimming, filtering, etc.).
# Learn error rates and perform sequence inference (DADA2).
# Generate and save ASV tables for 16S and ITS.
# Outputs X16S and XITS, which are list objects that contain all filtering steps for both barcodes
# Outputs TSE_16S and TSE_ITS, which are TreeSummarizedExperiment objects for both barcodes

library("tidyverse")
library("data.table")
library("dada2")
library("gridExtra")  # For combining plots
library("doParallel")
library("TreeSummarizedExperiment")
library("S4Vectors")
library("SummarizedExperiment")


### 1. Reading in data, setting up XITS and X16S objects
########################################################

# Some small custom functions:
source("./src/functions/functions.R")

# Read in metadata
META <- fread("./data/intermediate/cleaned_metadata.csv")
rownames(META) <- META$DNA_extract

# Here I define 5 colors that we will use throughout for site, used for the levels 
# "Arakawa", "Futaba", "Okuma", "Tsushima" and "University" 
# I re-level META to put them in a more logical order as well. 
META$sampling_site   <- factor(META$sampling_site,levels=c("University", "Arakawa", "Tsushima", "Okuma", "Futaba"))
sitecolors  <- c("#ADD8E6","#2E8B57","#FF7F50","#FF91A4","#800020")

# We will for now save everything that has to do with ITS in the list object XITS, and everything 16S in the list object X16S
XITS <- list(
  path = "./data/raw/ITS_reads",
  path_intermediate = "./data/intermediate/ITS_reads"
)
X16S <- list(
  path = "./data/raw/16S_reads",
  path_intermediate = "./data/intermediate/16S_reads"
)

# Make sure the levels and order of META and samples in X objects are the same:
XITS$META <- META %>%
  filter(!is.na(lib_ITS_forward) & !is.na(lib_ITS_reverse)) %>%  # Keep rows where ITS libraries are available
  select(
    DNA_extract, soil_sample, sampling_data, sampling_site, soil_type,
    Cs_137, Cs_134, "lib_forward"=lib_ITS_forward, "lib_reverse"=lib_ITS_reverse,
    "primer_forward"=primer_ITS_forward, "primer_reverse"=primer_ITS_reverse, mock, "Nreads"=Nreads_ITS
  )
XITS$sample.names <- XITS$META$DNA_extract

X16S$META <- META %>%
  filter(!is.na(lib_16S_forward) & !is.na(lib_16S_reverse)) %>%  # Keep rows where 16S libraries are available
  select(
    DNA_extract, soil_sample, sampling_data, sampling_site, soil_type,
    Cs_137, Cs_134, "lib_forward"=lib_16S_forward, "lib_reverse"=lib_16S_reverse,
    "primer_forward"=primer_16S_forward, "primer_reverse"=primer_16S_reverse, mock, "Nreads"=Nreads_16S
  )
X16S$sample.names <- X16S$META$DNA_extract

# The entire list of libraries present in the paths:
XITS$fnFs <- sort(list.files(XITS$path, pattern = "_1.fq.gz", full.names = TRUE))
XITS$fnFs <- XITS$fnFs[match(XITS$sample.names, sapply(XITS$fnFs, FUNS$get.sample.name))]
#length(XITS$fnFs)
XITS$fnRs <- sort(list.files(XITS$path, pattern = "_2.fq.gz", full.names = TRUE))
XITS$fnRs <- XITS$fnRs[match(XITS$sample.names, sapply(XITS$fnRs, FUNS$get.sample.name))]
#length(XITS$fnRs)

X16S$fnFs <- sort(list.files(X16S$path, pattern = "_1.fq.gz", full.names = TRUE))
X16S$fnFs <- X16S$fnFs[match(X16S$sample.names, sapply(X16S$fnFs, FUNS$get.sample.name))]
#length(X16S$fnFs)
X16S$fnRs <- sort(list.files(X16S$path, pattern = "_2.fq.gz", full.names = TRUE))
X16S$fnRs <- X16S$fnRs[match(X16S$sample.names, sapply(X16S$fnRs, FUNS$get.sample.name))]
#length(X16S$fnRs)

# Remember samples A5T and T3T were discarded earlier (see script 01_data_preparation.R) 
# So here we just want to double check if everything is matching between data and metadata
sample.names         <- unique(c(XITS$sample.names,X16S$sample.names))
setdiff(sample.names,META$DNA_extract)
setdiff(META$DNA_extract,sample.names)

### 2. Quality profiles
#######################

# Output directory
output_dir <- "results/figures/qualityprofiles"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save paired quality profiles for ITS
FUNS$save_paired_quality_profiles(XITS$fnFs, XITS$fnRs, "ITS")
# Save paired quality profiles for 16S
FUNS$save_paired_quality_profiles(X16S$fnFs, X16S$fnRs, "16S")


### 3. Filtering: removing N and phix
#####################################

# Our first filtering step using the filterAndTrim function is really just getting rid of sequences containing N nucleotides, and PhiX. 
# About the multithread parameter: only use this if your computer allows parallel processing. 

XITS$fnFs.filtN <- paste(XITS$path_intermediate,"/",XITS$sample.names,"_ITS_1.filtN.fq.gz",sep="")
XITS$fnRs.filtN <- paste(XITS$path_intermediate,"/",XITS$sample.names,"_ITS_2.filtN.fq.gz",sep="")

X16S$fnFs.filtN <- paste(X16S$path_intermediate,"/",X16S$sample.names,"_16S_1.filtN.fq.gz",sep="")
X16S$fnRs.filtN <- paste(X16S$path_intermediate,"/",X16S$sample.names,"_16S_2.filtN.fq.gz",sep="")


XITS$out <- filterAndTrim(XITS$fnFs, XITS$fnFs.filtN, XITS$fnRs, XITS$fnRs.filtN, maxN = 0, matchIDs=T, rm.phix=T, multithread = 12)
X16S$out <- filterAndTrim(X16S$fnFs, X16S$fnFs.filtN, X16S$fnRs, X16S$fnRs.filtN, maxN = 0, matchIDs=T, rm.phix=T, multithread = 12)



### 4. Filtering: removing primers and investigating readthrough
################################################################

# Primer definitions
# 
# Primers are present and already confirmed by fastQC.
# All primers are stored in their standard 5′–3′ orientation, as written in lab protocols and publications.
# This is the orientation you would order from a supplier, NOT how they appear in sequencing reads.
# 
# IMPORTANT: Reverse primers will appear as reverse-complemented sequences in the actual reads.
# Therefore, you MUST reverse-complement primers before using them for:
#   - Matching/removal in sequencing reads (e.g., cutadapt trimming)
#   - Searching for primer hits in raw sequences
#
# To assist with this, we generate `.FWD.RC` and `.REV.RC` versions of all primers for internal use.
# 
# This design choice:
# * Keeps primer definitions readable and consistent with protocols
# * Matches conventions used in tools like QIIME2, MetaAmp, and DADA2 documentation
# * Allows easy comparison to published primer sets and metadata


#        "5′–GTGAATCATCGAATCTTTGAA–3′"
XITS$FWD <- "GTGARTCATCGARTCTTTGAA"
#         5′–TCCTCCGCTTATTGATATGC–3′
XITS$REV <- "TCCTCCGCTTATTGATATGC" 
#         5′–CCTACGGGNGGCWGCAG–3′
X16S$FWD <- "CCTACGGGNGGCWGCAG"
#         5′–GACTACHVGGGTATCTAATCC–3′
X16S$REV <- "GACTACHVGGGTATCTAATCC" 


# At first sight, there seems to be no need to turn them in all directions, but I agree it’s always better safe than sorry. 
# So we do check for primers in all orientations. 

XITS$FWD.orients <- FUNS$allOrients(XITS$FWD)
XITS$REV.orients <- FUNS$allOrients(XITS$REV)
X16S$FWD.orients <- FUNS$allOrients(X16S$FWD)
X16S$REV.orients <- FUNS$allOrients(X16S$REV)


# This function generates a table in which we count the number of reads that contain a primer in any direction. 
# The column names are a bit of a mess:
# Upper case FW or RV: whether we are looking for the forward or the reverse primers
# Lower case fw or rv: whether we are looking into the forward or reverse reads (you can for example find forward primers in reverse reads in case of readthrough)
# The direction of the primer. Whether it's oriented the way it should, reversed, complemented or reverse complemented. 

# For ITS:
XITS$primersummary <- FUNS$generate_primer_summary(XITS, "fnFs.filtN", "fnRs.filtN", multithread = 12)
# For 16S:
X16S$primersummary <- FUNS$generate_primer_summary(X16S, "fnFs.filtN", "fnRs.filtN", multithread = 12)


# Relative to number of reads:
XITS$primersummary_rel <- t(t(XITS$primersummary) / XITS$META$Nreads) 

# For ITS
XITS$proportion_revcomp_fw <- XITS$primersummary[,8]/XITS$primersummary[,1]
XITS$proportion_revcomp_rv <- XITS$primersummary[,12]/XITS$primersummary[,13]

# For 16S
X16S$proportion_revcomp_fw <- X16S$primersummary[,8]/X16S$primersummary[,1]
X16S$proportion_revcomp_rv <- X16S$primersummary[,12]/X16S$primersummary[,13]

# Diagnostic plots for % of readthrough:

pdf("results/figures/02_ITS_readthrough_diagnostics.pdf", width = 7, height = 7)
plot(XITS$proportion_revcomp_fw, XITS$proportion_revcomp_rv, log = "xy", col = "white", cex = 0.75,
     main = "ITS: proportion of reads with full primer readthrough")
text(XITS$proportion_revcomp_fw, XITS$proportion_revcomp_rv, labels = XITS$sample.names, cex = 0.75)
dev.off()

pdf("results/figures/02_16S_readthrough_diagnostics.pdf", width = 7, height = 7)
plot(X16S$proportion_revcomp_fw + 1e-06, X16S$proportion_revcomp_rv + 1e-06, log = "xy", col = "white", cex = 0.75,
     main = "16S: proportion of reads with full primer readthrough")
text(X16S$proportion_revcomp_fw + 1e-06, X16S$proportion_revcomp_rv + 1e-06, labels = X16S$sample.names, cex = 0.75)
dev.off()

# There is virtually no readthrough for the 16S samples, and when there is some, 
# it’s correlated between forward and reverse. 
# For plotting, I had to add a small quantity to the 16S counts in order to correct for zeros.

# For ITS, in the forward libraries, the proportion ranges from 0.5% to 25%. 
# In the reverse libraries, it’s between 0.9% and 42%. 
# The plot shows that forward and reverse are correlated, as expected.
# Samples C5T and C6T are clearly outliers in this respect, with very few detected readthroughs. 
# What we can tentatively conclude for ITS at this point:
#   * Some libraries are dominated by shorter inserts and others by longer ones.
#   * Given the substantial proportion of readthrough in most libraries, 
#     we can conclude that most inserts are shorter than 300 bp, 
#     which is lower than the expected fragment length of ~370 bp for basidiomycetes. 
#     This should be checked.
#   * Something is off in C5T and C6T. Perhaps base call quality is so bad that primers are not detected?
  

# Next step is cutadapt. I installed this on my local machine, and ran it via a system call:

cutadapt <- "/usr/bin/cutadapt" 
system2(cutadapt, args = "--version")

# Generate a subdirectory that will contain the cutadapt trimmed files
XITS$path.cut <- file.path(XITS$path, "cutadapt"); if(!dir.exists(XITS$path.cut)) dir.create(XITS$path.cut)
X16S$path.cut <- file.path(X16S$path, "cutadapt"); if(!dir.exists(X16S$path.cut)) dir.create(X16S$path.cut)

XITS$fnFs.cut <- file.path(XITS$path.cut, basename(XITS$fnFs))
XITS$fnRs.cut <- file.path(XITS$path.cut, basename(XITS$fnRs))

X16S$fnFs.cut <- file.path(X16S$path.cut, basename(X16S$fnFs))
X16S$fnRs.cut <- file.path(X16S$path.cut, basename(X16S$fnRs))

XITS$FWD.RC <- dada2:::rc(XITS$FWD)
XITS$REV.RC <- dada2:::rc(XITS$REV)

X16S$FWD.RC <- dada2:::rc(X16S$FWD)
X16S$REV.RC <- dada2:::rc(X16S$REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
XITS$R1.flags <- paste("-g", XITS$FWD, "-a", XITS$REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
XITS$R2.flags <- paste("-G", XITS$REV, "-A", XITS$FWD.RC)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
X16S$R1.flags <- paste("-g", X16S$FWD, "-a", X16S$REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
X16S$R2.flags <- paste("-G", X16S$REV, "-A", X16S$FWD.RC)


# Run Cutadapt for ITS
for(i in seq_along(XITS$fnFs)) {
  message("Running cutadapt on: ", XITS$fnFs[i])
  system2("cutadapt", args = c(XITS$R1.flags, XITS$R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", XITS$fnFs.cut[i], "-p", XITS$fnRs.cut[i], # output files
                               XITS$fnFs.filtN[i], XITS$fnRs.filtN[i])) # input files
}


# Run Cutadapt for 16S
for(i in seq_along(X16S$fnFs)) {
  message("Running cutadapt on: ", X16S$fnFs[i])
  system2("cutadapt", args = c(X16S$R1.flags, X16S$R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", X16S$fnFs.cut[i], "-p", X16S$fnRs.cut[i], # output files
                               X16S$fnFs.filtN[i], X16S$fnRs.filtN[i])) # input files
}


# For ITS:
XITS$primersummary.2 <- FUNS$generate_primer_summary(XITS, "fnFs.cut", "fnRs.cut")
# For 16S:
X16S$primersummary.2 <- FUNS$generate_primer_summary(X16S, "fnFs.cut", "fnRs.cut")


### 5. Filtering for base call quality
######################################

### Now we move on to read filtering. We split the analysis between 16S and ITS
### Because for 16S (virtually constant insert size) we can use a fixed read length cut-off
### While for ITS, inserts are of variable size and hence it makes less sense to cut them off at a fixed length
### (see here: https://benjjneb.github.io/dada2/ITS_workflow.html)

### 5A. Filtering the 16S reads for base call quality
#####################################################

X16S$filtFs <- file.path(X16S$path_intermediate, basename(X16S$fnFs.cut))
X16S$filtRs <- file.path(X16S$path_intermediate, basename(X16S$fnRs.cut))

X16S$out <- filterAndTrim(
  fwd = X16S$fnFs.cut,
  filt = X16S$filtFs,
  rev = X16S$fnRs.cut,
  filt.rev = X16S$filtRs,
  truncLen = c(250, 230),
  truncQ = 5,
  maxEE = c(2, 4),
  maxN = 0,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = 20,
  verbose = TRUE
)

# Learn error rates (using all filtered reads)
X16S$errF <- FUNS$learnErrorsBalanced(X16S$filtFs, n_reads_per_sample = 12000, nbases = 1.5e8, multithread = 20, MAX_CONSIST = 30)
X16S$errR <- FUNS$learnErrorsBalanced(X16S$filtRs, n_reads_per_sample = 12000, nbases = 1.5e8, multithread = 20, MAX_CONSIST = 30)

pdf("results/figures/02_errorprofile_estimates_16s.pdf", width = 10, height = 10)
par(mfrow = c(1, 2))
plotErrors(X16S$errF, nominalQ = TRUE)
plotErrors(X16S$errR, nominalQ = TRUE)
dev.off()

# Dereplicate
X16S$derepFs <- derepFastq(X16S$filtFs)
X16S$derepRs <- derepFastq(X16S$filtRs)

# Sample inference (DADA)
X16S$dadaFs <- dada(X16S$derepFs, err = X16S$errF, multithread = 20, pool = "pseudo")
X16S$dadaRs <- dada(X16S$derepRs, err = X16S$errR, multithread = 20, pool = "pseudo")

# Merging paired reads
X16S$mergers <- mergePairs(X16S$dadaFs, X16S$derepFs, X16S$dadaRs, X16S$derepRs, verbose = TRUE, maxMismatch=3)

# Construct ASV table
X16S$seqtab <- makeSequenceTable(X16S$mergers)
rownames(X16S$seqtab) <- X16S$sample.names

# Chimera removal
X16S$seqtab.nochim <- removeBimeraDenovo(X16S$seqtab, method = "consensus", multithread = 12, verbose = TRUE)
rownames(X16S$seqtab.nochim) <- X16S$sample.names

# Summary stats
cat("Reads after filtering:\n")
print(rowSums(X16S$out))

cat("Percentage of non-chimeric reads:\n")
print(sum(X16S$seqtab.nochim) / sum(X16S$seqtab))

pdf("results/figures/02_16S_ASV_abundance_vs_length.pdf", width = 7, height = 6)
FUNS$plot_asv_length_abundance(X16S$seqtab.nochim, X16S$META)
dev.off()

### Important: monitor drop-off at each step of the pipeline !!!
X16S$track <- cbind(X16S$out, sapply(X16S$dadaFs, FUNS$getN), sapply(X16S$dadaRs, FUNS$getN), sapply(X16S$mergers, FUNS$getN), rowSums(X16S$seqtab.nochim))
colnames(X16S$track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

FUNS$plot_read_retention(X16S$track, pdf_file = "results/figures/02_read_retention_16S_pipeline.pdf", title_prefix = "16S")


# Sequences between 380 and 390 seem to be archaons (BLASTn 100% identity to "uncultured archaeon 16S")
# The shorter asv's are a bit more vague

# save.image(file='session20250401.RData')


### 5B. Filtering the ITS reads for base call quality
#####################################################

XITS$filtFs <- file.path(XITS$path_intermediate, basename(XITS$fnFs.cut))
XITS$filtRs <- file.path(XITS$path_intermediate, basename(XITS$fnRs.cut))

XITS$out <- filterAndTrim(
  XITS$fnFs, XITS$filtFs, XITS$fnRs, XITS$filtRs,
  truncLen = c(0, 0),  # No truncation
  maxN = 0, maxEE = c(2, 4), truncQ = 8, minLen = 100, rm.phix = TRUE,
  compress = TRUE, multithread = 20
)

# Here we can already see that samples U3T_ITS_1.fq.gz (55899/67812 = 0.82) and C5T_ITS_1.fq.gz (66541/82142 = 0.81)     
# fall well outside the distribution of % reads retained that we see in the other libraries
# which is between 0.9 and 0.95. They are no outliers in absolute numbers of reads. 

# Learn error rates (using all filtered reads)
XITS$errF <- FUNS$learnErrorsBalanced(XITS$filtFs, n_reads_per_sample = 12000, nbases = 1.5e8, multithread = 20, MAX_CONSIST = 30)
XITS$errR <- FUNS$learnErrorsBalanced(XITS$filtRs, n_reads_per_sample = 12000, nbases = 1.5e8, multithread = 20, MAX_CONSIST = 30)

pdf("results/figures/02_errorprofile_estimates_its.pdf", width = 10, height = 10)
par(mfrow = c(1, 2))
plotErrors(XITS$errF, nominalQ = TRUE)
plotErrors(XITS$errR, nominalQ = TRUE)
dev.off()


# Dereplicate
XITS$derepFs <- derepFastq(XITS$filtFs)
XITS$derepRs <- derepFastq(XITS$filtRs)

# Sample inference (DADA)
XITS$dadaFs <- dada(XITS$derepFs, err = XITS$errF, multithread = 20, pool = "pseudo")
XITS$dadaRs <- dada(XITS$derepRs, err = XITS$errR, multithread = 20, pool = "pseudo")

# Merging paired reads
XITS$mergers <- mergePairs(XITS$dadaFs, XITS$derepFs, XITS$dadaRs, XITS$derepRs, verbose = TRUE, maxMismatch=3)

# Construct ASV table
XITS$seqtab <- makeSequenceTable(XITS$mergers)
rownames(XITS$seqtab) <- XITS$sample.names

# Chimera removal
XITS$seqtab.nochim <- removeBimeraDenovo(XITS$seqtab, method = "consensus", multithread = 12, verbose = TRUE)
rownames(XITS$seqtab.nochim) <- XITS$sample.names

### Monitor drop-off at each step of the pipeline !!!
XITS$track <- cbind(XITS$out, sapply(XITS$dadaFs, FUNS$getN), sapply(XITS$dadaRs, FUNS$getN), sapply(XITS$mergers, FUNS$getN), rowSums(XITS$seqtab.nochim))
colnames(XITS$track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

FUNS$plot_read_retention(XITS$track, pdf_file = "results/figures/02_read_retention_ITS_pipeline.pdf", title_prefix = "ITS")


pdf("results/figures/02_ITS_ASV_abundance_vs_length.pdf", width = 7, height = 6)
FUNS$plot_asv_length_abundance(XITS$seqtab.nochim, XITS$META)
dev.off()



### 6. Collapsing technical replicates
######################################

### Finally, we need to merge multiple extracts from the same sample. 
# See version 1 of our pipeline, where I demonstrate that there are no big differences between technical duplicates. 

# 16S: aggregate counts
x16s_counts <- FUNS$aggregate_counts_by_sample(X16S$seqtab.nochim, X16S$META)
x16s_meta   <- unique(X16S$META[, .(soil_sample, sampling_data, sampling_site, soil_type, Cs_137, Cs_134, mock)])
stopifnot(nrow(x16s_meta) == nrow(x16s_counts))
setkey(x16s_meta, soil_sample)
x16s_meta <- x16s_meta[rownames(x16s_counts)]

# Same for ITS:
ITS_counts <- FUNS$aggregate_counts_by_sample(XITS$seqtab.nochim, XITS$META)
ITS_meta   <- unique(XITS$META[, .(soil_sample, sampling_data, sampling_site, soil_type, Cs_137, Cs_134, mock)])
stopifnot(nrow(ITS_meta) == nrow(ITS_counts))
setkey(ITS_meta, soil_sample)
ITS_meta <- ITS_meta[rownames(ITS_counts)]


### 7. Create TreeSummarizedExperiment (TreeSE) objects
#######################################################

# First separate mock samples from the rest based on the metadata
# Filter out ASVs that are exclusively found in mocks
# Construct TreeSummarizedExperiment (TSE) objects for both ITS and 16S datasets
# Assign simple ASV identifiers (ASV1, ASV2, ...) while storing full sequences in rowData

# Resulting TSEs:
# - Rows: ASVs (identified as ASV1, ASV2, ...)  
# - Columns: Samples (soil samples or mock samples)  
# - assay(): matrix of raw counts (ASVs x samples)  
# - colData(): per-sample metadata (e.g. site, radiation level, mock status)  
# - rowData(): contains the full ASV DNA sequence corresponding to each ASV ID  


# ---- ITS ----
ITS_split      <- FUNS$split_TSE_remove_mockASVs(XITS$seqtab.nochim, XITS$META)
ITS_counts     <- ITS_split$main$counts
ITS_meta       <- ITS_split$main$meta
mockITS_counts <- ITS_split$mock$counts
mockITS_meta   <- ITS_split$mock$meta

TSE_ITS        <- FUNS$make_TSE(ITS_counts, ITS_meta, prefix = "its")
TSE_ITS_mock   <- FUNS$make_TSE(mockITS_counts, mockITS_meta, prefix = "mock_its")


# ---- 16S ----
X16S_split     <- FUNS$split_TSE_remove_mockASVs(X16S$seqtab.nochim, X16S$META)
x16s_counts    <- X16S_split$main$counts
x16s_meta      <- X16S_split$main$meta
mock16S_counts <- X16S_split$mock$counts
mock16S_meta   <- X16S_split$mock$meta

TSE_16S        <- FUNS$make_TSE(x16s_counts, x16s_meta, prefix = "16s")
TSE_16S_mock   <- FUNS$make_TSE(mock16S_counts, mock16S_meta, prefix = "mock_16s")









# Create TSE for 16S
X16S_TSE <- TreeSummarizedExperiment::TreeSummarizedExperiment(
  assays = list(counts = t(x16s_counts)),
  colData = DataFrame(x16s_meta),
  rowData = DataFrame(sequence = colnames(x16s_counts))
)

rownames(rowData(X16S_TSE)) <- paste0("ASV", seq_len(nrow(X16S_TSE)))
# Get current sequences
asv_seqs <- rownames(X16S_TSE)
# Create new simple ASV names
asv_ids <- paste0("ASV", seq_along(asv_seqs))
# Replace rownames in assay and rowData
rownames(X16S_TSE) <- asv_ids
# Store actual sequences as a rowData field
rowData(X16S_TSE)$sequence <- asv_seqs


# Create TSE for ITS
XITS_TSE <- TreeSummarizedExperiment::TreeSummarizedExperiment(
  assays = list(counts = t(ITS_counts)),
  colData = DataFrame(ITS_meta),
  rowData = DataFrame(sequence = colnames(ITS_counts))
)

rownames(rowData(XITS_TSE)) <- paste0("ASV", seq_len(nrow(XITS_TSE)))
asv_seqs <- rownames(XITS_TSE)
asv_ids <- paste0("ASV", seq_along(asv_seqs))
rownames(XITS_TSE) <- asv_ids
rowData(XITS_TSE)$sequence <- asv_seqs

pdf("results/figures/02_ASV_counts.pdf", width = 7, height = 5)
par(mfrow=c(1,2))
plot(rowSums(assay(XITS_TSE)),log="y", main = "Sum of ASV counts for ITS", ylab = "counts per ASV",type="l", lwd = 2, ylim = c(1,5e05))
plot(rowSums(assay(X16S_TSE)),log="y", main = "Sum of ASV counts for 16S", ylab = "counts per ASV",type="l", lwd = 2, ylim = c(1,5e05))
dev.off()



### 8. Save all files: 
######################
# Save the denoising environment (XITS and X16S contain all intermediate steps)
save(XITS, X16S, file = "results/objects/02_denoising_workspace.RData")

# Save the final TreeSE objects for downstream analysis
save(TSE_ITS, TSE_ITS_mock, TSE_16S, TSE_16S_mock, file = "results/objects/02_TreeSE_objects.RData")


