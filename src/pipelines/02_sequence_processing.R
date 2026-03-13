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
library("dplyr")
library("ShortRead")


### 0. Read in necessary functions. 
#############################################################
source("src/functions/functions.R")

### 1. Reading in data, setting up XITS and X16S objects
########################################################

# Read in metadata
meta_16S <- fread("./data/metadata_16S.csv")
meta_ITS <- fread("./data/metadata_ITS.csv")

rownames(meta_16S) <- meta_16S$DNA_extract
rownames(meta_ITS) <- meta_ITS$DNA_extract



# Here I define 5 colors that we will use throughout for site, used for the levels 
# "Arakawa", "Futaba", "Okuma", "Tsushima" and "University" 
# I re-level META to put them in a more logical order as well. 
meta_16S$sampling_site   <- factor(meta_16S$sampling_site,levels=c("University", "Arakawa", "Tsushima", "Okuma", "Futaba"))
meta_ITS$sampling_site   <- factor(meta_ITS$sampling_site,levels=c("University", "Arakawa", "Tsushima", "Okuma", "Futaba"))

sitecolors  <- c("#ADD8E6","#2E8B57","#FF7F50","#FF91A4","#800020")

# We will for now save everything that has to do with ITS in the list object XITS, and everything 16S in the list object X16S

X16S <- list(
  path.raw = "./data/reads_16S/raw",
  path.cutadapt = "./data/reads_16S/cutadapt",
  path.filtered = "./data/reads_16S/filtered",
  meta = meta_16S, 
  sample.names = meta_16S$DNA_extract
)

X16S$fnFs <- sort(list.files(X16S$path.raw, pattern = "_1.fq.gz", full.names = TRUE))
X16S$fnFs <- X16S$fnFs[match(X16S$sample.names, sapply(X16S$fnFs, FUNS$get.sample.name))]
#length(X16S$fnFs)
X16S$fnRs <- sort(list.files(X16S$path.raw, pattern = "_2.fq.gz", full.names = TRUE))
X16S$fnRs <- X16S$fnRs[match(X16S$sample.names, sapply(X16S$fnRs, FUNS$get.sample.name))]
#length(X16S$fnRs)

XITS <- list(
  path.raw = "./data/reads_ITS/raw",
  path.cutadapt = "./data/reads_ITS/cutadapt",
  path.filtered = "./data/reads_ITS/filtered",
  meta = meta_ITS,
  sample.names = meta_ITS$DNA_extract
)

# The entire list of libraries present in the paths:
XITS$fnFs <- sort(list.files(XITS$path.raw, pattern = "_1.fq.gz", full.names = TRUE))
XITS$fnFs <- XITS$fnFs[match(XITS$sample.names, sapply(XITS$fnFs, FUNS$get.sample.name))]
#length(XITS$fnFs)
XITS$fnRs <- sort(list.files(XITS$path.raw, pattern = "_2.fq.gz", full.names = TRUE))
XITS$fnRs <- XITS$fnRs[match(XITS$sample.names, sapply(XITS$fnRs, FUNS$get.sample.name))]
#length(XITS$fnRs)


# Remember samples A5T and T3T were discarded earlier (see script 01_data_preparation.R) 
# So here we just want to double check if everything is matching between data and metadata
sample.names         <- unique(c(XITS$sample.names,X16S$sample.names))
setdiff(sample.names, meta_16S$DNA_extract)
setdiff(meta_16S$DNA_extract,sample.names)
setdiff(sample.names, meta_16S$DNA_extract)
setdiff(meta_16S$DNA_extract,sample.names)


### 2. Quality profiles
#######################

# Output directory
output_dir <- "results/figures/qualityprofiles"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save paired quality profiles for 16S
FUNS$save_paired_quality_profiles(X16S$fnFs, X16S$fnRs, "16S")
# Save paired quality profiles for ITS
FUNS$save_paired_quality_profiles(XITS$fnFs, XITS$fnRs, "ITS")




### 3. Primers
##############
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

#                5′–CCTACGGGNGGCWGCAG–3′
X16S$FWD     <-    "CCTACGGGNGGCWGCAG"
X16S$FWD.RC  <- dada2:::rc(X16S$FWD)
#                5′–GACTACHVGGGTATCTAATCC–3′
X16S$REV     <-    "GACTACHVGGGTATCTAATCC" 
X16S$REV.RC  <- dada2:::rc(X16S$REV)

#               "5′–GTGAATCATCGAATCTTTGAA–3′"
XITS$FWD     <-    "GTGARTCATCGARTCTTTGAA"
XITS$FWD.RC  <- dada2:::rc(XITS$FWD)
#                5′–TCCTCCGCTTATTGATATGC–3′
XITS$REV     <-    "TCCTCCGCTTATTGATATGC" 
XITS$REV.RC  <- dada2:::rc(XITS$REV)

# At first sight, there seems to be no need to turn them in all directions, but I agree it’s always better safe than sorry. 
# So we do check for primers in all orientations. 
XITS$FWD.orients <- FUNS$allOrients(XITS$FWD)
XITS$REV.orients <- FUNS$allOrients(XITS$REV)
X16S$FWD.orients <- FUNS$allOrients(X16S$FWD)
X16S$REV.orients <- FUNS$allOrients(X16S$REV)


# This function generates a table in which we count the number of reads that contain a primer in any direction. 
# The column names:
# R1 or R2: whether we are looking into the first (forward) or second (reverse) reads 
# (e.g. you can find forward primers in reverse reads in case of readthrough)
# FP or RP: whether we are looking for the forward primer or the reverse primer
# The direction of the primer: 
# Whether it's oriented the way it should (5` to 3`), complement, reversed (3` to 5`), or reverse complemented. 

# For 16S:
X16S$primersummary <- FUNS$generate_primer_summary(X16S, "fnFs", "fnRs", multithread = 12)
# For ITS:
XITS$primersummary <- FUNS$generate_primer_summary(XITS, "fnFs", "fnRs", multithread = 12)

# Relative to number of reads:
X16S$primersummary_rel <- X16S$primersummary/ X16S$meta$Nreads 
XITS$primersummary_rel <- XITS$primersummary/ XITS$meta$Nreads 

# Readthrough:

# For 16S
X16S$readthrough_R1 <- X16S$primersummary[,"R1_RP.rvcmp"]/X16S$meta$Nreads 
X16S$readthrough_R2 <- X16S$primersummary[,"R2_FP.rvcmp"]/X16S$meta$Nreads
# For ITS
XITS$readthrough_R1 <- XITS$primersummary[,"R1_RP.rvcmp"]/XITS$meta$Nreads 
XITS$readthrough_R2 <- XITS$primersummary[,"R2_FP.rvcmp"]/XITS$meta$Nreads


# Diagnostic plots for % of read through:
pdf("results/figures/02_16S_readthrough_diagnostics.pdf", width = 7, height = 7)
plot(X16S$readthrough_R2 + 1e-05, X16S$readthrough_R1 + 1e-05, log = "xy", col = "white", cex = 0.75,
     main = "16S: proportion of reads with full primer readthrough",xlab="proportion readthrough read 2", ylab="proportion readthrough read 1")
text(X16S$readthrough_R2 + 1e-05, X16S$readthrough_R1 + 1e-05, labels = X16S$sample.names, cex = 0.75)
dev.off()

pdf("results/figures/02_ITS_readthrough_diagnostics.pdf", width = 7, height = 7)
plot(XITS$readthrough_R2, XITS$readthrough_R1, log = "xy", col = "white", cex = 0.75,
     main = "ITS: proportion of reads with full primer readthrough",xlab="proportion readthrough read 2", ylab="proportion readthrough read 1")
text(XITS$readthrough_R2, XITS$readthrough_R1, labels = XITS$sample.names, cex = 0.75)
dev.off()


# There is virtually no readthrough for the 16S samples, and when there is some, 
# it’s correlated between forward and reverse. 
# For plotting, I had to add a small quantity to the 16S counts in order to correct for zeros.

# For ITS, both in the forward and reverse libraries, the proportion ranges from 0.5% to 25%. 
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

cutadapt <- "/home/rik/miniconda3/bin/cutadapt" 
system2(cutadapt, args = "--version")

X16S$fnFs.cut <- file.path(X16S$path.cutadapt, basename(X16S$fnFs)); dir.create(X16S$path.cutadapt, showWarnings = FALSE)
X16S$fnRs.cut <- file.path(X16S$path.cutadapt, basename(X16S$fnRs)); dir.create(X16S$path.cutadapt, showWarnings = FALSE)

XITS$fnFs.cut <- file.path(XITS$path.cutadapt, basename(XITS$fnFs)); dir.create(XITS$path.cutadapt, showWarnings = FALSE)
XITS$fnRs.cut <- file.path(XITS$path.cutadapt, basename(XITS$fnRs)); dir.create(XITS$path.cutadapt, showWarnings = FALSE)

# Next, generate the flags that will be part of the parameters of cutadapt:
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
X16S$R1.flags <- paste("-g", X16S$FWD, "-a", X16S$REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
X16S$R2.flags <- paste("-G", X16S$REV, "-A", X16S$FWD.RC)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
XITS$R1.flags <- paste("-g", XITS$FWD, "-a", XITS$REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
XITS$R2.flags <- paste("-G", XITS$REV, "-A", XITS$FWD.RC)


# Run Cutadapt for 16S

dir.create("./logs/cutadapt")
for(i in seq_along(X16S$fnFs)) {
  message("Running cutadapt on: ", X16S$fnFs[i])
  logfile <- paste0("./logs/cutadapt/", X16S$sample[i], "_16S.log")
  
  system2("/home/rik/miniconda3/bin/cutadapt", 
          args = c(X16S$R1.flags, X16S$R2.flags, 
          "-n", 2, "--overlap", 10, "--discard-untrimmed", "--pair-filter=any","-e", 0.1,# -n 2 required to remove FWD and REV from reads
          "-o", X16S$fnFs.cut[i],"-p",X16S$fnRs.cut[i], # output files
          X16S$fnFs[i], X16S$fnRs[i]),  # input files
          stdout = logfile, stderr = logfile
  )
}


dir.create("./logs/cutadapt")
for(i in seq_along(XITS$fnFs)) {
  message("Running cutadapt on: ", XITS$fnFs[i])
  logfile <- paste0("./logs/cutadapt/", XITS$sample[i], "_ITS.log")
  
  system2("/home/rik/miniconda3/bin/cutadapt", 
          args = c(XITS$R1.flags, XITS$R2.flags, 
                   "-n", 2, "--overlap", 10, "--discard-untrimmed", "--pair-filter=any","-e", 0.1,# -n 2 required to remove FWD and REV from reads
                   "-o", XITS$fnFs.cut[i],"-p",XITS$fnRs.cut[i], # output files
                   XITS$fnFs[i], XITS$fnRs[i]),  # input files
          stdout = logfile, stderr = logfile
  )
}

# Let's add the numbers of reads left after primer trimming:
X16S$meta$Nreads_cut <- ShortRead::countFastq(X16S$fnRs.cut)$records
XITS$meta$Nreads_cut <- ShortRead::countFastq(XITS$fnRs.cut)$records

X16S$meta$Nnucl_cut <- ShortRead::countFastq(X16S$fnRs.cut)$nucleotides
XITS$meta$Nnucl_cut <- ShortRead::countFastq(XITS$fnRs.cut)$nucleotides


# Now we can check the trimmed reads for their presence of primer sequences:

# For 16S:
X16S$primersummary2 <- FUNS$generate_primer_summary(X16S, "fnFs.cut", "fnRs.cut", multithread = 12)
# For ITS:
XITS$primersummary2 <- FUNS$generate_primer_summary(XITS, "fnFs.cut", "fnRs.cut", multithread = 12)

# Relative to number of reads:
X16S$primersummary_rel2 <- X16S$primersummary2/ X16S$meta$Nreads_cut 
XITS$primersummary_rel2 <- XITS$primersummary2/ XITS$meta$Nreads_cut 


# Readthrough:
# (You can also recover this information from the cutadapt logs, but this is a nice way to double check and also to get the readthrough relative to the number of reads after primer trimming)
# For 16S:
X16S$readthrough_R1.2 <- X16S$primersummary2[,"R1_RP.rvcmp"]/X16S$meta$Nreads_cut 
X16S$readthrough_R2.2 <- X16S$primersummary2[,"R2_FP.rvcmp"]/X16S$meta$Nreads_cut
# For ITS:
XITS$readthrough_R1.2 <- XITS$primersummary2[,"R1_RP.rvcmp"]/XITS$meta$Nreads_cut 
XITS$readthrough_R2.2 <- XITS$primersummary2[,"R2_FP.rvcmp"]/XITS$meta$Nreads_cut












### 4. Filtering for base call quality
#######################################

### Now we move on to read filtering. We split the analysis between 16S and ITS
### Because for 16S (virtually constant insert size) we can use a fixed read length cut-off
### Given that we know the length of the primers, 
### we consider it better to just clip off the lenght of the primers, 
### rather than have them removed by cutadapt, which may miss them. 

### While for ITS, inserts are of variable size and hence it makes less sense to cut them off at a fixed length
### (see here: https://benjjneb.github.io/dada2/ITS_workflow.html)



### 4A. Filtering for base call quality: 16S
############################################
X16S$filtFs <- file.path(X16S$path.filtered, basename(X16S$fnFs.cut))
X16S$filtRs <- file.path(X16S$path.filtered, basename(X16S$fnRs.cut))

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
  matchIDs = T,
  verbose = TRUE
)

X16S$meta$Nreads_filter <- X16S$out[,2]

derepFs <- derepFastq(X16S$filtFs[-c(40,41)]) # exclude the mocks here
derepRs <- derepFastq(X16S$filtRs[-c(40,41)])

#"pseudo",PSEUDO_PREVALENCE = 1, PSEUDO_ABUNDANCE = 5, OMEGA_A = 1e-20)
dadaFs  <- dada(derepFs, err = NULL, selfConsist = T, MAX_CONSIST = 15, multithread = 20, pool = T, OMEGA_A = 1e-100, OMEGA_C = 1e-10, MIN_ABUNDANCE = 2) #"pseudo",PSEUDO_PREVALENCE = 1, PSEUDO_ABUNDANCE = 5, OMEGA_A = 1e-20)
dadaRs  <- dada(derepRs, err = NULL, selfConsist = T, MAX_CONSIST = 15, multithread = 20, pool = T, OMEGA_A = 1e-100, OMEGA_C = 1e-10, MIN_ABUNDANCE = 2) #"pseudo",PSEUDO_PREVALENCE = 1, PSEUDO_ABUNDANCE = 5, OMEGA_A = 1e-20)


X16S$errF <- FUNS$learnErrorsBalanced(X16S$filtFs, n_reads_per_sample = 12000, nbases = 1e8, multithread = 20, MAX_CONSIST = 15)
X16S$errR <- FUNS$learnErrorsBalanced(X16S$filtRs, n_reads_per_sample = 12000, nbases = 1e8, multithread = 20, MAX_CONSIST = 15)



























# Our first filtering step using the filterAndTrim function is really just getting rid of sequences containing N nucleotides, and PhiX. 
# About the multithread parameter: only use this if your computer allows parallel processing. 

# Generate a subdirectory that will contain the files
XITS$path_filtN <- file.path(XITS$path_intermediate, "filtN"); if(!dir.exists(XITS$path_filtN)) dir.create(XITS$path_filtN)
X16S$path_filtN <- file.path(X16S$path_intermediate, "filtN"); if(!dir.exists(X16S$path_filtN)) dir.create(X16S$path_filtN)

XITS$fnFs.filtN <- paste(XITS$path_filtN,"/",basename(XITS$fnFs),sep="")
XITS$fnRs.filtN <- paste(XITS$path_filtN,"/",basename(XITS$fnRs),sep="")
X16S$fnFs.filtN <- paste(X16S$path_filtN,"/",basename(X16S$fnFs),sep="")
X16S$fnRs.filtN <- paste(X16S$path_filtN,"/",basename(X16S$fnRs),sep="")

XITS$out <- filterAndTrim(XITS$fnFs, XITS$fnFs.filtN, XITS$fnRs, XITS$fnRs.filtN, maxN = 0, matchIDs=T, rm.phix=T, multithread = 12)
X16S$out <- filterAndTrim(X16S$fnFs, X16S$fnFs.filtN, X16S$fnRs, X16S$fnRs.filtN, maxN = 0, matchIDs=T, rm.phix=T, multithread = 12)








### 4. Filtering: removing primers and investigating readthrough
################################################################





# This function generates a table in which we count the number of reads that contain a primer in any direction. 
# The column names:
# R1 or R2: whether we are looking into the first (forward) or second (reverse) reads 
# (e.g. you can find forward primers in reverse reads in case of readthrough)
# FP or RP: whether we are looking for the forward primer or the reverse primer
# The direction of the primer: 
# Whether it's oriented the way it should (5` to 3`), complement, reversed (3` to 5`), or reverse complemented. 

# For ITS:
XITS$primersummary <- FUNS$generate_primer_summary(XITS, "fnFs.filtN", "fnRs.filtN", multithread = 12)
# For 16S:
X16S$primersummary <- FUNS$generate_primer_summary(X16S, "fnFs.filtN", "fnRs.filtN", multithread = 12)


# Relative to number of reads:
XITS$primersummary_rel <- XITS$primersummary/ XITS$META$Nreads 
X16S$primersummary_rel <- X16S$primersummary/ X16S$META$Nreads 

# For ITS
XITS$readthrough_R1 <- XITS$primersummary[,"R1_RP.rvcmp"]/XITS$META$Nreads 
XITS$readthrough_R2 <- XITS$primersummary[,"R2_FP.rvcmp"]/XITS$META$Nreads

# For 16S
X16S$readthrough_R1 <- X16S$primersummary[,"R1_RP.rvcmp"]/X16S$META$Nreads 
X16S$readthrough_R2 <- X16S$primersummary[,"R2_FP.rvcmp"]/X16S$META$Nreads

# Diagnostic plots for % of read through:

pdf("results/figures/02_ITS_readthrough_diagnostics.pdf", width = 7, height = 7)
plot(XITS$readthrough_R2, XITS$readthrough_R1, log = "xy", col = "white", cex = 0.75,
     main = "ITS: proportion of reads with full primer readthrough",xlab="proportion readthrough read 2", ylab="proportion readthrough read 1")
text(XITS$readthrough_R2, XITS$readthrough_R1, labels = XITS$sample.names, cex = 0.75)
dev.off()

pdf("results/figures/02_16S_readthrough_diagnostics.pdf", width = 7, height = 7)
plot(X16S$readthrough_R2 + 1e-05, X16S$readthrough_R1 + 1e-05, log = "xy", col = "white", cex = 0.75,
     main = "16S: proportion of reads with full primer readthrough",xlab="proportion readthrough read 2", ylab="proportion readthrough read 1")
text(X16S$readthrough_R2 + 1e-05, X16S$readthrough_R1 + 1e-05, labels = X16S$sample.names, cex = 0.75)
dev.off()


# For ITS:
XITS$primersummary.2 <- FUNS$generate_primer_summary(XITS, "fnFs.cut", "fnRs.cut")
# For 16S:
X16S$primersummary.2 <- FUNS$generate_primer_summary(X16S, "fnFs.cut", "fnRs.cut")








X16S$fnFs <- sort(list.files("/home/rik/fukushima/metabarcoding/data/raw/16S_reads/cutadapt", pattern = "1_bbmerge_corrected.fq.gz", full.names = TRUE))
#X16S$fnFs <- X16S$fnFs[match(X16S$sample.names, sapply(X16S$fnFs, FUNS$get.sample.name))]
#length(X16S$fnFs)
X16S$fnRs <- sort(list.files("/home/rik/fukushima/metabarcoding/data/raw/16S_reads/cutadapt", pattern = "2_bbmerge_corrected.fq.gz", full.names = TRUE))
#X16S$fnRs <- X16S$fnRs[match(X16S$sample.names, sapply(X16S$fnRs, FUNS$get.sample.name))]
#length(X16S$fnRs)

X16S$path_cut <- file.path(X16S$path, "cutadapt"); if(!dir.exists(X16S$path_cut)) dir.create(X16S$path_cut)


X16S$fnFs.cut <- file.path(X16S$path_cut, basename(X16S$fnFs))
X16S$fnRs.cut <- file.path(X16S$path_cut, basename(X16S$fnRs))

X16S$filtFs <- file.path(X16S$path_filtN, basename(X16S$fnFs.cut))
X16S$filtRs <- file.path(X16S$path_filtN, basename(X16S$fnRs.cut))


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
  matchIDs = T,
  verbose = TRUE
)

#X16S$errF    <- learnErrors(X16S$filtFs, multithread = 20, verbose = FALSE)
#X16S$errR    <- learnErrors(X16S$filtRs, multithread = 20, verbose = FALSE)

# Learn error rates (using all filtered reads)
X16S$errF <- FUNS$learnErrorsBalanced(X16S$filtFs, n_reads_per_sample = 12000, nbases = 1e8, multithread = 20, MAX_CONSIST = 15)
X16S$errR <- FUNS$learnErrorsBalanced(X16S$filtRs, n_reads_per_sample = 12000, nbases = 1e8, multithread = 20, MAX_CONSIST = 15)
#save.image("session20250926.Rdata")
save.image("session20260123a.Rdata")

#derepFs <- derepFastq(X16S$filtFs[1:49]) # exclude the mocks here
#derepRs <- derepFastq(X16S$filtRs[1:49])

derepFs <- derepFastq(X16S$filtFs[-21]) # exclude the mocks here
derepRs <- derepFastq(X16S$filtRs[-21])


#"pseudo",PSEUDO_PREVALENCE = 1, PSEUDO_ABUNDANCE = 5, OMEGA_A = 1e-20)
dadaFs  <- dada(derepFs, err = NULL, selfConsist = T, MAX_CONSIST = 15, multithread = 20, pool = T, OMEGA_A = 1e-100, OMEGA_C = 1e-100, MIN_ABUNDANCE = 2) #"pseudo",PSEUDO_PREVALENCE = 1, PSEUDO_ABUNDANCE = 5, OMEGA_A = 1e-20)
dadaRs  <- dada(derepRs, err = NULL, selfConsist = T, MAX_CONSIST = 15, multithread = 20, pool = T, OMEGA_A = 1e-100, OMEGA_C = 1e-100, MIN_ABUNDANCE = 2) #"pseudo",PSEUDO_PREVALENCE = 1, PSEUDO_ABUNDANCE = 5, OMEGA_A = 1e-20)
save.image("session20260212.Rdata")


#"pseudo",PSEUDO_PREVALENCE = 1, PSEUDO_ABUNDANCE = 5, OMEGA_A = 1e-20)
#dadaFs  <- dada(derepFs, err = X16S$errF, multithread = 20, pool = "pseudo", OMEGA_A = 1e-100, OMEGA_C = 1e-3)#"pseudo",PSEUDO_PREVALENCE = 1, PSEUDO_ABUNDANCE = 5, OMEGA_A = 1e-20)
#dadaRs  <- dada(derepRs, err = X16S$errR, multithread = 20, pool = "pseudo", OMEGA_A = 1e-100, OMEGA_C = 1e-3)#"pseudo",PSEUDO_PREVALENCE = 1, PSEUDO_ABUNDANCE = 5, OMEGA_A = 1e-20)


#save.image("session20250926.Rdata")
#save.image("session20260123a.Rdata")

mergers  <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=12, maxMismatch = 1, verbose = TRUE)
seqtab   <- makeSequenceTable(mergers)

#save.image("session20250926.Rdata")
#save.image("session20260123a.Rdata")

#seqtab <- seqtab[,colSums(seqtab)>1]
nochim  <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
#save.image("session20250926.Rdata")
#save.image("session20260123a.Rdata")
save.image("session_20260210.RData")



dadaFs  <- dada(derepFs, err = X16S$errF, multithread = 20, pool = T)
dadaRs  <- dada(derepRs, err = X16S$errR, multithread = 20, pool = T)
#save.image("session20250926.Rdata")
save.image("session20260123b.Rdata")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
seqtab  <- makeSequenceTable(mergers)
#save.image("session20250926.Rdata")
save.image("session20260123b.Rdata")

seqtab <- seqtab[,colSums(seqtab)>1]
nochim  <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
#save.image("session20250926.Rdata")
save.image("session20260123b.Rdata")



dadaFs  <- dada(derepFs, err = X16S$errF, multithread = 20, pool = T, OMEGA_A = 1e-20)
dadaRs  <- dada(derepRs, err = X16S$errR, multithread = 20, pool = T, OMEGA_A = 1e-20)
#save.image("session20250926.Rdata")
save.image("session20260123c.Rdata")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
seqtab  <- makeSequenceTable(mergers)
#save.image("session20250926.Rdata")
save.image("session20260123c.Rdata")

seqtab <- seqtab[,colSums(seqtab)>1]
nochim  <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
#save.image("session20250926.Rdata")
save.image("session20260123c.Rdata")






#table(colSums(seqtab)==1)
#FALSE  TRUE 
#58556 37814 

#table(colSums(seqtab)<=2)
#FALSE  TRUE 
#46786 49584 


seqtab  <- seqtab[,colSums(seqtab)>1]


group_ids <- X16S$META$soil_sample
X <- as.data.frame(seqtab) %>%
  mutate(group = group_ids) %>%
  group_by(group) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("group")


nochim  <- removeBimeraDenovo(X, method = "pooled", multithread = TRUE)


# Extract sequences
seqs1 <- DNAStringSet(colnames(X))
seqs2 <- readDNAStringSet("data/mockbacteria_original.fasta")# Set k-mer size
k <- 7  # You can tweak this
# Compute k-mer frequency matrix

kmer1 <- oligonucleotideFrequency(seqs1, width = k, as.prob = TRUE)
kmer2 <- oligonucleotideFrequency(seqs2, width = k, as.prob = TRUE)
dist_mat <- as.matrix(proxy::dist(kmer1, kmer2, method = "euclidean"))

kmer_mat <- oligonucleotideFrequency(seqs, width = k, as.prob = TRUE)
# Compute Euclidean distance between k-mer profiles
dist_mat <- dist(kmer_mat, method = "euclidean")

dist_mat <- as.matrix(proxy::dist(kmer1, kmer2, method = "euclidean"))

save.image("session20250925.Rdata")


# Step 1: get grouping variable
group_ids <- X16S$META$soil_sample

# Step 2: aggregate counts by soil_sample
X <- as.data.frame(nochim) %>%
  mutate(group = group_ids) %>%
  group_by(group) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("group")




asv_seqs <- colnames(X)
asv_fasta <- setNames(asv_seqs, paste0("ASV", seq_along(asv_seqs)))
writeXStringSet(DNAStringSet(asv_fasta), filepath = "asvs.fasta")


# ~/cdhit/cd-hit-est -i asvs.fasta -o clustered_otus.fasta -c 0.99 -n 10 -T 10 -M 12000
# mafft --localpair toalign.fasta > aligned.fasta



parse_cdhit_clstr <- function(clstr_file) {
     lines <- readLines(clstr_file)
     cluster_id <- NA
     asv_to_otu <- list()
     
     for (line in lines) {
           if (startsWith(line, ">Cluster")) {
               cluster_id <- sub(">Cluster ", "OTU", line)
             } else {
                 # Extract ASV name from line (e.g. ">ASV42... *")
                   asv <- sub('.*>(ASV\\d+).*', '\\1', line)
                   
                       # If representative (has "*"), this is the OTU name
                       if (grepl("\\*$", line)) {
                           rep <- asv
                       }
                   asv_to_otu[[asv]] <- rep
                 }
       }
   
     return(unlist(asv_to_otu))
}
asv_to_otu <- parse_cdhit_clstr("clustered_otus.fasta.clstr")
names(which(asv_to_otu == "ASV58"))
asv_fasta[names(which(asv_to_otu == "ASV58"))]


# Extract sequences
seqs <- DNAStringSet(colnames(X))

# Set k-mer size
k <- 7  # You can tweak this

# Compute k-mer frequency matrix
kmer_mat <- oligonucleotideFrequency(seqs, width = k, as.prob = TRUE)

# Compute Euclidean distance between k-mer profiles
dist_mat <- dist(kmer_mat, method = "euclidean")

# Convert to a matrix if you like
dist_matrix    <- as.matrix(dist_mat)
dist_mat_upper <- dist_matrix
dist_mat_upper[lower.tri(dist_mat_upper)] <- 1   # or 0 if you prefer
dist_mat_lower <- dist_matrix
dist_mat_lower[upper.tri(dist_mat_lower)] <- 1







dim(dist_mat_upper)



XITS$out <- filterAndTrim(XITS$fnFs, XITS$fnFs.filtN, XITS$fnRs, XITS$fnRs.filtN, maxN = 0, matchIDs=T, rm.phix=T, multithread = 12)




library(dada2)
library(data.table)
library(Biostrings)

# Parameter grid
truncLen_opts <- list(c(230, 210), c(230, 220))
truncQ_opts   <- c(5, 7, 9, 11)
maxEE_opts <- list(
  c(0.5,1), c(0.5,2), c(1,1), c(1,2), c(1,3),
  c(2,2), c(2,3), c(2,4), c(3,3), c(3,4)
)

# Sample selection
samples <- c("A6T", "T5T", "O5T", "F1T", "HT")#,"mock1","mock2")
basedir <- "data/intermediate/16S_reads/cutadapt"
fwd_files <- file.path(basedir, paste0(samples, "_16S_1.fq.gz"))
rev_files <- file.path(basedir, paste0(samples, "_16S_2.fq.gz"))
names(fwd_files) <- samples
names(rev_files) <- samples

results <- list()
i <- 1

# Loop over parameter combinations
for (truncLen in truncLen_opts) {
  for (truncQ in truncQ_opts) {
    for (maxEE in maxEE_opts) {
      
      fw_ee <- maxEE[1]
      rv_ee <- maxEE[2]
      
      message(">>> Running: truncLen = ", paste(truncLen, collapse = ","),
              " | truncQ = ", truncQ,
              " | maxEE = ", fw_ee, "/", rv_ee)
      
      tmpdir <- tempfile("dada_test_")
      dir.create(tmpdir)
      
      filtFs <- file.path(tmpdir, paste0(samples, "_filtF.fastq.gz"))
      filtRs <- file.path(tmpdir, paste0(samples, "_filtR.fastq.gz"))
      
      
      
      out <- suppressWarnings(filterAndTrim(
        fwd = fwd_files,
        filt = filtFs,
        rev = rev_files,
        filt.rev = filtRs,
        truncLen = truncLen,
        truncQ = truncQ,
        maxEE = c(fw_ee, rv_ee),
        maxN = 0,
        rm.phix = TRUE,
        compress = TRUE,
        multithread = 20
      ))
      
      reads_in  <- out[, "reads.in"]
      reads_out <- out[, "reads.out"]
      pct_retained <- round(reads_out / reads_in * 100, 1)
      
      if (!any(c(file.exists(filtFs),file.exists(filtRs))==F)) {
        
      errF    <- suppressWarnings(learnErrors(filtFs, multithread = 20, verbose = FALSE))
      errR    <- suppressWarnings(learnErrors(filtRs, multithread = 20, verbose = FALSE))
      derepFs <- derepFastq(filtFs)
      derepRs <- derepFastq(filtRs)
      dadaFs  <- dada(derepFs, err = errF, multithread = 20, pool = F)
      dadaRs  <- dada(derepRs, err = errR, multithread = 20, pool = F)
      mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = FALSE)
      seqtab  <- makeSequenceTable(mergers)
      
      n_asvs      <- ncol(seqtab)
      asv_sums    <- colSums(seqtab)
      n_asvs_1    <- sum(asv_sums == 1)
      n_asvs_5    <- sum(asv_sums < 5)
      n_asvs_10   <- sum(asv_sums < 10)
      merge_rates <- unlist(lapply(mergers, function(x) sum(x$abundance))) / reads_out
      merge_pct   <- round(mean(merge_rates, na.rm = TRUE) * 100, 1)
      
      results[[i]] <- data.frame(
        truncLen_F = truncLen[1],
        truncLen_R = truncLen[2],
        truncQ = truncQ,
        maxEE_F = fw_ee,
        maxEE_R = rv_ee,
        mean_reads_in = mean(reads_in),
        mean_reads_out = mean(reads_out),
        mean_pct_retained = mean(pct_retained),
        mean_pct_merged = merge_pct,
        n_ASVs = n_asvs,
        ASVs_eq1 = n_asvs_1,
        ASVs_lt5 = n_asvs_5,
        ASVs_lt10 = n_asvs_10
      )
      }else{
        warning("Filtered files not found. Skipping.")

        results[[i]] <- data.frame(
          truncLen_F = truncLen[1],
          truncLen_R = truncLen[2],
          truncQ = truncQ,
          maxEE_F = fw_ee,
          maxEE_R = rv_ee,
          mean_reads_in = mean(reads_in),
          mean_reads_out = mean(reads_out),
          mean_pct_retained = mean(pct_retained),
          mean_pct_merged = NA,
          n_ASVs = NA,
          ASVs_eq1 = NA,
          ASVs_lt5 = NA,
          ASVs_lt10 = NA
        )  
        
      }
      
      i <- i + 1
      unlink(tmpdir, recursive = TRUE)
    }
  }
}

summary_df <- rbindlist(results)
data.table::fwrite(summary_df, file = "results/figures/03_param_tuning_summary2.csv")
print(summary_df)



plot(x$mean_pct_retained,col=x$truncQ,axes=1,xlab="",xaxt="n",ylim = c(c(-3,1)+range(x$mean_pct_retained)),main=c("Mean percentage retained after first filtering step"),ylab = "mean percentage retained")
abline(v=c(40,80))
text(x=c(12,52,92),y=c(42,42,42), labels= c("truncLen 240 220","truncLen 230 210", "truncLen 230 220"),cex=0.6)
legend("topright",col=truncQ_opts,legend=truncQ_opts, title = "truncQ", pch = 20)

plot(x$mean_pct_merged,col=x$truncQ,axes=1,xlab="",xaxt="n",ylim = c(c(-1,1)+range(x$mean_pct_merged)),main=c("Mean percentage merged"),ylab = "mean percentage merged")
abline(v=c(40,80))
text(x=c(12,52,92),y=c(73,73,73), labels= c("truncLen 240 220","truncLen 230 210", "truncLen 230 220"),cex=0.6)
legend("topleft",col=truncQ_opts,legend=truncQ_opts, title = "truncQ", pch = 20)

plot(x$n_ASVs,col=x$truncQ,axes=1,xlab="",xaxt="n",ylim = c(c(-500,200)+range(x$n_ASVs)),main=c("Total number of ASVs"),ylab = "# ASV's")
abline(v=c(40,80))
text(x=c(12,52,92),y=c(5000,5000,5000), labels= c("truncLen 240 220","truncLen 230 210", "truncLen 230 220"),cex=0.6)
legend("topleft",col=truncQ_opts,legend=truncQ_opts, title = "truncQ", pch = 20)


plot(x$n_ASVs,col=x$truncQ,axes=1,xlab="",xaxt="n",ylim = c(0,12000),main=c("Total number of ASVs"),ylab = "# ASV's")
abline(v=c(40,80))
points(x$ASVs_lt10,col=x$truncQ,pch="*")
points(x$ASVs_eq1,col=x$truncQ,pch="°")
text(x=c(12,52,92),y=c(100,100,100), labels= c("truncLen 240 220","truncLen 230 210", "truncLen 230 220"),cex=0.6)
legend("topleft",col=truncQ_opts,legend=truncQ_opts, title = "truncQ", pch = 20)



plot(x$ASVs_lt10,col=x$truncQ,axes=1,xlab="",xaxt="n",ylim = c(c(-500,200)+range(x$ASVs_lt10)),main=c("Number of ASVs with at least 10 counts"),ylab = "# ASV's >10")
abline(v=c(40,80))
text(x=c(12,52,92),y=c(2100,2100,2100), labels= c("truncLen 240 220","truncLen 230 210", "truncLen 230 220"),cex=0.6)
legend("topleft",col=truncQ_opts,legend=truncQ_opts, title = "truncQ", pch = 20)


plot(x$n_ASVs / x$mean_reads_out*x$mean_pct_merged,col=x$truncQ,axes=1,xlab="",xaxt="n",ylim=c(0,20),main="Percentage of ASV's vs total merged reads")
text(x=c(15,55,95),y=c(2,2,2), labels= c("truncLen 240 220","truncLen 230 210", "truncLen 230 220"),cex=0.6)
abline(v=c(40,80))
legend("topleft",col=truncQ_opts,legend=truncQ_opts, title = "truncQ", pch = 20)



















### 5A. Filtering the 16S reads for base call quality
#####################################################

X16S$filtFs <- file.path(X16S$path_intermediate, basename(X16S$fnFs.cut))
X16S$filtRs <- file.path(X16S$path_intermediate, basename(X16S$fnRs.cut))
X16S$filtFs_fw_only <- file.path("./data/intermediate/16S_reads/filtered_fw_only", basename(X16S$fnFs.cut))  # Adjust path
X16S$filtFs_fw_only <- file.path("./data/intermediate/16S_reads/filtered_fw_only120", basename(X16S$fnFs.cut))  # Adjust path
X16S$filtFs_fw_only <- file.path("./data/intermediate/16S_reads/filtered_fw_only60", basename(X16S$fnFs.cut))  # Adjust path


for(TL in c(160,170,180,190,200))
    {
      for(TQ in c(15,20,25))
      {
        X16S$filtFs_fw_only <- file.path("./data/intermediate/16S_reads/filtered_fw_tmp", basename(X16S$fnFs.cut))  # Adjust path
        filt_out <- filterAndTrim(
        fwd = X16S$fnFs.filtN,
        filt = X16S$filtFs_fw_only,
        trimLeft = 20,
        truncLen = (TL+20),
        truncQ = TQ,
        maxEE = 0.05,
        minLen = (TL-10),
        multithread = TRUE
        )
        
        X16S$derepF_trimmed   <- derepFastq(X16S$filtFs_fw_only)
        X16S$errF_trimmed     <- learnErrors(X16S$filtFs_fw_only, multithread = TRUE)
        
        X16S$dadaF_trimmed    <- dada(X16S$derepF_trimmed, err = X16S$errF_trimmed, pool=TRUE, multithread = TRUE)
        X16S$seqtabF_trimmed  <- makeSequenceTable(X16S$dadaF_trimmed)
        X16S$seqtabF_trimmed.nochim1 <- removeBimeraDenovo(X16S$seqtabF_trimmed, method = "consensus", multithread = TRUE)
        X16S$seqtabF_trimmed.nochim2 <- removeBimeraDenovo(X16S$seqtabF_trimmed, method = "pooled", multithread = TRUE)
        
        out_pool <- list("derep" = X16S$derepF_trimmed, "seqtab" = X16S$seqtabF_trimmed,  "nochim_consensus" = X16S$seqtabF_trimmed.nochim1, "nochim_pooled" = X16S$seqtabF_trimmed.nochim2)
        
        X16S$dadaF_trimmed    <- dada(X16S$derepF_trimmed, err = X16S$errF_trimmed, pool=FALSE, multithread = TRUE)
        X16S$seqtabF_trimmed  <- makeSequenceTable(X16S$dadaF_trimmed)
        X16S$seqtabF_trimmed.nochim1 <- removeBimeraDenovo(X16S$seqtabF_trimmed, method = "consensus", multithread = TRUE)
        X16S$seqtabF_trimmed.nochim2 <- removeBimeraDenovo(X16S$seqtabF_trimmed, method = "pooled", multithread = TRUE)
        
        out_nopool <- list("derep" = X16S$derepF_trimmed, "seqtab" = X16S$seqtabF_trimmed,  "nochim_consensus" = X16S$seqtabF_trimmed.nochim1, "nochim_pooled" = X16S$seqtabF_trimmed.nochim2)
        
        OUT <- list("pool" = out_pool, "nopool" = out_nopool)
  
        outname <- paste("out",TL, TQ, sep="_")
        
        assign(outname, OUT)
      }
    }

    
ssum <- function(outobj)
{    
NASV1                <- ncol(outobj[[2]])  
retained_reads_perc1 <- summary(round(rowSums(outobj[[2]])/X16S$META$Nreads*100))    
cor_duplicates1      <- round(100*cor(as.vector(log(outobj[[2]][which(duplicated(X16S$META$soil_sample)),]+1)),as.vector(log(outobj[[2]][which(duplicated(X16S$META$soil_sample))-1,]+1))))
#plot(as.vector(log(out_100_20$pool[[2]][which(duplicated(X16S$META$soil_sample)),]+1)),as.vector(log(out_100_20$pool[[2]][which(duplicated(X16S$META$soil_sample))-1,]+1)))
cor_control_site1    <- round(100*cor(as.vector(log(outobj[[2]][c(2,5,9,15,21,25,29,33,37,45,43,23),]+1)),as.vector(log(outobj[[2]][c(12,10,6,14,18,24,16,28,36,46,42,34),]+1))))
cor_control_random1  <- round(100*cor(as.vector(log(outobj[[2]][c(1,4,8,16,20,24,28,32,36,40,44,48),]+1)),as.vector(log(outobj[[2]][c(21,35,37,39,7,15,43,3,1,27,29,25),]+1))))

NASV2                <- ncol(outobj[[3]])  
retained_reads_perc2 <- summary(round(rowSums(outobj[[3]])/X16S$META$Nreads*100))    
cor_duplicates2      <- round(100*cor(as.vector(log(outobj[[3]][which(duplicated(X16S$META$soil_sample)),]+1)),as.vector(log(outobj[[3]][which(duplicated(X16S$META$soil_sample))-1,]+1))))
cor_control_site2    <- round(100*cor(as.vector(log(outobj[[3]][c(2,5,9,15,21,25,29,33,37,45,43,23),]+1)),as.vector(log(outobj[[3]][c(12,10,6,14,18,24,16,28,36,46,42,34),]+1))))
cor_control_random2  <- round(100*cor(as.vector(log(outobj[[3]][c(1,4,8,16,20,24,28,32,36,40,44,48),]+1)),as.vector(log(outobj[[3]][c(21,35,37,39,7,15,43,3,1,27,29,25),]+1))))

NASV3                <- ncol(outobj[[4]])  
retained_reads_perc3 <- summary(round(rowSums(outobj[[4]])/X16S$META$Nreads*100))    
cor_duplicates3      <- round(100*cor(as.vector(log(outobj[[4]][which(duplicated(X16S$META$soil_sample)),]+1)),as.vector(log(outobj[[4]][which(duplicated(X16S$META$soil_sample))-1,]+1))))
cor_control_site3    <- round(100*cor(as.vector(log(outobj[[4]][c(2,5,9,15,21,25,29,33,37,45,43,23),]+1)),as.vector(log(outobj[[4]][c(12,10,6,14,18,24,16,28,36,46,42,34),]+1))))
cor_control_random3  <- round(100*cor(as.vector(log(outobj[[4]][c(1,4,8,16,20,24,28,32,36,40,44,48),]+1)),as.vector(log(outobj[[4]][c(21,35,37,39,7,15,43,3,1,27,29,25),]+1))))

OUT <- as.data.frame(rbind.data.frame(c(NASV1,cor_duplicates1,cor_control_site1,cor_control_random1,retained_reads_perc1),c(NASV2,cor_duplicates2,cor_control_site2,cor_control_random2,retained_reads_perc2),c(NASV3,cor_duplicates3,cor_control_site3,cor_control_random3,retained_reads_perc3)))
OUT <- cbind.data.frame(as.data.frame(c("seqtab","nochim_consensus","nochim_pooled")),round(OUT))
colnames(OUT) <- c("bimeras","n_ASV","cor_dup","cor_site","cor_ctr","min","q1","median","mean","q3","max")

return(OUT)
}


bar<-ssum(out_100_20$nopool)    


SSUM <- function(TL,TQ)
{
  bar <- paste("out",TL,TQ,sep="_")
  foo <- get(bar)
  x_pool   <- cbind.data.frame(rep(TL,3),rep(TQ,3),rep("pool",3),ssum(foo$pool))
  colnames(x_pool) <- c("TL","TQ","dada","bimeras","n_ASV","cor_dup","cor_site","cor_ctr","min","q1","median","mean","q3","max")
  x_nopool <- cbind.data.frame(rep(TL,3),rep(TQ,3),rep("nopool",3),ssum(foo$nopool))
  colnames(x_nopool) <- c("TL","TQ","dada","bimeras","n_ASV","cor_dup","cor_site","cor_ctr","min","q1","median","mean","q3","max")
  
  OUT <- rbind.data.frame(x_pool,x_nopool,make.row.names=F)
  return(OUT)
}

SSUM(100,20)


LIST <- list()
counter = i
for(TL in c(50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200))
{
  for(TQ in c(15,20,25))
  {
    LIST[[i]] <- SSUM(TL,TQ)
    i=i+1
  }
}

outcomes <- do.call("rbind.data.frame",LIST)







   
    
    

filt_out <- filterAndTrim(
  fwd = X16S$fnFs.filtN,
  filt = X16S$filtFs_fw_only,
  trimLeft = 20,
  truncLen = 80,
  truncQ = 20,
  maxEE = 0.05,
  minLen = 50,
  multithread = TRUE
)

# 2. Dereplicate
X16S$derepF_trimmed <- derepFastq(X16S$filtFs_fw_only)
#names(X16S$derepF_trimmed) <- names(X16S$fqFs)

# 3. Learn error rates (optional: plot it)
X16S$errF_trimmed <- learnErrors(X16S$filtFs_fw_only, multithread = TRUE)

# 4. Denoise
X16S$dadaF_trimmed <- dada(X16S$derepF_trimmed, err = X16S$errF_trimmed, pool=TRUE, multithread = TRUE)

# 5. Build sequence table
X16S$seqtabF_trimmed <- makeSequenceTable(X16S$dadaF_trimmed)

# 6. Remove bimeras
X16S$seqtabF_trimmed.nochim1 <- removeBimeraDenovo(X16S$seqtabF_trimmed, method = "consensus", multithread = TRUE)
X16S$seqtabF_trimmed.nochim2 <- removeBimeraDenovo(X16S$seqtabF_trimmed, method = "pooled", multithread = TRUE)



X16S$out <- filterAndTrim(
  fwd = X16S$fnFs.cut,
  filt = X16S$filtFs,
  rev = X16S$fnRs.cut,
  filt.rev = X16S$filtRs,
  truncLen = c(270, 190),
  truncQ = 7,
  maxEE = c(0.5, 1),
  maxN = 0,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = 20,
  matchIDs = T,
  verbose = TRUE
)


X16S$out <- filterAndTrim(
  fwd = X16S$fnFs.cut,
  filt = X16S$filtFs,
  rev = X16S$fnRs.cut,
  filt.rev = X16S$filtRs,
  truncQ = 8,
  maxEE = c(0.5, 0.5),
  maxN = 0,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = 20,
  verbose = TRUE,
  trimRight = 30,
  minLen = 170,
  matchIDs = T
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
X16S$dadaFs <- dada(X16S$derepFs, err = X16S$errF, multithread = 20, pool = F)
X16S$dadaRs <- dada(X16S$derepRs, err = X16S$errR, multithread = 20, pool = F)

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


