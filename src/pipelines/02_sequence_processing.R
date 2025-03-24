# @ Rik Verdonck 20250121
# Purpose: Preprocess raw sequencing data (DADA2 pipeline).

# Quality filtering of reads (trimming, filtering, etc.).
# Learn error rates and perform sequence inference (DADA2).
# Generate and save ASV tables for 16S and ITS.
# Assign taxonomy to ASVs (e.g., using SILVA for 16S and UNITE for ITS).

library("tidyverse")
library("data.table")
library("dada2")
library("gridExtra")  # For combining plots
library("doParallel")



# We will for now save everything that has to do with ITS in the list object XITS, and everything 16S in the list object X16S
XITS <- list()
X16S <- list()

# Some small custom functions:
source("./src/functions/functions.R")

# The paths to the data:
XITS$path <-"./data/raw/ITS_reads"
X16S$path <-"./data/raw/16S_reads"

# The paths to intermediates:
XITS$path_intermediate <-"./data/intermediate/ITS_reads"
X16S$path_intermediate <-"./data/intermediate/16S_reads"

# The entire list of libraries present in the paths:
XITS$fnFs <- sort(list.files(XITS$path, pattern = "_1.fq.gz", full.names = TRUE))
length(XITS$fnFs)
XITS$fnRs <- sort(list.files(XITS$path, pattern = "_2.fq.gz", full.names = TRUE))
length(XITS$fnRs)

X16S$fnFs <- sort(list.files(X16S$path, pattern = "_1.fq.gz", full.names = TRUE))
length(X16S$fnFs)
X16S$fnRs <- sort(list.files(X16S$path, pattern = "_2.fq.gz", full.names = TRUE))
length(X16S$fnRs)


# Extract all sample names:
XITS$sample.names    <- unname(sapply(XITS$fnFs, FUNS$get.sample.name))
X16S$sample.names    <- unname(sapply(X16S$fnFs, FUNS$get.sample.name))
sample.names         <- unique(c(XITS$sample.names,X16S$sample.names))

# Read in metadata and double check whether sample names correspond
META <- fread("./data/intermediate/cleaned_metadata.csv")
rownames(META) <- META$DNA_extract

# Remember some samples were discarded earlier. 
setdiff(sample.names,META$DNA_extract)
setdiff(META$DNA_extract,sample.names)





# Here I define 5 colors that we will use throughout for site, used for the levels 
# "Arakawa", "Futaba", "Okuma", "Tsushima" and "University" 
# I re-level META to put them in a more logical order as well. 
META$sampling_site   <- factor(META$sampling_site,levels=c("University", "Arakawa", "Tsushima", "Okuma", "Futaba"))
sitecolors  <- c("#ADD8E6","#2E8B57","#FF7F50","#FF91A4","#800020")


# Make sure the levels and order of META and samples in X objects are the same:

XITS$META <- META %>%
  filter(!is.na(lib_ITS_forward) & !is.na(lib_ITS_reverse)) %>%  # Keep rows where ITS libraries are available
  select(
    DNA_extract, soil_sample, sampling_data, sampling_site, soil_type,
    Cs_137, Cs_134, "lib_forward"=lib_ITS_forward, "lib_reverse"=lib_ITS_reverse,
    "primer_forward"=primer_ITS_forward, "primer_reverse"=primer_ITS_reverse, mock, "Nreads"=Nreads_ITS
  )

its_indices       <- match(XITS$META$DNA_extract, XITS$sample.names)
XITS$fnFs         <- XITS$fnFs[its_indices]     
XITS$fnRs         <- XITS$fnRs[its_indices]     
XITS$sample.names <- XITS$sample.names[its_indices]
  

X16S$META   <- META[match(X16S$sample.names, META$DNA_extract), ]
X16S$META <- META %>%
  filter(!is.na(lib_16S_forward) & !is.na(lib_16S_reverse)) %>%  # Keep rows where 16S libraries are available
  select(
    DNA_extract, soil_sample, sampling_data, sampling_site, soil_type,
    Cs_137, Cs_134, "lib_forward"=lib_16S_forward, "lib_reverse"=lib_16S_reverse,
    "primer_forward"=primer_16S_forward, "primer_reverse"=primer_16S_reverse, mock, "Nreads"=Nreads_16S
  )

x16s_indices      <- match(XITS$META$DNA_extract, X16S$sample.names)
X16S$fnFs         <- X16S$fnFs[x16s_indices]     
X16S$fnRs         <- X16S$fnRs[x16s_indices]     
X16S$sample.names <- X16S$sample.names[x16s_indices]


# Output directory
output_dir <- "results/figures/qualityprofiles"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save paired quality profiles for ITS
FUNS$save_paired_quality_profiles(XITS$fnFs, XITS$fnRs, "ITS")
# Save paired quality profiles for 16S
FUNS$save_paired_quality_profiles(X16S$fnFs, X16S$fnRs, "16S")




# Our first filtering step using the filterAndTrim function is really just getting rid of sequences containing N nucleotides, and PhiX. 
# About the multithread parameter: only use this if your computer allows parallel processing. 

XITS$fnFs.filtN <- paste(XITS$path_intermediate,"/",XITS$sample.names,"_ITS_1.filtN.fq.gz",sep="")
XITS$fnRs.filtN <- paste(XITS$path_intermediate,"/",XITS$sample.names,"_ITS_2.filtN.fq.gz",sep="")

X16S$fnFs.filtN <- paste(X16S$path_intermediate,"/",X16S$sample.names,"_16S_1.filtN.fq.gz",sep="")
X16S$fnRs.filtN <- paste(X16S$path_intermediate,"/",X16S$sample.names,"_16S_2.filtN.fq.gz",sep="")


XITS$out <- filterAndTrim(XITS$fnFs, XITS$fnFs.filtN, XITS$fnRs, XITS$fnRs.filtN, maxN = 0, matchIDs=T, rm.phix=T, multithread = 12)
X16S$out <- filterAndTrim(X16S$fnFs, X16S$fnFs.filtN, X16S$fnRs, X16S$fnRs.filtN, maxN = 0, matchIDs=T, rm.phix=T, multithread = 12)


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

plot(XITS$proportion_revcomp_fw,XITS$proportion_revcomp_rv,log="xy",col="white",cex=0.75, main = "ITS: proportion of reads with readthrough of the entire primer")
text(XITS$proportion_revcomp_fw,XITS$proportion_revcomp_rv,labels=XITS$sample.names,cex=0.75)

plot(X16S$proportion_revcomp_fw + 1e-06 ,X16S$proportion_revcomp_rv + 1e-06 , log="xy",col="white",cex=0.75, main = "16S: proportion of reads with readthrough of the entire primer")
text(X16S$proportion_revcomp_fw + 1e-06 ,X16S$proportion_revcomp_rv + 1e-06 , labels=X16S$sample.names,cex=0.75)


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




### Now we move on to read filtering. We split the analysis between 16S and ITS
### Because for 16S (virtually constant insert size) we can use a fixed read length cut-off
### While for ITS, inserts are of variable size and hence it makes less sense to cut them off at a fixed length
### (see here: https://benjjneb.github.io/dada2/ITS_workflow.html)

### So, first for 16S: 

filtFs <- file.path(X16S$path_intermediate, basename(X16S$fnFs.cut))
filtRs <- file.path(X16S$path_intermediate, basename(X16S$fnRs.cut))


X16S$out <- filterAndTrim(
  fwd = X16S$fnFs.cut,
  filt = filtFs,
  rev = X16S$fnRs.cut,
  filt.rev = filtRs,
  truncLen = c(250, 230),
  truncQ = 5,
  maxEE = c(2, 4),
  maxN = 0,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = 20,
  verbose = TRUE
)


# 2. Learn error rates (using all filtered reads)
errF <- FUNS$learnErrorsBalanced(filtFs, n_reads_per_sample = 10000, nbases = 1e8, multithread = 20, MAX_CONSIST = 20)
errR <- FUNS$learnErrorsBalanced(filtRs, n_reads_per_sample = 10000, nbases = 1e8, multithread = 20, MAX_CONSIST = 20)


# 3. Dereplicate
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)


# 4. Sample inference (DADA)
dadaFs <- dada(derepFs, err = errF, multithread = 20, pool = "pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = 20, pool = "pseudo")


# 5. Merging paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE, maxMismatch=3)

# 6. Construct ASV table
seqtab <- makeSequenceTable(mergers)

# 7. Chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = 12, verbose = TRUE)

# Save outputs to X16S
X16S$filtered.F <- filtFs ; rm(filtFs)
X16S$filtered.R <- filtRs ; rm(filtRs)
X16S$errF <- errF ; rm(errF)
X16S$errR <- errR ; rm(errR)
X16S$dadaFs <- dadaFs ; rm(dadaFs)
X16S$dadaRs <- dadaRs ; rm(dadaRs)
X16S$mergers <- mergers ; rm(mergers)
X16S$seqtab <- seqtab ; rm(seqtab)
X16S$seqtab.nochim <- seqtab.nochim ; rm(seqtab.nochim)

# Optional: summary stats
cat("Reads after filtering:\n")
print(rowSums(X16S$out))

cat("Percentage of non-chimeric reads:\n")
print(sum(seqtab.nochim) / sum(seqtab))



