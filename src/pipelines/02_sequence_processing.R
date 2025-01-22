# @ Rik Verdonck 20250121
# Purpose: Preprocess raw sequencing data (DADA2 pipeline).

# Quality filtering of reads (trimming, filtering, etc.).
# Learn error rates and perform sequence inference (DADA2).
# Generate and save ASV tables for 16S and ITS.
# Assign taxonomy to ASVs (e.g., using SILVA for 16S and UNITE for ITS).


library("dada2")



# We will for now save everything that has to do with ITS in the list object XITS, and everything 16S in the list object X16S
XITS <- list()
X16S <- list()

# Some small custom functions:
FUNS <- list()
source("./src/functions/functions.R")

# The paths to the data:
XITS$path <-"./data/raw/ITS_reads"
X16S$path <-"./data/raw/16S_reads"


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
XITS$sample.names[XITS$sample.names == "mock"] <- c("mockITS","mockITS1")
X16S$sample.names    <- unname(sapply(X16S$fnFs, FUNS$get.sample.name))
X16S$sample.names[X16S$sample.names == "mock"] <- c("mock16S","mock16S1")
sample.names         <- unique(c(XITS$sample.names,X16S$sample.names))

# Read in metadata and double check whether sample names correspond
META <- fread("./data/intermediate/cleaned_metadata.csv")
rownames(META) <- META$DNA_extract

setdiff(sample.names,META$DNA_extract)
setdiff(META$DNA_extract,sample.names)


# Here I define 5 colors that we will use throughout for site, used for the levels 
# "Arakawa", "Futaba", "Okuma", "Tsushima" and "University" 
# I re-level META to put them in a more logical order as well. 
META$site <- factor(META$site,levels=c("University", "Arakawa", "Tsushima", "Okuma", "Futaba"))
sitecolors     <- c("#ADD8E6","#2E8B57","#FF7F50","#FF91A4","#800020")

XITS$META <- META[match(XITS$sample.names, META$library), ]
X16S$META <- META[match(X16S$sample.names, META$library), ]








plotQualityProfile(c(XITS$fnFs[16],XITS$fnRs[16]))


