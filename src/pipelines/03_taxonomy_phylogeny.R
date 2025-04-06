# @ Rik Verdonck 20250403
# Purpose: 
# * Read in the TreeSE objects and assign taxonomy
# * Check the taxonomy assingments of the mock data as a quality check
# * Generate phylogenies for all our ASV's
# * Remove plastid-derived sequences if necessary
# * Cluster sequences into OTU's
# * Output a first series of visualisations of taxonomic composition

library("dada2")
library("phyloseq")
library("mia")
library("TreeSummarizedExperiment")

### 1. Reading in data, setting up TreeSE objects
########################################################

# Some small custom functions:
source("./src/functions/functions.R")

# Data:
load("results/objects/02_TreeSE_objects.RData")

# Define TSE objects for more lenient taxonomy assignment as well. 
TSE_16S_lenient        <- TSE_16S     
TSE_16S_mock_lenient   <- TSE_16S_mock
TSE_ITS_lenient        <- TSE_ITS     
TSE_ITS_mock_lenient   <- TSE_ITS_mock


# ITS mock taxonomy
TSE_ITS_mock_tax <- assignTaxonomy(
  refFasta = "db/sh_general_release_dynamic_s_all_19.02.2025.fasta",
  seqs = rowData(TSE_ITS_mock)$sequence,
  multithread = 16,
  tryRC = TRUE,
  verbose = TRUE
)
rowData(TSE_ITS_mock) <- cbind(rowData(TSE_ITS_mock), TSE_ITS_mock_tax)
rm(TSE_ITS_mock_tax)

# ITS taxonomy
TSE_ITS_tax <- assignTaxonomy(
  refFasta = "db/sh_general_release_dynamic_s_all_19.02.2025.fasta",
  seqs = rowData(TSE_ITS)$sequence,
  multithread = 16,
  tryRC = TRUE,
  verbose = TRUE
)
rowData(TSE_ITS) <- cbind(rowData(TSE_ITS), TSE_ITS_tax)
rm(TSE_ITS_tax)

# 16S mock taxonomy
TSE_16S_mock_tax <- assignTaxonomy(
  refFasta = "db/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
  seqs = rowData(TSE_16S_mock)$sequence,
  multithread =16,
  tryRC = TRUE,
  verbose = TRUE
)
rowData(TSE_16S_mock) <- cbind(rowData(TSE_16S_mock), TSE_16S_mock_tax)
rm(TSE_16S_mock_tax)

# 16S taxonomy
TSE_16S_tax <- assignTaxonomy(
  refFasta = "db/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
  seqs = rowData(TSE_16S)$sequence,
  multithread = 16,
  tryRC = TRUE,
  verbose = TRUE
)
rowData(TSE_16S) <- cbind(rowData(TSE_16S), TSE_16S_tax)
rm(TSE_16S_tax)



tax_summary <- data.frame(
  `TSE_16S`       = FUNS$summarize_taxonomy(TSE_16S),
  `TSE_ITS`       = FUNS$summarize_taxonomy(TSE_ITS),
  `TSE_16S_mock`  = FUNS$summarize_taxonomy(TSE_16S_mock),
  `TSE_ITS_mock`  = FUNS$summarize_taxonomy(TSE_ITS_mock)
)




### Now turns out we had rather low ID rates. So we try again with minBoot = 30

# ITS mock taxonomy
TSE_ITS_mock_tax_lenient <- assignTaxonomy(
  refFasta = "db/sh_general_release_dynamic_s_all_19.02.2025.fasta",
  seqs = rowData(TSE_ITS_mock)$sequence,
  multithread = 16,
  minBoot = 30,
  tryRC = TRUE,
  verbose = TRUE
)
rowData(TSE_ITS_mock_lenient) <- cbind(rowData(TSE_ITS_mock_lenient), TSE_ITS_mock_tax_lenient)
rm(TSE_ITS_mock_tax)

# ITS taxonomy
TSE_ITS_tax_lenient <- assignTaxonomy(
  refFasta = "db/sh_general_release_dynamic_s_all_19.02.2025.fasta",
  seqs = rowData(TSE_ITS)$sequence,
  multithread = 16,
  minBoot = 30,
  tryRC = TRUE,
  verbose = TRUE
)
rowData(TSE_ITS_lenient) <- cbind(rowData(TSE_ITS_lenient), TSE_ITS_tax_lenient)
rm(TSE_ITS_tax_lenient)

# 16S mock taxonomy
TSE_16S_mock_tax_lenient <- assignTaxonomy(
  refFasta = "db/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
  seqs = rowData(TSE_16S_mock)$sequence,
  multithread =16,
  minBoot = 30,
  tryRC = TRUE,
  verbose = TRUE
)
rowData(TSE_16S_mock_lenient) <- cbind(rowData(TSE_16S_mock_lenient), TSE_16S_mock_tax_lenient)
rm(TSE_16S_mock_tax_lenient)

# 16S taxonomy
TSE_16S_tax_lenient <- assignTaxonomy(
  refFasta = "db/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
  seqs = rowData(TSE_16S)$sequence,
  multithread = 16,
  minBoot = 30,
  tryRC = TRUE,
  verbose = TRUE
)
rowData(TSE_16S_lenient) <- cbind(rowData(TSE_16S_lenient), TSE_16S_tax_lenient)
rm(TSE_16S_tax_lenient)



tax_summary_lenient <- data.frame(
  `TSE_16S`       = FUNS$summarize_taxonomy(TSE_16S_lenient),
  `TSE_ITS`       = FUNS$summarize_taxonomy(TSE_ITS_lenient),
  `TSE_16S_mock`  = FUNS$summarize_taxonomy(TSE_16S_mock_lenient),
  `TSE_ITS_mock`  = FUNS$summarize_taxonomy(TSE_ITS_mock_lenient)
)


# We generate some fasta files that contain sequences that did receive a taxonomy assignment 
# at a given level for minBoot = 30, but not for minBoot = 50. 
# Next, we can investigate these sequences, run them through BLASTn, and see whether the 
# taxonomy assignments made with the more lenient parameters make sense. 

set.seed(314) # for consistency
FUNS$sample_taxonomic_gain_ASVs(TSE_ITS, TSE_ITS_lenient, tax_level = "Species", n = 20, outfile = "results/fasta/03_sampled_ITS_species.fasta")
FUNS$sample_taxonomic_gain_ASVs(TSE_ITS, TSE_ITS_lenient, tax_level = "Genus"  , n = 20, outfile = "results/fasta/03_sampled_ITS_genus.fasta")
FUNS$sample_taxonomic_gain_ASVs(TSE_ITS, TSE_ITS_lenient, tax_level = "Family" , n = 20, outfile = "results/fasta/03_sampled_ITS_family.fasta")
FUNS$sample_taxonomic_gain_ASVs(TSE_ITS, TSE_ITS_lenient, tax_level = "Order"  , n = 20, outfile = "results/fasta/03_sampled_ITS_order.fasta")
FUNS$sample_taxonomic_gain_ASVs(TSE_16S, TSE_16S_lenient, tax_level = "Species", n = 20, outfile = "results/fasta/03_sampled_16S_species.fasta")
FUNS$sample_taxonomic_gain_ASVs(TSE_16S, TSE_16S_lenient, tax_level = "Genus"  , n = 20, outfile = "results/fasta/03_sampled_16S_genus.fasta")
FUNS$sample_taxonomic_gain_ASVs(TSE_16S, TSE_16S_lenient, tax_level = "Family" , n = 20, outfile = "results/fasta/03_sampled_16S_family.fasta")
FUNS$sample_taxonomic_gain_ASVs(TSE_16S, TSE_16S_lenient, tax_level = "Order"  , n = 20, outfile = "results/fasta/03_sampled_16S_order.fasta")






Next steps: 
  * Carefully inspect the BLAST outputs
  * Next decide which level of stringency. For ITS I think 30, but for bacteria perhaps 50
  * I can always add a column with a warning for the 30pct. Like a column that says the level of taxonomy down to which it went with 50pct identification
  * Add identifications from other sources? Flag them? 
  * Check other kingdoms in the ITS barcodes, perhaps kick them out
  * Look for AMF
  * remove plastids from bacteria
plastid_hits <- grepl("Chloroplast|Mitochondria", rowData(TSE_16S)$Order, ignore.case = TRUE) |
  grepl("Chloroplast|Mitochondria", rowData(TSE_16S)$Family, ignore.case = TRUE)

TSE_16S_plastids <- TSE_16S[plastid_hits, ]
TSE_16S_filtered <- TSE_16S[!plastid_hits, ]
plastid_counts_per_lib <- colSums(assay(TSE_16S[plastid_hits, ]))
plastid_pct_per_lib <- colSums(assay(TSE_16S[plastid_hits, ])) / colSums(assay(TSE_16S)) * 100

 * Phylogeny: guided or de novo. 

