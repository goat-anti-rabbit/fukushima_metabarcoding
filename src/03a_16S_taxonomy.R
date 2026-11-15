# @ Rik Verdonck 20261028
# Purpose: taxonomy of the 16S ASVs, removal of the ASVs that are not bacteria or archaea, and OTUs.
#
# Taxonomy is assigned with dada2::assignTaxonomy (naive Bayesian classifier) against two reference databases:
#   - Greengenes2 (GG2) 2024.09: the main taxonomy. GG2 follows the genome-based GTDB taxonomy and also names groups
#     that have no Latin name yet (placeholder names such as UBA2421), so it gives a name to the most ASVs and reads.
#   - SILVA 138.2: kept next to it, to look up ASVs with more classical names.
# The two are not combined into one taxonomy: their names differ in many places (e.g. 40 % of the class names).
# For every rank, the bootstrap confidence (0-100) is stored next to the name. 
# Names with a confidence below 50 are NA (the default of assignTaxonomy). 
# The downstream analyses can apply a stricter cutoff with the confidence columns.
# Taxonomy stops at genus: the V3-V4 region (about 420 nt) is often identical between related species. 
# Even the mock strains mostly match several species exactly 
# (e.g. the E. coli ASV 9 Escherichia/Shigella species, B. subtilis 14, Bacillus species), 
# and only 1.3% of the ASVs match any species in SILVA 138.2 exactly (check of 12 Nov 2026).
#
# After taxonomy assignment, we remove the ASVs that are not bacteria or archaea: chloroplasts and mitochondria (according to GG2 or SILVA),
# and ASVs that neither database assigns to Bacteria or Archaea. The removed ASVs are listed in a table.
#
# Finally, the remaining ASVs are clustered into OTUs at 97, 98 and 99 % identity with vsearch (FUNS$cluster_otus).
# The rowData columns otu97, otu98 and otu99 give the OTU of each ASV, named after its centroid (its most abundant
# ASV), e.g. 16S_OTU97_00012 is the 97 % OTU around 16S_ASV00012.
# The counts stay per ASV. An OTU table is made when needed, e.g. mia::agglomerateByVariable(tse, by = "rows", f = "otu97").
# The mocks are clustered on their own.
#
# Input:  results/objects/TSE_16S.rds, TSE_16S_mock.rds (from 02a_16S_denoise.R), reference databases in db/ (see db/README.md)
# Output: results/objects/TSE_16S_taxonomy.rds       TSE of the soil samples with the taxonomy and OTUs in rowData, non-target ASVs removed
#         results/objects/TSE_16S_mock_taxonomy.rds  the same for the mocks (nothing removed)
#         results/objects/TAX16S.rds                 all working objects of this script, incl. the taxonomy of all ASVs
#         results/tables/taxonomy_16S_removed.csv    the removed ASVs, grouped by reason and taxon: number of ASVs and reads
#         results/tables/taxonomy_16S_rates.csv      share of ASVs and reads with a name, per database, rank and confidence cutoff
#         results/figures/16S_taxonomy_rates.pdf     the same as a figure
# In the R session, only TAX16S (plus FUNS and PATH) is left at the end.
# Runtime: about 30 minutes on 20 threads.

source("src/config.R")
source("src/functions/functions.R")
library("dada2")
suppressPackageStartupMessages(library("TreeSummarizedExperiment"))
suppressPackageStartupMessages(library("dplyr"))


### 1. ASVs and reference databases
####################################

TAX16S <- list()

TAX16S$tse      <- readRDS(file.path(PATH$objects, "TSE_16S.rds"))
TAX16S$tse.mock <- readRDS(file.path(PATH$objects, "TSE_16S_mock.rds"))

# Soil and mock ASVs are classified in one go (the classifier is trained only once per database)
TAX16S$sequences <- unique(c(as.character(referenceSeq(TAX16S$tse)), as.character(referenceSeq(TAX16S$tse.mock))))

# DADA2-formatted training files (https://benjjneb.github.io/dada2/training.html), checksums in db/README.md
TAX16S$reference <- c(gg2   = file.path(PATH$db, "gg2_2024_09_toGenus_trainset.fa.gz"),
                      silva = file.path(PATH$db, "silva_nr99_v138.2_toGenus_trainset.fa.gz"))
TAX16S$ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")


### 2. Taxonomy assignment
##########################

# One table per database: a column per rank with the name, and a column per rank with the confidence
# (e.g. family, family_confidence) one row per sequence in TAX16S$sequences
TAX16S$gg2   <- FUNS$assign_taxonomy(TAX16S$sequences, TAX16S$reference["gg2"],   TAX16S$ranks, min_boot = 50, multithread = 20)
TAX16S$silva <- FUNS$assign_taxonomy(TAX16S$sequences, TAX16S$reference["silva"], TAX16S$ranks, min_boot = 50, multithread = 20)


### 3. ASVs that are not bacteria or archaea
############################################

soil  <- match(as.character(referenceSeq(TAX16S$tse)), TAX16S$sequences)   # rows of the soil ASVs in the taxonomy tables
gg2   <- TAX16S$gg2[soil, ]
silva <- TAX16S$silva[soil, ]

# Chloroplasts are a class in GG2 and an order in SILVA, mitochondria are a family in both.
TAX16S$reason <- case_when(
  gg2$class %in% "Chloroplast" | silva$order %in% "Chloroplast"    ~ "chloroplast",
  gg2$family %in% "Mitochondria" | silva$family %in% "Mitochondria" ~ "mitochondria",
  !(gg2$kingdom %in% c("Bacteria", "Archaea") | silva$kingdom %in% c("Bacteria", "Archaea")) ~ "not bacteria or archaea",
  TRUE ~ NA_character_
)
TAX16S$keep <- is.na(TAX16S$reason)

TAX16S$removed <- data.frame(
  reason      = TAX16S$reason,
  gg2_taxon   = paste(gg2$kingdom, gg2$phylum, gg2$class, gg2$order, gg2$family, sep = ";"),
  silva_taxon = paste(silva$kingdom, silva$phylum, silva$class, silva$order, silva$family, sep = ";"),
  reads       = rowSums(assay(TAX16S$tse, "counts"))
) %>%
  filter(!is.na(reason)) %>%
  group_by(reason, gg2_taxon, silva_taxon) %>%
  summarise(asvs = n(), reads = sum(reads), .groups = "drop") %>%
  arrange(reason, desc(reads))
print(as.data.frame(TAX16S$removed))
data.table::fwrite(TAX16S$removed, file.path(PATH$tables, "taxonomy_16S_removed.csv"))
cat("removed:", sum(!TAX16S$keep), "of", length(TAX16S$keep), "ASVs,",
    round(100 * sum(TAX16S$removed$reads) / sum(assay(TAX16S$tse, "counts")), 2), "% of the reads\n")


### 4. Assignment rates of the remaining ASVs
#############################################

reads <- rowSums(assay(TAX16S$tse, "counts"))[TAX16S$keep]
TAX16S$rates <- rbind(
  data.frame(database = "gg2",   FUNS$taxonomy_rates(gg2[TAX16S$keep, ],   reads, tolower(TAX16S$ranks))),
  data.frame(database = "silva", FUNS$taxonomy_rates(silva[TAX16S$keep, ], reads, tolower(TAX16S$ranks)))
)
print(TAX16S$rates[TAX16S$rates$cutoff %in% c(50, 80), ])
data.table::fwrite(TAX16S$rates, file.path(PATH$tables, "taxonomy_16S_rates.csv"))

pdf(file.path(PATH$figures, "16S_taxonomy_rates.pdf"), width = 10, height = 5)
FUNS$plot_taxonomy_rates(TAX16S$rates[TAX16S$rates$database == "gg2", ],   "16S, Greengenes2 2024.09")
FUNS$plot_taxonomy_rates(TAX16S$rates[TAX16S$rates$database == "silva", ], "16S, SILVA 138.2")
dev.off()


### 5. Output: taxonomy in the rowData of the TreeSummarizedExperiment objects
#############################################################################

# Column names as in the 16S object shared in March 2026: GG2 names without prefix (kingdom ... genus) and
# confidences gg2_kingdom_confidence ... gg2_genus_confidence. The SILVA columns all start with silva_.
rename_columns <- function(tax) {
  gg2   <- TAX16S$gg2[match(tax, TAX16S$sequences), ]
  silva <- TAX16S$silva[match(tax, TAX16S$sequences), ]
  names(gg2)   <- sub("^(.*)_confidence$", "gg2_\\1_confidence", names(gg2))
  names(silva) <- paste0("silva_", names(silva))
  S4Vectors::DataFrame(gg2, silva)
}
rowData(TAX16S$tse)      <- cbind(rowData(TAX16S$tse),      rename_columns(as.character(referenceSeq(TAX16S$tse))))
rowData(TAX16S$tse.mock) <- cbind(rowData(TAX16S$tse.mock), rename_columns(as.character(referenceSeq(TAX16S$tse.mock))))
TAX16S$tse      <- FUNS$update_reads(TAX16S$tse[TAX16S$keep, ])   # "reads": the totals after removing the non-target ASVs
TAX16S$tse.mock <- FUNS$update_reads(TAX16S$tse.mock)

# OTUs, clustered after the removal of the non-target ASVs
rowData(TAX16S$tse)      <- cbind(rowData(TAX16S$tse),      S4Vectors::DataFrame(FUNS$cluster_otus(TAX16S$tse)))
rowData(TAX16S$tse.mock) <- cbind(rowData(TAX16S$tse.mock), S4Vectors::DataFrame(FUNS$cluster_otus(TAX16S$tse.mock)))
cat("Number of OTUs (97, 98, 99 %), soil:", sapply(rowData(TAX16S$tse)[c("otu97", "otu98", "otu99")], function(x) length(unique(x))),
    "; mocks:", sapply(rowData(TAX16S$tse.mock)[c("otu97", "otu98", "otu99")], function(x) length(unique(x))), "\n")

# How the taxonomy was made, stored with the objects
TAX16S$info <- list(method = paste("dada2::assignTaxonomy, dada2", packageVersion("dada2")), min_boot = 50,
                    references = TAX16S$reference, md5 = tools::md5sum(TAX16S$reference), date = Sys.Date())
metadata(TAX16S$tse)$taxonomy      <- TAX16S$info
metadata(TAX16S$tse.mock)$taxonomy <- TAX16S$info
metadata(TAX16S$tse)$otus <- metadata(TAX16S$tse.mock)$otus <- "vsearch --cluster_size at 0.97, 0.98 and 0.99, default search settings (FUNS$cluster_otus)"

TAX16S$tse
saveRDS(TAX16S$tse,      file.path(PATH$objects, "TSE_16S_taxonomy.rds"))
saveRDS(TAX16S$tse.mock, file.path(PATH$objects, "TSE_16S_mock_taxonomy.rds"))
saveRDS(TAX16S, file.path(PATH$objects, "TAX16S.rds"))

rm(soil, gg2, silva, reads, rename_columns)
