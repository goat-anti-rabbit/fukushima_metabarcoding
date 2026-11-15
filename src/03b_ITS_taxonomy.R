# @ Rik Verdonck 20261028
# Purpose: taxonomy of the ITS ASVs, removal of the ASVs that are not fungi, and OTUs.
#
# Taxonomy is assigned with dada2::assignTaxonomy (naive Bayesian classifier) against UNITE 10.0 (general FASTA
# release of 19.02.2025, dynamic species hypotheses, all eukaryotes). We use the reference with all eukaryotes, not
# only fungi, so that ASVs of plants, protists and animals are recognised as such instead of being forced into a fungal name.
# For every rank, the bootstrap confidence (0-100) is stored next to the name. 
# Names with a confidence below 50 are NA (the default of assignTaxonomy). 
# The analyses can apply a stricter cutoff with the confidence columns.
#
# UNITE comes in two versions: without and with the "global and 3 % distance singletons" (groups known from a single
# sequence). Both are run, separately, and kept side by side in the rowData (as GG2 and SILVA for 16S):
#   - without singletons: the main taxonomy, down to species
#   - with singletons: down to genus only (columns unite_s_...). Down to species this version has 88,989 different
#     taxonomies, and above 65,536 assignTaxonomy (dada2 1.40.0) gives wrong results because of an integer overflow
#     in its C++ code (see FUNS$unite_to_genus). Cut after the genus it has 19,808. (The version without singletons
#     has 64,976 down to species and is fine.) FUNS$assign_taxonomy checks this number for every reference.
# Section 4 compares the two versions on an equal footing, both cut after the genus.
#
# We then remove the ASVs that are not fungi according to the main taxonomy (other eukaryotes, or kingdom unknown).
# They are described in a table.
#
# Finally, the fungal ASVs are clustered into OTUs at 97, 98 and 99 % identity with vsearch (FUNS$cluster_otus).
# The rowData columns otu97, otu98 and otu99 give the OTU of each ASV, named after its centroid (its most abundant
# ASV), e.g. ITS_OTU97_00012 is the 97 % OTU around ITS_ASV00012.
# The counts stay per ASV. An OTU table is made when needed, e.g. mia::agglomerateByVariable(tse, by = "rows", f = "otu97").
# The mocks are clustered on their own.
#
# Input:  results/objects/TSE_ITS.rds, TSE_ITS_mock.rds (from 02b_ITS_denoise.R), reference databases in db/ (see db/README.md)
# Output: results/objects/TSE_ITS_taxonomy.rds        TSE of the soil samples with the taxonomy and OTUs in rowData, non-fungal ASVs removed
#         results/objects/TSE_ITS_mock_taxonomy.rds   the same for the mocks (nothing removed)
#         results/objects/TAXITS.rds                  all working objects of this script, incl. the taxonomy of all ASVs
#         results/tables/taxonomy_ITS_removed.csv     the removed ASVs, grouped by kingdom and phylum: number of ASVs and reads
#         results/tables/taxonomy_ITS_rates.csv       share of ASVs and reads with a name, per UNITE version, rank and confidence cutoff
#         results/tables/taxonomy_ITS_singletons_kingdoms.csv, taxonomy_ITS_singletons_rates.csv,
#         taxonomy_ITS_singletons_agreement.csv       comparison of UNITE without and with singletons, both cut after the genus
#         results/figures/ITS_taxonomy_rates.pdf      the rates as a figure
# In the R session, only TAXITS (plus FUNS and PATH) is left at the end.
# Runtime: about 70 minutes on 20 threads.

source("src/config.R")
source("src/functions/functions.R")
library("dada2")
suppressPackageStartupMessages(library("TreeSummarizedExperiment"))
suppressPackageStartupMessages(library("dplyr"))


### 1. ASVs and reference databases
####################################

TAXITS <- list()

TAXITS$tse      <- readRDS(file.path(PATH$objects, "TSE_ITS.rds"))
TAXITS$tse.mock <- readRDS(file.path(PATH$objects, "TSE_ITS_mock.rds"))

# Soil and mock ASVs are classified in one go (the classifier is trained only once per database)
TAXITS$sequences <- unique(c(as.character(referenceSeq(TAXITS$tse)), as.character(referenceSeq(TAXITS$tse.mock))))

# UNITE general FASTA releases (https://unite.ut.ee/repository.php), checksums in db/README.md
TAXITS$reference <- c(unite   = file.path(PATH$db, "sh_general_release_dynamic_all_19.02.2025.fasta"),
                      unite_s = file.path(PATH$db, "sh_general_release_dynamic_s_all_19.02.2025.fasta"))
TAXITS$ranks       <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAXITS$genus_ranks <- TAXITS$ranks[1:6]


### 2. Taxonomy assignment
##########################

# Each run gives a table with a column per rank with the name, and a column per rank with the confidence
# (e.g. family, family_confidence) one row per sequence in TAXITS$sequences.

# Main taxonomy: UNITE without singletons, down to species
TAXITS$unite <- FUNS$assign_taxonomy(TAXITS$sequences, TAXITS$reference["unite"], TAXITS$ranks, min_boot = 50, multithread = 20)

# Both versions cut after the genus (temporary copies of the reference files): the one with singletons is kept
# next to the main taxonomy, the one without singletons is only used for the comparison in section 4
for (v in c("unite", "unite_s")) {
  ref_genus <- tempfile(fileext = ".fasta")
  FUNS$unite_to_genus(TAXITS$reference[v], ref_genus)
  TAXITS[[paste0(v, "_genus")]] <- FUNS$assign_taxonomy(TAXITS$sequences, ref_genus, TAXITS$genus_ranks, min_boot = 50, multithread = 20)
  unlink(ref_genus)
}


### 3. ASVs that are not fungi
##############################

soil  <- match(as.character(referenceSeq(TAXITS$tse)), TAXITS$sequences)   # rows of the soil ASVs in the taxonomy tables
unite <- TAXITS$unite[soil, ]
reads <- rowSums(assay(TAXITS$tse, "counts"))

TAXITS$keep <- unite$kingdom %in% "Fungi"

TAXITS$removed <- data.frame(
  kingdom = ifelse(is.na(unite$kingdom), "unknown", unite$kingdom),
  phylum  = unite$phylum,
  keep    = TAXITS$keep,
  reads   = reads
) %>%
  filter(!keep) %>%
  group_by(kingdom, phylum) %>%
  summarise(asvs = n(), reads = sum(reads), .groups = "drop") %>%
  arrange(kingdom, desc(reads))
print(as.data.frame(TAXITS$removed))
data.table::fwrite(TAXITS$removed, file.path(PATH$tables, "taxonomy_ITS_removed.csv"))
cat("removed:", sum(!TAXITS$keep), "of", length(TAXITS$keep), "ASVs,",
    round(100 * sum(TAXITS$removed$reads) / sum(reads), 2), "% of the reads\n")


### 4. UNITE without and with singletons, both cut after the genus
##################################################################

a <- TAXITS$unite_genus[soil, ]
b <- TAXITS$unite_s_genus[soil, ]

# Kingdom of all soil ASVs according to both versions
TAXITS$singletons_kingdoms <- data.frame(
  without_singletons = ifelse(is.na(a$kingdom), "unknown", a$kingdom),
  with_singletons    = ifelse(is.na(b$kingdom), "unknown", b$kingdom),
  reads              = reads
) %>%
  group_by(without_singletons, with_singletons) %>%
  summarise(asvs = n(), reads = sum(reads), .groups = "drop") %>%
  arrange(desc(reads))
print(as.data.frame(TAXITS$singletons_kingdoms))
data.table::fwrite(TAXITS$singletons_kingdoms, file.path(PATH$tables, "taxonomy_ITS_singletons_kingdoms.csv"))

# For the ASVs that both versions call fungi: share with a name at each rank (confidence >= 50 and >= 80), and how
# often the two versions give the same name where both give a real name with at least that confidence
fungi <- a$kingdom %in% "Fungi" & b$kingdom %in% "Fungi"
TAXITS$singletons_rates <- rbind(
  data.frame(version = "without singletons", FUNS$taxonomy_rates(a[fungi, ], reads[fungi], tolower(TAXITS$genus_ranks), cutoffs = c(50, 80))),
  data.frame(version = "with singletons",    FUNS$taxonomy_rates(b[fungi, ], reads[fungi], tolower(TAXITS$genus_ranks), cutoffs = c(50, 80)))
)
TAXITS$singletons_agreement <- FUNS$name_agreement(a[fungi, ], b[fungi, ], tolower(TAXITS$genus_ranks), cutoffs = c(50, 80))
print(TAXITS$singletons_rates)
print(TAXITS$singletons_agreement)
data.table::fwrite(TAXITS$singletons_rates,     file.path(PATH$tables, "taxonomy_ITS_singletons_rates.csv"))
data.table::fwrite(TAXITS$singletons_agreement, file.path(PATH$tables, "taxonomy_ITS_singletons_agreement.csv"))


### 5. Assignment rates of the remaining ASVs
#############################################

TAXITS$rates <- rbind(
  data.frame(database = "unite",   FUNS$taxonomy_rates(unite[TAXITS$keep, ],                reads[TAXITS$keep], tolower(TAXITS$ranks))),
  data.frame(database = "unite_s", FUNS$taxonomy_rates(TAXITS$unite_s_genus[soil, ][TAXITS$keep, ], reads[TAXITS$keep], tolower(TAXITS$genus_ranks)))
)
print(TAXITS$rates[TAXITS$rates$cutoff %in% c(50, 80), ])
data.table::fwrite(TAXITS$rates, file.path(PATH$tables, "taxonomy_ITS_rates.csv"))

pdf(file.path(PATH$figures, "ITS_taxonomy_rates.pdf"), width = 10, height = 5)
FUNS$plot_taxonomy_rates(TAXITS$rates[TAXITS$rates$database == "unite", ],   "ITS, UNITE 10.0 without singletons")
FUNS$plot_taxonomy_rates(TAXITS$rates[TAXITS$rates$database == "unite_s", ], "ITS, UNITE 10.0 with singletons (to genus)")
dev.off()


### 6. Output: taxonomy in the rowData of the TreeSummarizedExperiment objects
#############################################################################

# Main taxonomy: UNITE names without prefix (kingdom ... species) and confidences unite_kingdom_confidence ...
# unite_species_confidence (as for 16S, where the confidences are gg2_..._confidence).
# UNITE with singletons: all columns start with unite_s_ (as silva_ for 16S).
add_taxonomy <- function(tse) {
  rows <- match(as.character(referenceSeq(tse)), TAXITS$sequences)
  u   <- TAXITS$unite[rows, ]
  u_s <- TAXITS$unite_s_genus[rows, ]
  names(u)   <- sub("^(.*)_confidence$", "unite_\\1_confidence", names(u))
  names(u_s) <- paste0("unite_s_", names(u_s))
  rowData(tse) <- cbind(rowData(tse), S4Vectors::DataFrame(u, u_s))
  metadata(tse)$taxonomy <- list(method = paste("dada2::assignTaxonomy, dada2", packageVersion("dada2")), min_boot = 50,
                                 references = TAXITS$reference, md5 = tools::md5sum(TAXITS$reference),
                                 note = "unite_s_ columns: reference with singletons, taxonomy cut after the genus",
                                 date = Sys.Date())
  tse
}
TAXITS$tse      <- FUNS$update_reads(add_taxonomy(TAXITS$tse)[TAXITS$keep, ])   # "reads": the totals after removing the non-fungal ASVs
TAXITS$tse.mock <- FUNS$update_reads(add_taxonomy(TAXITS$tse.mock))

# OTUs, clustered after the removal of the non-fungal ASVs
rowData(TAXITS$tse)      <- cbind(rowData(TAXITS$tse),      S4Vectors::DataFrame(FUNS$cluster_otus(TAXITS$tse)))
rowData(TAXITS$tse.mock) <- cbind(rowData(TAXITS$tse.mock), S4Vectors::DataFrame(FUNS$cluster_otus(TAXITS$tse.mock)))
metadata(TAXITS$tse)$otus <- metadata(TAXITS$tse.mock)$otus <- "vsearch --cluster_size at 0.97, 0.98 and 0.99, default search settings (FUNS$cluster_otus)"
cat("Number of OTUs (97, 98, 99 %), soil:", sapply(rowData(TAXITS$tse)[c("otu97", "otu98", "otu99")], function(x) length(unique(x))),
    "; mocks:", sapply(rowData(TAXITS$tse.mock)[c("otu97", "otu98", "otu99")], function(x) length(unique(x))), "\n")

TAXITS$tse
saveRDS(TAXITS$tse,      file.path(PATH$objects, "TSE_ITS_taxonomy.rds"))
saveRDS(TAXITS$tse.mock, file.path(PATH$objects, "TSE_ITS_mock_taxonomy.rds"))
saveRDS(TAXITS, file.path(PATH$objects, "TAXITS.rds"))

rm(soil, unite, reads, v, ref_genus, a, b, fungi, add_taxonomy)
