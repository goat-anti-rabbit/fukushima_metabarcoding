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
library("DECIPHER")
library("ShortRead")
library("phangorn")
library("Biostrings")
library("ggplot2")
library("patchwork")
library("tidyverse")


### 1. Reading in data, setting up TreeSE objects
########################################################

# Some small custom functions:
source("./src/functions/functions.R")

# Data:
#load("results/objects/02_TreeSE_objects.RData")
TSE_16S <- readRDS("results/objects/TSE_16S.rds")

# Read in trees:
tree_gg_path    <- "exported-tree-gg/tree.nwk"        # pruned or full GG tree
tree_silva_path <- "exported-tree-silva128/tree.nwk"  # pruned or full SILVA tree
assay_name <- "counts"   # change if your count assay has another name
group_var  <- "site"     # metadata column to color by
shape_var  <- NULL       # e.g. "contamination_group" or NULL
title_var  <- "sample_id" # only used if you want labels later; not used now

prune_tree_to_tse <- function(tree, tse) {
  feat_ids <- rownames(tse)
  
  if (is.null(feat_ids)) {
    stop("rownames(tse) are NULL. They must contain the ASV IDs.")
  }
  
  missing_in_tree <- setdiff(feat_ids, tree$tip.label)
  if (length(missing_in_tree) > 0) {
    stop(
      sprintf(
        "Tree is missing %d features from the TSE. Example: %s",
        length(missing_in_tree),
        paste(head(missing_in_tree, 10), collapse = ", ")
      )
    )
  }
  
  tree2 <- keep.tip(tree, feat_ids)
  tree2 <- collapse.singles(tree2)
  
  # Reorder tips to match rownames(tse)
  tree2 <- keep.tip(tree2, feat_ids)
  
  if (!setequal(tree2$tip.label, feat_ids)) {
    stop("After pruning, tree tip labels do not match TSE rownames.")
  }
  
  tree2
}

rownames(TSE_16S) <- paste("seq",1:19055,sep="_")
colnames(TSE_16S) <- colData(TSE_16S)$DNA_extract
TSE_16S <- TSE_16S[, !colnames(TSE_16S) %in% c("mock1", "mock2")]

# Caution: this removes rare features!!!
# TSE_16S <- rarefyAssay(TSE_16S , assay.type = "counts", name = "rarefied",sample=20000, replace=F)
TSE_16S <- transformAssay(TSE_16S , method = "relabundance") 
TSE_16S <- transformAssay(TSE_16S , assay.type = "counts", method = "clr", pseudocount = TRUE)

tree_gg_full    <- read.tree(tree_gg_path)
tree_silva_full <- read.tree(tree_silva_path)
tree_gg    <- prune_tree_to_tse(tree_gg_full, TSE_16S)
tree_silva <- prune_tree_to_tse(tree_silva_full, TSE_16S)

cat("GG tips after pruning:", Ntip(tree_gg), "\n")
cat("SILVA tips after pruning:", Ntip(tree_silva), "\n")
cat("GG rooted:", is.rooted(tree_gg), "\n")
cat("SILVA rooted:", is.rooted(tree_silva), "\n")
cat("GG negative branches:", any(tree_gg$edge.length < 0), "\n")
cat("SILVA negative branches:", any(tree_silva$edge.length < 0), "\n")


TSE_16S_gg    <- TSE_16S
TSE_16S_silva <- TSE_16S


rowTree(TSE_16S_gg)    <- tree_gg
rowTree(TSE_16S_silva) <- tree_silva

identical(sort(rowTree(TSE_16S_gg)$tip.label), sort(rownames(TSE_16S_gg)))
identical(sort(rowTree(TSE_16S_silva)$tip.label), sort(rownames(TSE_16S_silva)))


# Unifrac for silva:
TSE_16S_silva <- addDissimilarity(
  x = TSE_16S_silva,
  assay.type = "counts", # Whether we use abundance or counts, it doesn't matter because the unifrac function uses normalized abundances anyway (in case of weighted)
  method = "unifrac",
  weighted = TRUE,
  name = "w_unifrac"
)

TSE_16S_silva <- addDissimilarity(
  x = TSE_16S_silva,
  assay.type = "counts",
  method = "unifrac",
  weighted = F,
  name = "u_unifrac_counts"
)

TSE_16S_silva <- addDissimilarity(
  x = TSE_16S_silva,
  assay.type = "counts",
  method = "unifrac",
  weighted = FALSE,
  niter = 100,      # Averages distances over 100 subsampling iterations
  sample = 10000,   # Set your target depth here
  name = "u_unifrac_rarefied"
)

# Unifrac for gg:
TSE_16S_gg<- addDissimilarity(
  x = TSE_16S_gg,
  assay.type = "counts", 
  method = "unifrac",
  weighted = TRUE,
  name = "w_unifrac"
)

TSE_16S_gg <- addDissimilarity(
  x = TSE_16S_gg,
  assay.type = "counts",
  method = "unifrac",
  weighted = F,
  name = "u_unifrac_counts"
)

TSE_16S_gg <- addDissimilarity(
  x = TSE_16S_gg,
  assay.type = "counts",
  method = "unifrac",
  weighted = FALSE,
  niter = 100,     
  sample = 10000,
  name = "u_unifrac_rarefied"
)



pcoa_silva <- cmdscale(as.dist(metadata(TSE_16S_silva)$u_unifrac_counts), k = 2, eig = TRUE)
pcoa_gg    <- cmdscale(as.dist(metadata(TSE_16S_gg)$u_unifrac_counts), k = 2, eig = TRUE)

df_silva <- data.frame(Axis1 = pcoa_silva$points[, 1], Axis2 = pcoa_silva$points[, 2], sampling_site = colData(TSE_16S_silva)$sampling_site,sample_id = colnames(TSE_16S_silva))
df_gg    <- data.frame(Axis1 = pcoa_gg$points[, 1], Axis2 = pcoa_gg$points[, 2], sampling_site = colData(TSE_16S_gg)$sampling_site, sample_id = colnames(TSE_16S_gg))

var_silva <- pcoa_silva$eig / sum(pcoa_silva$eig[pcoa_silva$eig > 0])
var_gg    <- pcoa_gg$eig / sum(pcoa_gg$eig[pcoa_gg$eig > 0])

p_silva <- ggplot(df_silva, aes(Axis1, Axis2, color = sampling_site)) +
  geom_point(size = 2.5) +
  labs(title = "Unweighted UniFrac - SEPP SILVA128 (counts)", x = paste0("PCoA1 (", round(100 * var_silva[1], 1), "%)"), y = paste0("PCoA2 (", round(100 * var_silva[2], 1), "%)")) +
  theme_bw()

p_gg   <- ggplot(df_gg, aes(Axis1, Axis2, color = sampling_site)) +
  geom_point(size = 2.5) +
  labs(
    title = "Unweighted UniFrac - SEPP GG13_8 (counts)", x = paste0("PCoA1 (", round(100 * var_gg[1], 1), "%)"), y = paste0("PCoA2 (", round(100 * var_gg[2], 1), "%)"))+
  theme_bw()

p_silva + p_gg



pcoa_silva <- cmdscale(as.dist(metadata(TSE_16S_silva)$u_unifrac_rarefied), k = 2, eig = TRUE)
pcoa_gg    <- cmdscale(as.dist(metadata(TSE_16S_gg)$u_unifrac_rarefied), k = 2, eig = TRUE)

df_silva <- data.frame(Axis1 = pcoa_silva$points[, 1], Axis2 = pcoa_silva$points[, 2], sampling_site = colData(TSE_16S_silva)$sampling_site,sample_id = colnames(TSE_16S_silva))
df_gg    <- data.frame(Axis1 = pcoa_gg$points[, 1], Axis2 = pcoa_gg$points[, 2], sampling_site = colData(TSE_16S_gg)$sampling_site, sample_id = colnames(TSE_16S_gg))

var_silva <- pcoa_silva$eig / sum(pcoa_silva$eig[pcoa_silva$eig > 0])
var_gg    <- pcoa_gg$eig / sum(pcoa_gg$eig[pcoa_gg$eig > 0])

p_silva <- ggplot(df_silva, aes(Axis1, Axis2, color = sampling_site)) +
  geom_point(size = 2.5) +
  labs(title = "Unweighted UniFrac - SEPP SILVA128 (rarified)", x = paste0("PCoA1 (", round(100 * var_silva[1], 1), "%)"), y = paste0("PCoA2 (", round(100 * var_silva[2], 1), "%)")) +
  theme_bw()

p_gg   <- ggplot(df_gg, aes(Axis1, Axis2, color = sampling_site)) +
  geom_point(size = 2.5) +
  labs(
    title = "Unweighted UniFrac - SEPP GG13_8 (rarified)", x = paste0("PCoA1 (", round(100 * var_gg[1], 1), "%)"), y = paste0("PCoA2 (", round(100 * var_gg[2], 1), "%)"))+
  theme_bw()

p_silva + p_gg





pcoa_silva <- cmdscale(as.dist(metadata(TSE_16S_silva)$w_unifrac), k = 2, eig = TRUE)
pcoa_gg    <- cmdscale(as.dist(metadata(TSE_16S_gg)$w_unifrac), k = 2, eig = TRUE)

df_silva <- data.frame(Axis1 = pcoa_silva$points[, 1], Axis2 = pcoa_silva$points[, 2], sampling_site = colData(TSE_16S_silva)$sampling_site,sample_id = colnames(TSE_16S_silva))
df_gg    <- data.frame(Axis1 = pcoa_gg$points[, 1], Axis2 = pcoa_gg$points[, 2], sampling_site = colData(TSE_16S_gg)$sampling_site, sample_id = colnames(TSE_16S_gg))

var_silva <- pcoa_silva$eig / sum(pcoa_silva$eig[pcoa_silva$eig > 0])
var_gg    <- pcoa_gg$eig / sum(pcoa_gg$eig[pcoa_gg$eig > 0])

p_silva <- ggplot(df_silva, aes(Axis1, Axis2, color = sampling_site)) +
  geom_point(size = 2.5) +
  labs(title = "Weighted UniFrac - SEPP SILVA128", x = paste0("PCoA1 (", round(100 * var_silva[1], 1), "%)"), y = paste0("PCoA2 (", round(100 * var_silva[2], 1), "%)")) +
  theme_bw()

p_gg   <- ggplot(df_gg, aes(Axis1, Axis2, color = sampling_site)) +
  geom_point(size = 2.5) +
  labs(
    title = "Weighted UniFrac - SEPP GG13_8", x = paste0("PCoA1 (", round(100 * var_gg[1], 1), "%)"), y = paste0("PCoA2 (", round(100 * var_gg[2], 1), "%)"))+
  theme_bw()

p_silva + p_gg






































# taxonomydown to genus level:
# taxonomy down to genus level, greengenes2:
TSE_16S_tax_gg2 <- assignTaxonomy(
  refFasta = "db/gg2_2024_09_toGenus_trainset.fa.gz",
  seqs = rowData(TSE_16S)$sequence,
  taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  multithread = 16,
  tryRC = TRUE,
  verbose = TRUE,
  minBoot = 30,
  outputBootstraps=T
)

TSE_16S_tax_gg2$tax[,1] <- gsub("^.__", "", TSE_16S_tax_gg2$tax[,1])
TSE_16S_tax_gg2$tax[,2] <- gsub("^.__", "", TSE_16S_tax_gg2$tax[,2])
TSE_16S_tax_gg2$tax[,3] <- gsub("^.__", "", TSE_16S_tax_gg2$tax[,3])
TSE_16S_tax_gg2$tax[,4] <- gsub("^.__", "", TSE_16S_tax_gg2$tax[,4])
TSE_16S_tax_gg2$tax[,5] <- gsub("^.__", "", TSE_16S_tax_gg2$tax[,5])
TSE_16S_tax_gg2$tax[,6] <- gsub("^.__", "", TSE_16S_tax_gg2$tax[,6])


# taxonomy down to genus level, silva
TSE_16S_tax_silva <- assignTaxonomy(
  refFasta = "db/silva_nr99_v138.2_toGenus_trainset.fa.gz",
  seqs = rowData(TSE_16S)$sequence,
  taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  multithread = 16,
  tryRC = TRUE,
  verbose = TRUE,
  minBoot = 30,
  outputBootstraps=T
)

save.image("session20260318_taxonomy.RData")

# taxonomy down to genus level, gtdb
TSE_16S_tax_gtdb <- assignTaxonomy(
  refFasta = "db/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz",
  seqs = rowData(TSE_16S)$sequence,
  taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  multithread = 16,
  tryRC = TRUE,
  verbose = TRUE,
  minBoot = 30,
  outputBootstraps=T
)

save.image("session20260318_taxonomy.RData")

# taxonomy down to genus level, rdp
TSE_16S_tax_rdp <- assignTaxonomy(
  refFasta = "db/rdp_19_toGenus_trainset.fa.gz",
  seqs = rowData(TSE_16S)$sequence,
  taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  multithread = 16,
  tryRC = TRUE,
  verbose = TRUE,
  minBoot = 30,
  outputBootstraps=T
)

save.image("session20260318_taxonomy.RData")


# taxonomy down to genus level, refseq
TSE_16S_tax_refseq <- assignTaxonomy(
  refFasta = "db/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz",
  seqs = rowData(TSE_16S)$sequence,
  taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  multithread = 16,
  tryRC = TRUE,
  verbose = TRUE,
  minBoot = 30,
  outputBootstraps=T
)

save.image("session20260318_taxonomy.RData")

dna <- DNAStringSet(rowData(TSE_16S)$sequence) # Create a DNAStringSet from the ASVs

load("db/SILVA_SSU_r138.2_v2.RData")
ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=FALSE, threshold=30)
ranks <- c("domain", "phylum", "class", "order", "family", "genus") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid_silva <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
taxConf_silva <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$confidence[m]
  taxa
}))
save.image("session20260318_taxonomy.RData")


load("db/RDP_TrainingSet_v18-mod.RData")
ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=FALSE, threshold=30, )
ranks <- c("domain", "phylum", "class", "order", "family", "genus") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid_rdp <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
taxConf_rdp <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$confidence[m]
  taxa
}))
save.image("session20260318_taxonomy.RData")


load("db/GTDB_r226_classifier.RData")
ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=FALSE, threshold=30, )
ranks <- c("domain", "phylum", "class", "order", "family", "genus") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid_gtdb <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
taxConf_gtdb <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$confidence[m]
  taxa
}))
save.image("session20260318_taxonomy.RData")


# For Kingdom:
kingdom.tax <- cbind.data.frame(
  "dada_gg2"   = unname(TSE_16S_tax_gg2$tax[,"Kingdom"]),
  "dada_silva" = unname(TSE_16S_tax_silva$tax[,"Kingdom"]),
  "dada_gtdb"  = unname(TSE_16S_tax_gtdb$tax[,"Kingdom"]),
  "dada_rdp"   = unname(TSE_16S_tax_rdp$tax[,"Kingdom"]),
  "dada_refseq"= unname(TSE_16S_tax_refseq$tax[,"Kingdom"]),
  "dec_silva"  = unname(taxid_silva[,1]) ,
  "dec_rdp"    = unname(taxid_rdp[,1]),
  "dec_gtdb"    = unname(taxid_gtdb[,1]))
kingdom.conf <- cbind.data.frame(
  "dada_gg2"   = unname(TSE_16S_tax_gg2$boot[,"Kingdom"]),
  "dada_silva" = unname(TSE_16S_tax_silva$boot[,"Kingdom"]),
  "dada_gtdb"  = unname(TSE_16S_tax_gtdb$boot[,"Kingdom"]),
  "dada_rdp"   = unname(TSE_16S_tax_rdp$boot[,"Kingdom"]),
  "dada_refseq"= unname(TSE_16S_tax_refseq$boot[,"Kingdom"]),
  "dec_silva"  = unname(taxConf_silva[,1]) ,
  "dec_rdp"    = unname(taxConf_rdp[,1]),
  "dec_gtdb"    = unname(taxConf_gtdb[,1]))

kingdom.real <- sapply(X=kingdom.tax, function(x){c(is.na(x) + grepl("Incertae", x) + c(nchar(x,keepNA=F) == 0))==0})



# For phylum:
phylum.tax <- cbind.data.frame(
  "dada_gg2"   = unname(TSE_16S_tax_gg2$tax[,"Phylum"]),
  "dada_silva" = unname(TSE_16S_tax_silva$tax[,"Phylum"]),
  "dada_gtdb"  = unname(TSE_16S_tax_gtdb$tax[,"Phylum"]),
  "dada_rdp"   = unname(TSE_16S_tax_rdp$tax[,"Phylum"]),
  "dada_refseq"= unname(TSE_16S_tax_refseq$tax[,"Phylum"]),
  "dec_silva"  = unname(taxid_silva[,2]) ,
  "dec_rdp"    = unname(taxid_rdp[,2]),
  "dec_gtdb"    = unname(taxid_gtdb[,2]))
phylum.conf <- cbind.data.frame(
  "dada_gg2"   = unname(TSE_16S_tax_gg2$boot[,"Phylum"]),  
  "dada_silva" = unname(TSE_16S_tax_silva$boot[,"Phylum"]),
  "dada_gtdb"  = unname(TSE_16S_tax_gtdb$boot[,"Phylum"]),
  "dada_rdp"   = unname(TSE_16S_tax_rdp$boot[,"Phylum"]),
  "dada_refseq"= unname(TSE_16S_tax_refseq$boot[,"Phylum"]),
  "dec_silva"  = unname(taxConf_silva[,2]) ,
  "dec_rdp"    = unname(taxConf_rdp[,2]),
  "dec_gtdb"    = unname(taxConf_gtdb[,2]))

phylum.real <- sapply(X=phylum.tax, function(x){c(is.na(x) + grepl("Incertae", x) + c(nchar(x,keepNA=F) == 0))==0})
#phylum.conf[!phylum.real] <- 0

# Now for class:
class.tax <- cbind.data.frame(
  "dada_gg2"   = unname(TSE_16S_tax_gg2$tax[,"Class"]),  
  "dada_silva" = unname(TSE_16S_tax_silva$tax[,"Class"]),
  "dada_gtdb"  = unname(TSE_16S_tax_gtdb$tax[,"Class"]),
  "dada_rdp"   = unname(TSE_16S_tax_rdp$tax[,"Class"]),
  "dada_refseq"= unname(TSE_16S_tax_refseq$tax[,"Class"]),
  "dec_silva"  = unname(taxid_silva[,3]) ,
  "dec_rdp"    = unname(taxid_rdp[,3]),
  "dec_gtdb"    = unname(taxid_gtdb[,3]))
class.conf <- cbind.data.frame(
  "dada_gg2"   = unname(TSE_16S_tax_gg2$boot[,"Class"]),
  "dada_silva" = unname(TSE_16S_tax_silva$boot[,"Class"]),
  "dada_gtdb"  = unname(TSE_16S_tax_gtdb$boot[,"Class"]),
  "dada_rdp"   = unname(TSE_16S_tax_rdp$boot[,"Class"]),
  "dada_refseq"= unname(TSE_16S_tax_refseq$boot[,"Class"]),
  "dec_silva"  = unname(taxConf_silva[,3]) ,
  "dec_rdp"    = unname(taxConf_rdp[,3]),
  "dec_gtdb"    = unname(taxConf_gtdb[,3]))
class.real <- sapply(X=class.tax, function(x){c(is.na(x) + grepl("Incertae", x) + c(nchar(x,keepNA=F) == 0))==0})
#class.conf[!class.real] <- 0


# Now for order:
order.tax <- cbind.data.frame(
  "dada_gg2"   = unname(TSE_16S_tax_gg2$tax[,"Order"]),
  "dada_silva" = unname(TSE_16S_tax_silva$tax[,"Order"]),
  "dada_gtdb"  = unname(TSE_16S_tax_gtdb$tax[,"Order"]),
  "dada_rdp"   = unname(TSE_16S_tax_rdp$tax[,"Order"]),
  "dada_refseq"= unname(TSE_16S_tax_refseq$tax[,"Order"]),
  "dec_silva"  = unname(taxid_silva[,4]) ,
  "dec_rdp"    = unname(taxid_rdp[,4]),
  "dec_gtdb"    = unname(taxid_gtdb[,4]))
order.conf <- cbind.data.frame(
  "dada_gg2"   = unname(TSE_16S_tax_gg2$boot[,"Order"]),
  "dada_silva" = unname(TSE_16S_tax_silva$boot[,"Order"]),
  "dada_gtdb"  = unname(TSE_16S_tax_gtdb$boot[,"Order"]),
  "dada_rdp"   = unname(TSE_16S_tax_rdp$boot[,"Order"]),
  "dada_refseq"= unname(TSE_16S_tax_refseq$boot[,"Order"]),
  "dec_silva"  = unname(taxConf_silva[,4]) ,
  "dec_rdp"    = unname(taxConf_rdp[,4]),
  "dec_gtdb"    = unname(taxConf_gtdb[,4]))
order.real <- sapply(X=order.tax, function(x){c(is.na(x) + grepl("Incertae", x) + c(nchar(x,keepNA=F) == 0))==0})
#order.conf[!order.real] <- 0



# For family:
family.tax <- cbind.data.frame(
  "dada_gg2"   = unname(TSE_16S_tax_gg2$tax[,"Family"]),
  "dada_silva" = unname(TSE_16S_tax_silva$tax[,"Family"]),
  "dada_gtdb"  = unname(TSE_16S_tax_gtdb$tax[,"Family"]),
  "dada_rdp"   = unname(TSE_16S_tax_rdp$tax[,"Family"]),
  "dada_refseq"= unname(TSE_16S_tax_refseq$tax[,"Family"]),
  "dec_silva"  = unname(taxid_silva[,5]) ,
  "dec_rdp"    = unname(taxid_rdp[,5]),
  "dec_gtdb"    = unname(taxid_gtdb[,5]))
family.conf <- cbind.data.frame(
  "dada_gg2"   = unname(TSE_16S_tax_gg2$boot[,"Family"]),
  "dada_silva" = unname(TSE_16S_tax_silva$boot[,"Family"]),
  "dada_gtdb"  = unname(TSE_16S_tax_gtdb$boot[,"Family"]),
  "dada_rdp"   = unname(TSE_16S_tax_rdp$boot[,"Family"]),
  "dada_refseq"= unname(TSE_16S_tax_refseq$boot[,"Family"]),
  "dec_silva"  = unname(taxConf_silva[,5]) ,
  "dec_rdp"    = unname(taxConf_rdp[,5]),
  "dec_gtdb"    = unname(taxConf_gtdb[,5]))
family.real <- sapply(X=family.tax, function(x){c(is.na(x) + grepl("Incertae", x) + c(nchar(x,keepNA=F) == 0))==0})
#family.conf[!family.real] <- 0

# And finally for genus:

genus.tax <- cbind.data.frame(
  "dada_gg2"   = unname(TSE_16S_tax_gg2$tax[,"Genus"]),
  "dada_silva" = unname(TSE_16S_tax_silva$tax[,"Genus"]),
  "dada_gtdb"  = unname(TSE_16S_tax_gtdb$tax[,"Genus"]),
  "dada_rdp"   = unname(TSE_16S_tax_rdp$tax[,"Genus"]),
  "dada_refseq"= unname(TSE_16S_tax_refseq$tax[,"Genus"]),
  "dec_silva"  = unname(taxid_silva[,6]) ,
  "dec_rdp"    = unname(taxid_rdp[,6]),
  "dec_gtdb"    = unname(taxid_gtdb[,6]))
genus.conf <- cbind.data.frame(
  "dada_gg2"   = unname(TSE_16S_tax_gg2$boot[,"Genus"]),
  "dada_silva" = unname(TSE_16S_tax_silva$boot[,"Genus"]),
  "dada_gtdb"  = unname(TSE_16S_tax_gtdb$boot[,"Genus"]),
  "dada_rdp"   = unname(TSE_16S_tax_rdp$boot[,"Genus"]),
  "dada_refseq"= unname(TSE_16S_tax_refseq$boot[,"Genus"]),
  "dec_silva"  = unname(taxConf_silva[,6]) ,
  "dec_rdp"    = unname(taxConf_rdp[,6]),
  "dec_gtdb"    = unname(taxConf_gtdb[,6]))
genus.real <- sapply(X=genus.tax, function(x){c(is.na(x) + grepl("Incertae", x) + c(nchar(x,keepNA=F) == 0))==0})
#genus.conf[!genus.real] <- 0



tax_long <- kingdom.tax %>%
  mutate(ASV = row_number()) %>%
  mutate(sequence = as.character(dna)) %>% 
  pivot_longer(-c(ASV,sequence), names_to = "method", values_to = "taxon") %>% 
  mutate(level="kingdom") 

conf_long <- kingdom.conf %>%
  mutate(ASV = row_number()) %>%
  pivot_longer(-ASV, names_to = "method", values_to = "confidence")

df.kingdom <- dplyr::left_join(tax_long, conf_long, by = c("ASV", "method"))%>%
  relocate(sequence, .after=confidence)


tax_long <- phylum.tax %>%
  mutate(ASV = row_number()) %>%
  mutate(sequence = as.character(dna)) %>% 
  pivot_longer(-c(ASV,sequence), names_to = "method", values_to = "taxon") %>% 
  mutate(level="phylum") 

conf_long <- phylum.conf %>%
  mutate(ASV = row_number()) %>%
  pivot_longer(-ASV, names_to = "method", values_to = "confidence")

df.phylum <- dplyr::left_join(tax_long, conf_long, by = c("ASV", "method"))%>%
  relocate(sequence, .after=confidence)

tax_long <- class.tax %>%
  mutate(ASV = row_number()) %>%
  mutate(sequence = as.character(dna)) %>% 
  pivot_longer(-c(ASV,sequence), names_to = "method", values_to = "taxon") %>% 
  mutate(level="class")

conf_long <- class.conf %>%
  mutate(ASV = row_number()) %>%
  pivot_longer(-ASV, names_to = "method", values_to = "confidence")

df.class <- dplyr::left_join(tax_long, conf_long, by = c("ASV", "method"))%>%
  relocate(sequence, .after=confidence)

tax_long <- order.tax %>%
  mutate(ASV = row_number()) %>%
  mutate(sequence = as.character(dna)) %>% 
  pivot_longer(-c(ASV,sequence), names_to = "method", values_to = "taxon") %>% 
  mutate(level="order")

conf_long <- order.conf %>%
  mutate(ASV = row_number()) %>%
  pivot_longer(-ASV, names_to = "method", values_to = "confidence")

df.order <- dplyr::left_join(tax_long, conf_long, by = c("ASV", "method"))%>%
  relocate(sequence, .after=confidence)

tax_long <- family.tax %>%
  mutate(ASV = row_number()) %>%
  mutate(sequence = as.character(dna)) %>% 
  pivot_longer(-c(ASV,sequence), names_to = "method", values_to = "taxon") %>% 
  mutate(level="family")

conf_long <- family.conf %>%
  mutate(ASV = row_number()) %>%
  pivot_longer(-ASV, names_to = "method", values_to = "confidence")

df.family <- dplyr::left_join(tax_long, conf_long, by = c("ASV", "method"))%>%
  relocate(sequence, .after=confidence)

tax_long <- genus.tax %>%
  mutate(ASV = row_number()) %>%
  mutate(sequence = as.character(dna)) %>% 
  pivot_longer(-c(ASV,sequence), names_to = "method", values_to = "taxon") %>% 
  mutate(level="genus")

conf_long <- genus.conf %>%
  mutate(ASV = row_number()) %>%
  pivot_longer(-ASV, names_to = "method", values_to = "confidence")

df.genus <- dplyr::left_join(tax_long, conf_long, by = c("ASV", "method"))%>%
  relocate(sequence, .after=confidence)



tse_sums <- data.frame(
  sequence = rowData(TSE_16S)$sequence,
  total_abundance = rowSums(assay(TSE_16S))
)


DF <- bind_rows(df.genus, df.family, df.order, df.class, df.phylum, df.kingdom) %>%
  separate(method, into = c("method", "database"), sep = "_") %>% 
  mutate(
    usable_taxon = case_when(
      is.na(taxon) ~ FALSE,
      taxon == "" ~ FALSE,
      grepl("Incertae", taxon) ~ FALSE,
      TRUE ~ TRUE
    ),
    confidence_effective = ifelse(usable_taxon, confidence, NA_real_)
  ) %>%
  mutate(
    taxon_type = case_when(
      is.na(taxon) ~ "missing",
      taxon == "" ~ "missing",
      grepl("Incertae", taxon) ~ "incertae_sedis",
      TRUE ~ "named"
    )
  )  %>%
  dplyr::left_join(tse_sums, by = "sequence")








DF %>%
  mutate(assigned = !is.na(taxon)) %>%
  group_by(level, method, database) %>%
  summarise(
    assigned_frac = mean(assigned),
    mean_conf = mean(confidence, na.rm = TRUE),
    .groups = "drop"
  )


wide_80 <- DF %>%
  filter(confidence_effective > 80) %>%
  select(ASV, method, database, level, taxon) %>%
  unite(method_db, method, database) %>%
  pivot_wider(names_from = method_db, values_from = taxon)


my_countfunction<-function(LEVEL,DATABASE,METHOD,CONF=80,strict_family=F)
{  
  if(strict_family==F){
  N_identifications <- DF %>% 
  filter(level==LEVEL, usable_taxon==T, confidence_effective>CONF, database==DATABASE, method==METHOD) %>% 
  select(ASV) %>% n_distinct()
  
  total_reads <- DF %>% 
    filter(level==LEVEL, usable_taxon==T, confidence_effective>CONF, database==DATABASE, method==METHOD) %>% 
    pull(total_abundance) %>% sum()  
  
  }else{
  N_identifications <- DF %>% 
  filter(level==LEVEL, usable_taxon==T, confidence_effective>CONF, database==DATABASE, method==METHOD) %>% 
  filter(grepl("aceae",taxon))%>%  
  select(ASV) %>% n_distinct()  
  
  total_reads <- DF %>% 
    filter(level==LEVEL, usable_taxon==T, confidence_effective>CONF, database==DATABASE, method==METHOD) %>% 
    filter(grepl("aceae",taxon))%>%  
    pull(total_abundance) %>% sum()  
  }
  
  return(list("perc_asv" = N_identifications/n_distinct(DF$ASV), "perc_reads" = total_reads/sum(tse_sums$total_abundance)))
}


my_countfunction.plot <- function(taxonlevel,strict_family=F){

cutoff_vec <- 30:100
rates_gg2     <- NULL
rates_silva   <- NULL
rates_gtdb    <- NULL
rates_refseq  <- NULL
rates_rdp     <- NULL

reads_gg2     <- NULL
reads_silva   <- NULL
reads_gtdb    <- NULL
reads_refseq  <- NULL
reads_rdp     <- NULL

for (i in cutoff_vec){
  xgg2    <- my_countfunction(taxonlevel,"gg2","dada",i,strict_family)   
  xsilva  <- my_countfunction(taxonlevel,"silva","dada",i,strict_family) 
  xgtdb   <- my_countfunction(taxonlevel,"gtdb","dada",i,strict_family)  
  xrefseq <- my_countfunction(taxonlevel,"refseq","dada",i,strict_family)
  xrdp    <- my_countfunction(taxonlevel,"rdp","dada",i,strict_family)     
  
  rates_gg2    <- c(rates_gg2,     xgg2$perc_asv   )
  rates_silva  <- c(rates_silva,   xsilva$perc_asv )
  rates_gtdb   <- c(rates_gtdb,    xgtdb$perc_asv  )
  rates_refseq <- c(rates_refseq,  xrefseq$perc_asv)
  rates_rdp    <- c(rates_rdp,     xrdp$perc_asv   )
  
  reads_gg2    <- c(reads_gg2,     xgg2$perc_reads   )
  reads_silva  <- c(reads_silva,   xsilva$perc_reads )
  reads_gtdb   <- c(reads_gtdb,    xgtdb$perc_reads  )
  reads_refseq <- c(reads_refseq,  xrefseq$perc_reads)
  reads_rdp    <- c(reads_rdp,     xrdp$perc_reads   )  
}

par(mfrow=c(1,2))
plot(rates_gtdb ~cutoff_vec, type="l",ylim=c(0,1),main=taxonlevel,xlab="bootstrap cutoff",lwd=2,ylab="proportion of ASVs")
lines(rates_silva ~cutoff_vec,col="darkblue",lwd=2)
lines(rates_gg2 ~cutoff_vec,col="darkred",lwd=2)
lines(rates_refseq ~cutoff_vec,col="darkgreen",lwd=2)
lines(rates_rdp ~cutoff_vec,col="purple",lwd=2)

legend("bottomleft",c("gtdb","silva","gg2","refseq","rdp"),col=c("black","darkblue","darkred","darkgreen","purple"),pch=15)

plot(reads_gtdb ~cutoff_vec, type="l",ylim=c(0,1),main=taxonlevel,xlab="bootstrap cutoff",lwd=2, ylab="proportion of reads")
lines(reads_silva ~cutoff_vec,col="darkblue",lwd=2)
lines(reads_gg2 ~cutoff_vec,col="darkred",lwd=2)
lines(reads_refseq ~cutoff_vec,col="darkgreen",lwd=2)
lines(reads_rdp ~cutoff_vec,col="purple",lwd=2)

legend("bottomleft",c("gtdb","silva","gg2","refseq","rdp"),col=c("black","darkblue","darkred","darkgreen","purple"),pch=15)

}


my_countfunction.plot("genus")
my_countfunction.plot("family")
my_countfunction.plot("family",strict_family = T)
my_countfunction.plot("order")
my_countfunction.plot("class")







gg2_tax <- DF %>%
  filter(database == "gg2", method == "dada") %>%
  select(sequence, level, taxon, confidence, usable_taxon, confidence_effective, taxon_type)

gg2_tax_wide <- gg2_tax %>%
  select(sequence, level, taxon) %>%
  mutate(level = level) %>%
  pivot_wider(
    names_from = level,
    values_from = taxon
  )

gg2_conf_wide <- gg2_tax_best %>%
  select(sequence, level, confidence) %>%
  mutate(level = paste0("gg2_", level, "_confidence")) %>%
  pivot_wider(
    names_from = level,
    values_from = confidence
  )

gg2_all_wide <- gg2_tax_wide %>%
  left_join(gg2_conf_wide, by = "sequence")

rd <- as.data.frame(rowData(TSE_16S_silva))
rd2 <- rd %>%
  left_join(gg2_all_wide, by = "sequence")

rowData(TSE_16S_silva) <- DataFrame(rd2, row.names = rownames(TSE_16S_silva))


rd <- rowData(TSE_16S_silva)
keep <- (is.na(rd$family) | rd$family != "Mitochondria") & (is.na(rd$class)  | rd$class  != "Chloroplast")
TSE_16S_silva_filt <- TSE_16S_silva[keep, ]

tree <- rowTree(TSE)
# keep only ASVs present in TSE
tree_pruned <- keep.tip(tree, rownames(TSE))
# sanity checks
length(tree_pruned$tip.label)
setdiff(tree_pruned$tip.label, rownames(TSE))  # should be empty
setdiff(rownames(TSE), tree_pruned$tip.label)  # should be empty
sum(tree_pruned$edge.length == 0, na.rm = TRUE)
any(tree_pruned$edge.length < 0)

rowTree(TSE) <- tree_pruned

tree_sub <- ladderize(rowTree(TSE))
plot(tree_sub, cex = 0.4)



# Unifrac for silva:
TSE <- addDissimilarity(
  x = TSE,
  assay.type = "counts", # Whether we use abundance or counts, it doesn't matter because the unifrac function uses normalized abundances anyway (in case of weighted)
  method = "unifrac",
  weighted = TRUE,
  name = "w_unifrac"
)

TSE <- addDissimilarity(
  x = TSE,
  assay.type = "counts",
  method = "unifrac",
  weighted = F,
  name = "u_unifrac_counts"
)

TSE <- addDissimilarity(
  x = TSE,
  assay.type = "counts",
  method = "unifrac",
  weighted = FALSE,
  niter = 100,      # Averages distances over 100 subsampling iterations
  sample = 10000,   # Set your target depth here
  name = "u_unifrac_rarefied"
)



saveRDS(TSE, "results/objects/TSE_16S_taxonomy_gg2_phylogeny_silva.rds")





pcoa_w1<- cmdscale(as.dist(metadata(TSE)$w_unifrac), k = 2, eig = TRUE)
pcoa_u1<- cmdscale(as.dist(metadata(TSE)$u_unifrac_counts), k = 2, eig = TRUE)
pcoa_u2<- cmdscale(as.dist(metadata(TSE)$u_unifrac_rarefied), k = 2, eig = TRUE)

df_w1 <- data.frame(Axis1 = pcoa_w1$points[, 1], Axis2 = pcoa_w1$points[, 2], sampling_site = colData(TSE)$sampling_site,sample_id = colnames(TSE))
df_u1 <- data.frame(Axis1 = pcoa_u1$points[, 1], Axis2 = pcoa_u1$points[, 2], sampling_site = colData(TSE)$sampling_site,sample_id = colnames(TSE))
df_u2 <- data.frame(Axis1 = pcoa_u2$points[, 1], Axis2 = pcoa_u2$points[, 2], sampling_site = colData(TSE)$sampling_site,sample_id = colnames(TSE))

var_w1 <- pcoa_w1$eig / sum(pcoa_w1$eig[pcoa_w1$eig > 0])
var_u1 <- pcoa_u1$eig / sum(pcoa_u1$eig[pcoa_u1$eig > 0])
var_u2 <- pcoa_u2$eig / sum(pcoa_u2$eig[pcoa_u2$eig > 0])

p_w1<- ggplot(df_w1, aes(Axis1, Axis2, color = sampling_site)) +
  geom_point(size = 2.5) +
  labs(title = "Weighted UniFrac - SEPP SILVA128", x = paste0("PCoA1 (", round(100 * var_w1[1], 1), "%)"), y = paste0("PCoA2 (", round(100 * var_w1[2], 1), "%)")) +
  theme_bw()

p_u1<- ggplot(df_u1, aes(Axis1, Axis2, color = sampling_site)) +
  geom_point(size = 2.5) +
  labs(title = "Unweighted UniFrac (counts) - SEPP SILVA128", x = paste0("PCoA1 (", round(100 * var_u1[1], 1), "%)"), y = paste0("PCoA2 (", round(100 * var_u1[2], 1), "%)")) +
  theme_bw()

p_u2<- ggplot(df_u2, aes(Axis1, Axis2, color = sampling_site)) +
  geom_point(size = 2.5) +
  labs(title = "Unweighted UniFrac (rarefied) - SEPP SILVA128", x = paste0("PCoA1 (", round(100 * var_u2[1], 1), "%)"), y = paste0("PCoA2 (", round(100 * var_u2[2], 1), "%)")) +
  theme_bw()


library(cowplot)
library(patchwork)

leg <- get_legend(
  p_w1 + theme(legend.position = "right")
)

p_u1_noleg <- p_u1 + theme(legend.position = "none", plot.title = element_text(size = 10))
p_u2_noleg <- p_u2 + theme(legend.position = "none", plot.title = element_text(size = 10))
p_w1_noleg <- p_w1 + theme(legend.position = "none", plot.title = element_text(size = 10))

(p_u1_noleg | p_u2_noleg) /
  (p_w1_noleg | wrap_elements(leg))















DF %>% 
  filter(database=="gg2") %>%
  filter(level=="family") %>%
  select(taxon) %>%
  table()








































phylum_list <- list()
wide_phylum <- wide[wide$level=="phylum",]
for(i in unique(wide_phylum$dada_silva))
{
  slice      <- wide_phylum[wide_phylum$dada_silva == i,]
  tablenames <- sort(table(as.vector(as.matrix(slice[,4:9]))), decreasing = TRUE)
  phylum_list[[i]] <- tablenames
  #tablenames <- tablenames[names(tablenames) != i]
}  



# Let's follow this simple routine: 
# For each ASV, we go from kingdom down to genus
# At each stage, if silva has an effective confidence above a certain value, 
# we keep the taxon, and the silva assignment. 
# If not, we look whether the other databases do have an effective confidence above the threshold
# and if so, we keep those.


for(asv in unique(DF$ASV))
{
  x <- DF[DF$ASV == asv,]
  # kingdom_level
  x_kingdom <- x[x$level == "kingdom",]
  if(n_distinct(x_kingdom$level)==1){x_kingdom <- x_kingdom[1,]}else{(print(asv))}
  # phylum_level
  x_phylum <- x[x$level == "phylum",]
  if(n_distinct(x_phylum$taxon)==1){x_phylum <- x_phylum[1,] }else{(print(asv))}
  # class_level
}

































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






# Let's select 1000 non-redundant representative sequences:
dna       <- DNAStringSet(rowData(TSE_16S)$sequence)
km        <- oligonucleotideFrequency(dna, width = 6)
d         <- dist(km, method = "euclidean")
hc        <- hclust(d, method = "average")
cl        <- cutree(hc, k = 1000)


REPRESENTATIVES <- which(!duplicated(cl))


# taxonomy down to genus level, silva
TSE_16S_tax_silva <- assignTaxonomy(
  refFasta = "db/silva_nr99_v138.2_toGenus_trainset.fa.gz",
  seqs = rowData(TSE_16S)$sequence[REPRESENTATIVES],
  taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  multithread = 16,
  tryRC = TRUE,
  verbose = TRUE,
  minBoot = 30,
  outputBootstraps=T
)



# taxonomy down to genus level, gtdb
TSE_16S_tax_gtdb <- assignTaxonomy(
  refFasta = "db/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz",
  seqs = rowData(TSE_16S)$sequence[REPRESENTATIVES],
  taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  multithread = 16,
  tryRC = TRUE,
  verbose = TRUE,
  minBoot = 30,
  outputBootstraps=T
)


# taxonomy down to genus level, rdp
TSE_16S_tax_rdp <- assignTaxonomy(
  refFasta = "db/rdp_19_toGenus_trainset.fa.gz",
  seqs = rowData(TSE_16S)$sequence[REPRESENTATIVES],
  taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  multithread = 16,
  tryRC = TRUE,
  verbose = TRUE,
  minBoot = 30,
  outputBootstraps=T
)

# taxonomy down to genus level, refseq
TSE_16S_tax_refseq <- assignTaxonomy(
  refFasta = "db/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz",
  seqs = rowData(TSE_16S)$sequence[REPRESENTATIVES],
  taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  multithread = 16,
  tryRC = TRUE,
  verbose = TRUE,
  minBoot = 30,
  outputBootstraps=T
)


dna <- DNAStringSet(rowData(TSE_16S)$sequence[REPRESENTATIVES]) # Create a DNAStringSet from the ASVs
load("db/SILVA_SSU_r138.2_v2.RData")
ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=FALSE, threshold=30)
ranks <- c("domain", "phylum", "class", "order", "family", "genus") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid_silva <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
taxConf_silva <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$confidence[m]
  taxa
}))


load("db/RDP_TrainingSet_v18-mod.RData")
ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=FALSE, threshold=30, )
ranks <- c("domain", "phylum", "class", "order", "family", "genus") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid_rdp <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
taxConf_rdp <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$confidence[m]
  taxa
}))


load("db/GTDB_r226_classifier.RData")
ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=FALSE, threshold=30, )
ranks <- c("domain", "phylum", "class", "order", "family", "genus") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid_gtdb <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
taxConf_gtdb <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$confidence[m]
  taxa
}))






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



# I did a BLASTn search and got the top 20 results concatenated for all 160 query asv's that came out of the above. 
blasthits <- read.table("results/fasta/blast_summary_table.tsv",header=T,sep="\t", quote="")
blasthits <- blasthits[blasthits$Percent_identity>95,]

hits1 <- rep(F,nrow(blasthits))
hits2 <- rep(F,nrow(blasthits))
for (i in 1 : nrow (blasthits)){
  line <- blasthits[i,]
  foo <- unlist(gsub(".__","",unlist(line[2:8]))); foo <- foo[!is.na(foo)]
  hits1[i] <- grepl(tail(foo,1), line$Description, ignore.case=T)
  hits2[i] <- grepl(tail(foo,2)[1], line$Description, ignore.case=T)
  }


blasthits<- cbind.data.frame(blasthits, hits1, hits2)
bar <- blasthits[hits1|hits2,]

# In this manner, we find the name of the identification in 33 out of 160 asv's. This is just an indication. 
# It works a lot better in ITS (22) then in Bacteria (11)!
# If you look at less specific taxonomic levels, you see even more hits, especially also in bacteria (18). 
# In the class of Agaricomycetes, we find 19 hits (out of 80 ASV's)









### Next step: DECIPHER

load("db/SILVA_SSU_r138_2_2024.RData")
trainingset_silva <- trainingSet
load("db/UNITE_vAll_April2024.RData")
trainingset_unite <- trainingSet


dna_16S <- DNAStringSet(rowData(TSE_16S)$sequence)
tax_16S_decipher <- IdTaxa(dna_16S, trainingset_silva, strand="both", threshold=60, verbose=TRUE, processors = 15)







### CONTINUATION OF 7 April !!!

expected_seqs <- readDNAStringSet("data/raw/mockbacteria.fasta")
asv_seqs      <- DNAStringSet(rowData(TSE_16S_mock)$sequence)

# Taxonomy info
tax_df <- as.data.frame(rowData(TSE_16S_mock)[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")])
asv_ids <- rownames(TSE_16S_mock)

# Abundance
counts <- rowSums(assay(TSE_16S_mock))

# Label extractor
labels <- apply(tax_df, 1, function(x) {
  if (!is.na(x["Species"]) && !is.na(x["Genus"])) {
    paste0(x["Genus"], " ", x["Species"])
  } else {
    last_non_na <- tail(na.omit(x), 1)
    if (length(last_non_na) > 0) last_non_na else "Unclassified"
  }
})

# Final names
names(asv_seqs) <- paste0(asv_ids, " | ", labels, " | ", counts, " reads")


FWD <- "CCTACGGGNGGCWGCAG"
REV <- reverseComplement(DNAString("GACTACHVGGGTATCTAATCC"))
fwd_hits <- vmatchPattern(FWD, expected_seqs, fixed = F, max.mismatch=3)
rev_hits <- vmatchPattern(REV, expected_seqs, fixed = F, max.mismatch=3)

has_both <- (elementNROWS(fwd_hits) > 0) & (elementNROWS(rev_hits) > 0)

trimmed_expected <- NULL
for (i in which(has_both)) {
start_pos <- start(fwd_hits[[i]])[1]
end_pos <- end(rev_hits[[i]])[1]
if(!is.na(start_pos) && !is.na(end_pos) && end_pos > start_pos) {
trimmed_expected <- c(trimmed_expected, DNAStringSet(subseq(expected_seqs[i], start = start_pos+nchar(FWD), end = end_pos-nchar(REV))))
}}

trimmed_expected <- do.call(c, trimmed_expected)

# Now the fasta that we used for expected sequences contained many versions
# But after cutting out the region we are interested in, most look the same. 
# Hence, we collapse identical sequences in trimmed_expected
# Collapse identical sequences in trimmed_expected
unique_seqs <- unique(trimmed_expected)  # returns unique sequences
# Find indices of the first occurrence of each unique sequence
first_occurrences <- match(unique_seqs, trimmed_expected)
# Keep names from first occurrence
names(unique_seqs) <- names(trimmed_expected)[first_occurrences]
# Replace trimmed_expected with deduplicated version
trimmed_expected <- unique_seqs

N_expected <- length(trimmed_expected)

all_seqs <- c(trimmed_expected, asv_seqs)  # Both should be DNAStringSet object

# Step 1: Compute k-mer distance (k = 4)
kmer_counts <- oligonucleotideFrequency(all_seqs, width = 4)
kmer_dist   <- dist(kmer_counts, method = "manhattan")
attr(kmer_dist, "Labels") <- names(all_seqs)  # Crucial step!

# 4. Hierarchical clustering
#hc <- hclust(kmer_dist, method = "average")
nj <- ape::nj(kmer_dist)

# 5. Convert to ape::phylo object
phylo_tree <- as.phylo(nj)
tipcolors  <- rep("darkgreen",length(all_seqs)) 
tipcolors[1:10] <- "darkred"

# Match tip labels to ASV IDs
tip_labels <- phylo_tree$tip.label
abundance_vec <- counts[match(tip_labels, names(asv_seqs))]  # NA for expected seqs

# Expected sequences get highest size
is_expected <- is.na(abundance_vec)
abundance_vec[is_expected] <- max(counts, na.rm = TRUE)

# Scale for label size
scaled_cex <- log1p(abundance_vec)  # log scale
scaled_cex <- 0.4 + 0.6 * (scaled_cex / max(scaled_cex))  # range ~0.4–1.0




pdf("results/figures/03_mock_alignment_tree.pdf", width = 10, height = 20)
plot(
  phylo_tree,
  type = "phylogram",
  direction = "rightwards",
  cex = scaled_cex,
  label.offset = 1,
  tip.color = tipcolors,
  edge.width = 1  # uniform edge width again
)
title("Mock ASVs and Expected Sequences\n(Label size = log abundance)")
dev.off()


plot(
  phylo_tree,
  type = "phylogram",
  direction = "rightwards",
  cex = 0.7,
  label.offset = 1,
  tip.color = tipcolors,
  edge.width = edge_lwd
)


plot(
  phylo_tree,
  type = "phylogram",
  direction = "rightwards",
  cex = 0.7,
  label.offset = 1,
  tip.color = tipcolors
)
title("Mock ASVs and Expected Sequences")
dev.off()


# Step 3: Reorder the sequences based on clustering
ordered_seqs <- all_seqs[tree$order]

# Step 4: Perform multiple sequence alignment
alignment <- DECIPHER::AlignSeqs(ordered_seqs, verbose = TRUE)

# Optional: Save to FASTA if needed
Biostrings::writeXStringSet(alignment, filepath = "results/fasta/mock_alignment_ordered.fasta")


plot(tree)


rownames(dist_matrix) <- names(trimmed_expected)
colnames(dist_matrix) <- rownames(TSE_16S_mock)


apply(dist_matrix, 1, function(x) {
  best_idx <- which.min(x)
  paste0("Best match: ", names(x)[best_idx], " (dist = ", round(x[best_idx], 3), ")")
})


expected_kmers_counts <- oligonucleotideFrequency(trimmed_expected, width = 4, as.prob = FALSE)
observed_kmers_counts <- oligonucleotideFrequency(observed_asvs, width = 4, as.prob = FALSE)

dist_matrix_raw <- as.matrix(proxy::dist(observed_kmers_counts, observed_kmers_counts, method = "Manhattan"))
all_seqs <- c(trimmed_expected, observed_asvs)  # Both should be DNAStringSet object
names(observed_asvs) <- rownames(TSE_16S_mock)
fwd_start <- as.integer(start(fwd_hits)[,1])
rev_end   <- as.integer(end(rev_hits)[,1])




fwd_start <- as.integer(start(fwd_hits)[,1])+nchar(FWD)
rev_end   <- as.integer(end(rev_hits)[,1])-nchar(REV)


dists <- DistanceMatrix(alignment, type = "dist")

hc <- hclust(dists, method = "complete")
 ordered_indices <- hc$order
all_seqs_ordered <- all_seqs[ordered_indices]

alignment <- DECIPHER::AlignSeqs(all_seqs_ordered)
names(alignment) <- names(all_seqs_ordered)

Biostrings::writeXStringSet(alignment, "results/fasta/aligned_ordered_sequences.fasta")
















Next steps: 
  * Carefully inspect the BLAST outputs
  * Next decide which level of stringency. For ITS I think 30, but for bacteria perhaps 50
  * I can always add a column with a warning for the 30pct. Like a column that says the level of taxonomy down to which it went with 50pct identification
  * Add identifications from other sources? Flag them? 
  * Check other kingdoms in the ITS barcodes, perhaps kick them out
  * Look for AMF
  * remove plastids from bacteria
  * Phylogeny: guided or de novo. 


plastid_hits <- grepl("Chloroplast|Mitochondria", rowData(TSE_16S)$Order, ignore.case = TRUE) |
  grepl("Chloroplast|Mitochondria", rowData(TSE_16S)$Family, ignore.case = TRUE)

TSE_16S_plastids <- TSE_16S[plastid_hits, ]
TSE_16S_filtered <- TSE_16S[!plastid_hits, ]
plastid_counts_per_lib <- colSums(assay(TSE_16S[plastid_hits, ]))
plastid_pct_per_lib <- colSums(assay(TSE_16S[plastid_hits, ])) / colSums(assay(TSE_16S)) * 100



Check for perfect substrings:
  
  
  indices <- order(width(all_seqs))

SEQS    <- all_seqs[indices]
CLUST   <- ktest6$cluster[indices]

HITS    <- list()
LENGTHS <- list()

for (i in 1:length(SEQS))
{
  whichclust     <- CLUST[i]
  target_indices <- CLUST==whichclust & c(1:length(SEQS))>i 
  targets        <- SEQS[target_indices]
  hits           <- vcountPattern(SEQS[[i]],targets,max.mismatch=1,with.indels = T)
  
  if(sum(hits)>0){
    
    recalculated_indices <- c(indices[i],which(all_seqs %in% targets[hits>0]))
    HITS    <- append(HITS, list(recalculated_indices))
    LENGTHS <- append(LENGTHS, list(width(all_seqs[recalculated_indices]))) 
    print(i)
  }
}


