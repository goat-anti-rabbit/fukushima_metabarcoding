FUNS$get.sample.name      <- function(fname) strsplit(basename(fname), "_")[[1]][1]



FUNS$allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)  
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  
}


FUNS$primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- Biostrings::vcountPattern(primer, ShortRead::sread(ShortRead::readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

FUNS$substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}



FUNS$generate_primer_summary <- function(dataset, fnFs_col, fnRs_col, fwd_orients_col = "FWD.orients", rev_orients_col = "REV.orients", multithread = 4) 
{
  cl <- makeCluster(multithread)
  registerDoParallel(cl)
  # Perform primer hits calculation in parallel
  primer_summary <- foreach(i = 1:length(dataset[[fnFs_col]]), .combine = rbind, .export=c("FUNS")) %dopar% {
    c( FW.fw = sapply(dataset$FWD.orients, FUNS$primerHits, fn = dataset[[fnFs_col]][[i]]),
       FW.rv = sapply(dataset$FWD.orients, FUNS$primerHits, fn = dataset[[fnRs_col]][[i]]),
       RV.fw = sapply(dataset$REV.orients, FUNS$primerHits, fn = dataset[[fnFs_col]][[i]]),
       RV.rv = sapply(dataset$REV.orients, FUNS$primerHits, fn = dataset[[fnRs_col]][[i]])
    )
  }
  
  rownames(primer_summary) <- dataset$sample.names
  stopCluster(cl)
  return(primer_summary)
}


FUNS$getN     <- function(x) sum(getUniques(x))




FUNS$otu_cluster <- function(x, k = 5, initial_h = 0.05, final_h = 0.025)
{
  distk <- kdistance(sapply(colnames(x),strsplit,split=""),k=k); cat("k distance matrix ready")
  HC    <- hclust(distk, method="average"); cat ("clustering ready")
  otus  <- cutree(HC,h=initial_h)
  for(i in which(table(otus)>1))
  {
    alignment <- msa(DNAStringSet(names(which(otus==i))))
    hc        <- hclust(as.dist(DistanceMatrix(DNAStringSet(alignment),verbose=F)))
    clusters  <- cutree(hc,h=final_h)
    if(length(unique(clusters))>1)
    {
      otus[otus == i] <- paste(i,clusters,sep="_")
    }
  }
  return(unname(otus))
}




# This function generates an NMDS based on the distance matrix provided.
# Next, it plots the data (sites) onto the two first NMDS axes. 
# Next, it draws pink ellipses around all members of a site, and red ellipses based on the standard deviation around the centroid.
FUNS$betaplot <- function(DIST, k = 2, main)
{
  #pcoa_results     <- cmdscale(DIST,k=k)
  nmds_results    <- vegan::metaMDS(DIST, k = k, trymax = 100)
  scores          <- scores(nmds_results)
  centroids       <- aggregate(cbind(scores[,1], scores[,2]), by=list(site=META[rownames(as.matrix(DIST)),"site"]), mean); names(centroids)[2:3] <- c("x", "y")
  xlim            <- max(abs(range(scores[,1])))
  ylim            <- max(abs(range(scores[,2])))
  plot(scores[,1], scores[,2], xlab = "Coordinate 1", ylab = "Coordinate 2", main = main,col="white",xlim=c(-xlim,xlim),ylim=c(-ylim,ylim))
  points(centroids$x, centroids$y, pch=4, col="red", cex=1.5)
  with(centroids, ordiellipse(scores, groups = META[rownames(as.matrix(DIST)),"site"], draw="polygon", kind = "ehull", col = sitecolors, border="white", alpha=50,lwd = 1))
  with(centroids, ordiellipse(scores, groups = META[rownames(as.matrix(DIST)),"site"], kind = "sd", col = sitecolors, lwd = 3))
  text(centroids$x, centroids$y, labels=centroids$site, pos=3, col="red")
  text(scores[,1], scores[,2], labels=rownames(scores),cex=0.75)
} 

# This function draws two boxplots, side by side, comparing the first and second axis of NMDS between sites.  
FUNS$betaboxplot <- function(DIST, k=2, main)
{
  #pcoa_results     <- cmdscale(DIST,k=k)
  nmds_results     <- vegan::metaMDS(DIST, k = k, trymax = 100)
  scores           <- scores(nmds_results)
  dat              <- merge(scores,META,by="row.names")
  par(mfrow=c(1,2))
  boxplot(dat$NMDS1 ~ dat$site, xlab= "site", ylab= "NMDS 1",main="main")
  boxplot(dat$NMDS2 ~ dat$site, xlab ="site",ylab = "NMDS 2")
} 

# This function plots measured 137Cs levels as a function of NMDS axis 1 or 2
FUNS$radiationplot <- function(DIST, k=2, main)
{
  nmds_results     <- vegan::metaMDS(DIST, k = k, trymax = 100)
  scores           <- scores(nmds_results)
  dat              <- merge(scores,META,by="row.names")
  par(mfrow=c(1,2))
  plot(dat$NMDS1, dat$Cs_137, log="y",xlab="NMDS1", ylab = "137 Cs (Bq/kg)",col = )
  plot(dat$NMDS2, dat$Cs_137, log="y",xlab="NMDS2", ylab = "137 Cs (Bq/kg)",col = )
  
}

FUNS$collapse_duplicates <- function(df) {
  # Get row names
  row_names <- rownames(df)
  
  # Create a logical vector to check if a rowname ends with '1'
  duplicate_rows <- grepl("1$", row_names)
  
  # Loop through each row and sum duplicate rows with the preceding one
  for (i in seq_len(nrow(df))) {
    if (duplicate_rows[i]) {
      # Get the name of the row without the trailing "1"
      original_row_name <- sub("1$", "", row_names[i])
      
      # Find the row corresponding to the original (without '1')
      original_row_index <- which(row_names == original_row_name)
      
      # Sum the counts from the duplicate row with the original row
      df[original_row_index, ] <- df[original_row_index, ] + df[i, ]
    }
  }
  
  # Now remove the duplicate rows that end in "1"
  df <- df[!duplicate_rows, ]
  
  return(df)
}


FUNS$cosine_dissim <- function (DF)
{
  Matrix <- as.matrix(DF)
  sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
  sim <- sim %*% t(sim)
  D_sim <- as.dist(1 - sim)
  return(D_sim)
}


FUNS$generate_beta_asv <- function(input_data, k) {
  # Ensure the input object has the required structure
  if (!is.list(input_data) || !"X" %in% names(input_data)) {
    stop("The input object must have a component named 'X'")
  }
  
  # Initialize or update the beta_asv field with the appropriate structure
  if (!"beta_asv" %in% names(input_data)) {
    input_data$beta_asv <- list()
  }
  
  k_str <- paste0("k", k)
  input_data$beta_asv[[k_str]] <- list()
  
  # K-mer counts
  kcounts <- kcount(sapply(colnames(input_data$X), strsplit, split = ""), k = k, residues = "DNA")
  input_data$beta_asv[[k_str]]$kcounts <- kcounts
  
  # Euclidean Distance
  euc_dist <- dist(kcounts, method = "euclidean")
  input_data$beta_asv[[k_str]]$euc <- list(dist = euc_dist)
  
  # kdistance
  edg_dist <- kdistance(sapply(colnames(input_data$X), strsplit, split = ""), k = k)
  input_data$beta_asv[[k_str]]$edg <- list(dist = edg_dist)
  
  # Cosine dissimilarity
  cos_dist <- FUNS$cosine_dissim(kcounts)
  input_data$beta_asv[[k_str]]$cos <- list(dist = cos_dist)
  
  # Neighbor-joining trees
  input_data$beta_asv[[k_str]]$euc$nj <- ape::nj(euc_dist)
  input_data$beta_asv[[k_str]]$edg$nj <- ape::nj(edg_dist)
  input_data$beta_asv[[k_str]]$cos$nj <- ape::nj(cos_dist)
  
  # Phyloseq objects
  input_data$beta_asv[[k_str]]$euc$phyloseq <- phyloseq(otu_table(input_data$X, taxa_are_rows = FALSE), phy_tree(input_data$beta_asv[[k_str]]$euc$nj))
  input_data$beta_asv[[k_str]]$edg$phyloseq <- phyloseq(otu_table(input_data$X, taxa_are_rows = FALSE), phy_tree(input_data$beta_asv[[k_str]]$edg$nj))
  input_data$beta_asv[[k_str]]$cos$phyloseq <- phyloseq(otu_table(input_data$X, taxa_are_rows = FALSE), phy_tree(input_data$beta_asv[[k_str]]$cos$nj))
  
  # Weighted and Unweighted UniFrac distances
  input_data$beta_asv[[k_str]]$euc$WUF <- UniFrac(input_data$beta_asv[[k_str]]$euc$phyloseq, weighted = TRUE, normalized = TRUE)
  input_data$beta_asv[[k_str]]$edg$WUF <- UniFrac(input_data$beta_asv[[k_str]]$edg$phyloseq, weighted = TRUE, normalized = TRUE)
  input_data$beta_asv[[k_str]]$cos$WUF <- UniFrac(input_data$beta_asv[[k_str]]$cos$phyloseq, weighted = TRUE, normalized = TRUE)
  
  input_data$beta_asv[[k_str]]$euc$UUF <- UniFrac(input_data$beta_asv[[k_str]]$euc$phyloseq, weighted = FALSE, normalized = TRUE)
  input_data$beta_asv[[k_str]]$edg$UUF <- UniFrac(input_data$beta_asv[[k_str]]$edg$phyloseq, weighted = FALSE, normalized = TRUE)
  input_data$beta_asv[[k_str]]$cos$UUF <- UniFrac(input_data$beta_asv[[k_str]]$cos$phyloseq, weighted = FALSE, normalized = TRUE)
  
  return(input_data)
}


FUNS$barplot_function <- function(asvtable, taxonomy, taxonomic_level = "Order", samples, other = 10, include_na = TRUE, na_as_top = FALSE) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(stringr)  # For handling "unknown" taxa pooling
  
  # Convert ASV count table to long format
  X_long <- as.data.frame(asvtable[samples, , drop = FALSE]) %>%
    tibble::rownames_to_column("Sample") %>%
    pivot_longer(cols = -Sample, names_to = "ASV", values_to = "Count")
  
  # Ensure taxonomy has ASV as a column
  taxa_df <- as.data.frame(taxonomy) %>%
    tibble::rownames_to_column("ASV")
  
  # Join ASV counts with taxonomy
  combined_data <- left_join(X_long, taxa_df, by = "ASV")
  
  # Pool "unknown" taxa with "NA"
  combined_data <- combined_data %>%
    mutate(!!sym(taxonomic_level) := if_else(
      is.na(!!sym(taxonomic_level)) | str_detect(!!sym(taxonomic_level), regex("unknown", ignore_case = TRUE)),
      "NA",
      !!sym(taxonomic_level)
    ))
  
  # Separate counts for "NA" and assigned taxa
  na_counts <- combined_data %>%
    filter(!!sym(taxonomic_level) == "NA") %>%
    group_by(Sample) %>%
    summarise(NA_Count = sum(Count), .groups = 'drop')
  
  assigned_counts <- combined_data %>%
    filter(!!sym(taxonomic_level) != "NA") %>%
    group_by(Sample, !!sym(taxonomic_level)) %>%
    summarise(Count = sum(Count), .groups = 'drop')
  
  # Identify the top taxa by total abundance
  top_taxa <- assigned_counts %>%
    ungroup() %>%
    group_by(!!sym(taxonomic_level)) %>%
    summarise(TotalCount = sum(Count)) %>%
    arrange(desc(TotalCount)) %>%
    slice_head(n = other) %>%
    pull(!!sym(taxonomic_level))
  
  # Categorize "Other"
  assigned_counts <- assigned_counts %>%
    mutate(!!sym(taxonomic_level) := if_else(!!sym(taxonomic_level) %in% top_taxa, 
                                             as.character(!!sym(taxonomic_level)), "Other")) %>%
    group_by(Sample, !!sym(taxonomic_level)) %>%
    summarise(Count = sum(Count), .groups = 'drop')
  
  # Calculate total counts per sample
  total_counts_per_sample <- combined_data %>%
    group_by(Sample) %>%
    summarise(TotalCount = sum(Count), .groups = 'drop')
  
  # Normalize counts to proportions
  if (!include_na) {
    # Exclude "NA" and normalize to assigned taxa only
    assigned_counts <- assigned_counts %>%
      group_by(Sample) %>%
      mutate(Proportion = Count / sum(Count)) %>%
      ungroup()
    combined_plot_data <- assigned_counts
  } else if (na_as_top) {
    # Normalize assigned taxa to 100% (excluding "NA")
    assigned_counts <- assigned_counts %>%
      group_by(Sample) %>%
      mutate(Proportion = Count / sum(Count)) %>%
      ungroup()
    
    # Recalculate "NA" proportion relative to assigned taxa
    na_layer <- na_counts %>%
      left_join(
        assigned_counts %>% group_by(Sample) %>% summarise(AssignedTotal = sum(Count)),
        by = "Sample"
      ) %>%
      mutate(Proportion = NA_Count / AssignedTotal) %>%
      mutate(!!sym(taxonomic_level) := "NA") %>%
      select(Sample, !!sym(taxonomic_level), Proportion)
    
    # Combine the assigned taxa and "NA" layer
    combined_plot_data <- bind_rows(assigned_counts, na_layer)
  } else {
    # Normalize everything (including "NA") to total counts
    assigned_counts <- assigned_counts %>%
      left_join(total_counts_per_sample, by = "Sample") %>%
      mutate(Proportion = Count / TotalCount) %>%
      select(-TotalCount)
    
    na_layer <- na_counts %>%
      left_join(total_counts_per_sample, by = "Sample") %>%
      mutate(Proportion = NA_Count / TotalCount) %>%
      mutate(!!sym(taxonomic_level) := "NA") %>%
      select(Sample, !!sym(taxonomic_level), Proportion)
    
    combined_plot_data <- bind_rows(assigned_counts, na_layer)
  }
  
  # Reorder taxonomic levels for stacking (reverse order for correct stacking)
  taxon_order <- combined_plot_data %>%
    group_by(!!sym(taxonomic_level)) %>%
    summarise(TotalAbundance = sum(Proportion), .groups = 'drop') %>%
    arrange(desc(TotalAbundance)) %>%
    pull(!!sym(taxonomic_level))
  
  # Ensure "Other" and "NA" are at the top
  taxon_order <- c("NA", "Other", rev(setdiff(taxon_order, c("NA", "Other"))))
  
  combined_plot_data <- combined_plot_data %>%
    mutate(!!sym(taxonomic_level) := factor(!!sym(taxonomic_level), levels = taxon_order))
  
  # Define colors for "Other" and "NA"
  custom_colors <- c("Other" = "darkgrey", "NA" = "lightgrey")
  
  # Generate the plot
  ggplot(combined_plot_data, aes(x = Sample, y = Proportion, fill = !!sym(taxonomic_level))) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
      values = c(
        custom_colors,
        setNames(
          scales::hue_pal()(length(setdiff(taxon_order, names(custom_colors)))), 
          setdiff(taxon_order, names(custom_colors))
        )
      )
    ) +
    theme_minimal() +
    labs(x = "Sample", y = "Proportion", fill = taxonomic_level) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = guide_legend(reverse = FALSE))  # Legend order matches stacking order
}





