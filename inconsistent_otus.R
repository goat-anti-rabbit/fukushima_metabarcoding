library(dplyr)

# Combine OTUs with taxonomic data
otu_taxa_data <- cbind(OTU = X16S$otus97.5, as.data.frame(X16S$taxa))

# Check for inconsistencies within each OTU
inconsistent_otus <- otu_taxa_data %>%
    group_by(OTU) %>%
    summarise(across(everything(), ~ n_distinct(., na.rm = TRUE))) %>%
    filter(if_any(-OTU, ~ . > 1)) %>%
    pull(OTU)



i = 1
unname(X16S$taxa[X16S$otus97.5 == inconsistent_otus[i],])


