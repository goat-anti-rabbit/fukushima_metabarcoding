# @ Rik Verdonck 20250121

# Purpose: Perform exploratory data analyses (EDA) and preliminary visualizations.

# Alpha and beta diversity calculations.
# Visualization of sample distributions (e.g., ordination plots, rarefaction curves).
# Assess batch effects or biases in sequencing depth, metadata, etc.
# Save preliminary figures and summary tables.



# First get counts per family and genus. I remove the NA taxa here. You can also opt to include them as "unknown"

XITS$family_counts <- t(rowsum(t(XITS$x[, !is.na(XITS$taxa[, 5]) ]), group = XITS$taxa[!is.na(XITS$taxa[, 5])  ,5],na.rm=T))
XITS$genus_counts  <- t(rowsum(t(XITS$x[, !is.na(XITS$taxa[, 6]) ]), group = XITS$taxa[!is.na(XITS$taxa[, 6])  ,6],na.rm=T))


# Next I execute the pcoa
# Using the unfweighted unifrac here, but any distance measure can be used
pcoa_res <- ape::pcoa(XITS$beta_asv$k7$euc$UUF)


# Now you calculate the correlation between your family counts, and the two first principal coordinates. 
# You can do this with genus counts as well, obviously. 
cor_matrix <- cor(XITS$family_counts, pcoa_res$vectors[,1:2], method = "pearson")

# Select top 10 taxa with highest correlation to axis 1 or 2
top_taxa <- head(order(rowSums(abs(cor_matrix)), decreasing = TRUE), 10)

# These are the coordinates of your top 10 in the 2D PCoA plane:
arrow_coords <- cor_matrix[top_taxa, ]

# Plot the PCoA:
plot(pcoa_res$vectors[,1:2], col = sitecolors[XITS$META$site], pch = 19,xlab = "PCoA1", ylab = "PCoA2", main = "PCoA with genus arrow")

# You will need some scaling because the units of the arrows are a bit arbitrary
arrow_scale=0.25

# Add arrows:
arrows(0, 0, x1 = arrow_coords[,1] * arrow_scale, y1 = arrow_coords[,2] * arrow_scale, length = 0.1, col = "gray30", lwd = 1.5)

# Add labels:
text(arrow_coords[,1] * arrow_scale * 1.1, arrow_coords[,2] * arrow_scale * 1.1, labels = rownames(arrow_coords), col = "gray10", cex = 0.8)