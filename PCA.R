library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)

# Read in normalized counts and sample metadata
ncounts <- read.csv("E:/rna_seq/output/normalized_counts.csv", row.names = 1)
sample_metadata <- read.csv("E:/rna_seq/output/sample_metadata.csv")
de_results <- read.csv("E:/rna_seq/output/de_results_full.csv", row.names = 1)

# Adjust p-value threshold as needed
degs <- subset(de_results, padj < 0.05 & !is.na(padj))

# Subset normalized counts for DEGs only
degs_counts <- ncounts[rownames(ncounts) %in% degs$geneID,]

# Perform PCA on the counts of DEGs 
pca_res <- prcomp(t(degs_counts), scale. = TRUE)

# Extract PCA scores
scores <- as.data.frame(pca_res$x)

# Convert PCA results to a data frame
scores <- as.data.frame(pca_res$x)

scores$Sample <- rownames(scores)
sample_metadata <- rename(sample_metadata, Sample = X)
scores$Sample <- gsub("[.]", "_", scores$Sample)
sample_metadata$Sample <- gsub("-", "_", sample_metadata$Sample) 

print(head(scores$Sample))
print(head(sample_metadata$Sample))

# Join PCA scores with sample metadata
plot_data <- left_join(scores, sample_metadata, by = "Sample")

# Plot PCA
ggplot(plot_data, aes(x = PC1, y = PC2, color = cell_line, shape = prep)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Normalized Counts",
       x = paste("PC1 -", round(100 * var(pca_res$x[,1])/sum(pca_res$sdev^2), 1), "% Variance"),
       y = paste("PC2 -", round(100 * var(pca_res$x[,2])/sum(pca_res$sdev^2), 1), "% Variance")) +
  scale_color_manual(values = c("UHRR" = "blue", "HBR" = "red")) # Adjust colors as needed

