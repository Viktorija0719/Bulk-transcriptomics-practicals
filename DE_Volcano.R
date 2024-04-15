# Load required libraries
library(DESeq2)
library(tidyverse)
library(org.Hs.eg.db)
library(EnhancedVolcano)


# Step 1: Prepare count data ----------------
# Obtain sample information from filenames
files <- list.files(path = "E:/rna_seq/count/counts", pattern = "_counts.txt$", full.names = TRUE)
sampleNames <- gsub(".*\\/(.*)_counts.txt$", "\\1", files)

# Generate sample cell line and preparation methods from file names
sampleCell_line <- ifelse(grepl("UHRR", sampleNames), "UHRR", "HBR")
samplePreps <- ifelse(grepl("Collibri", sampleNames), "Collibri", "KAPA")


# Create colData DataFrame
colData <- DataFrame(
  cell_line = factor(sampleCell_line),
  prep = factor(samplePreps),
  row.names = sampleNames
)

# Load count data
countsList <- lapply(files, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t", row.names = 1)  # Read counts while preserving gene IDs as row names
  counts <- as.numeric(df[, ncol(df)])  # Extract count data and convert to numeric
  names(counts) <- rownames(df) # Preserve gene IDs as names
  return(counts)
})


# Check for any NAs in the counts data
if (any(sapply(countsList, function(c) any(is.na(c))))) {
  stop("NA values found in the counts data. Please check the input files.")
}

# Construct countsMatrix with gene IDs as row names
countsMatrix <- do.call(cbind, countsList)

colnames(countsMatrix) <- sampleNames # Ensure column names are correctly assigned
colData <- colData[match(colnames(countsMatrix), rownames(colData)), ] # Verify and ensure colData order matches countsMatrix columns
all(colnames(countsMatrix) %in% rownames(colData)) # making sure the row names in colData matches to column names in counts_data
all(colnames(countsMatrix) == rownames(colData)) # are they in the same order?


# set the factor level
head(colData)
colData$cell_line <- factor(colData$cell_line)
colData$prep <- factor(colData$prep)


# Step 2: DESeq2 Analysis ----------------
# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = countsMatrix,
  colData = colData,
  design = ~ cell_line + prep
)


# Pre-filtering: Remove rows with low counts
dds <- dds[rowSums(counts(dds)) >= 5, ]

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract results and maintain gene IDs
res <- results(dds, contrast = c("cell_line", "UHRR", "HBR"))

# Annotate results with gene symbols
res.df <- as.data.frame(res)
rownames(res.df) <- gsub("\\..*$", "", rownames(res.df))  # Remove version numbers from Ensembl IDs
res.df$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(res.df),
                        keytype = "ENSEMBL",
                        column = "SYMBOL",
                        multiVals = "first")

# Step 3: Volcano Plot ----------------
EnhancedVolcano(res.df,
                lab = res.df$symbol,
                x = "log2FoldChange",
                y = "padj",
                title = "Volcano Plot: Collibri vs KAPA",
                xlab = "Log2 Fold Change",
                ylab = "-Log10 Adjusted P-value",
                pCutoff = 0.05,
                FCcutoff = log2(2),
                pointSize = 2.0,
                labSize = 3.0)

# Step 4: MA Plot ----------------
plotMA(res)

# Save the annotated results
write.csv(res.df, "path/to/output/annotated_deseq2_results.csv", row.names = TRUE)
