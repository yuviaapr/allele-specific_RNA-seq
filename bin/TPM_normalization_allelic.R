#!/usr/bin/env Rscript

# Script to perform TPM normalization using a table with allelic and pooled read counts.
# Usage: Rscript TPM_normalization_allelic.R table_raw_counts.txt

# Read count table with the information of all libraries.

args <- commandArgs(trailingOnly = T)

fileName <- args[1]

# Load data

all_reads <- read.table(fileName, header = T)

# Calculate TPM values.

# Step 1. Divide read counts by gene length in kilobases (RPK).

rpk <- all_reads[,3:ncol(all_reads)] / (all_reads$Length/1000)

# Split table.

rpk_allelic <- rpk[,unlist(lapply(sort(grep("Gall", colnames(rpk), invert = T, value = T)), function(x) which(colnames(rpk) == x)))]

rpk_all <- rpk[,unlist(lapply(sort(grep("Gall", colnames(rpk), value = T)), function(x) which(colnames(rpk) == x))), drop = F]

# Step 2. Sum up RPK values by sample and divide by 1,000,000 (scaling factors).

scalingFactors <- colSums(rpk_all) / 1000000

# Step 3. Divide each RPK value by its corresponding scaling factor (TPM).

tpm_all <- sweep(rpk_all, 2, scalingFactors, FUN = '/')

tpm_allelic <- sweep(rpk_allelic, 2, rep(scalingFactors, each = 2), FUN = '/')

# Include additional gene information.

tpm_all <- cbind(all_reads[,1:2], tpm_all)
tpm_allelic <- cbind(all_reads[,1:2], tpm_allelic)

# Write annotated TPM values.

write.table(tpm_all, "TPM_normalized_counts_genes_all.txt", quote = F, row.names = F, sep = "\t")
write.table(tpm_allelic, "TPM_normalized_counts_genes_by_allele.txt", quote = F, row.names = F, sep = "\t")


