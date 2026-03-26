###############################################################################
# BULK RNA-SEQ INTEGRATION: Genotype × Ketamine (Mouse)
# Count files contain gene_id only (Ensembl). This script:
#   - Loads the two count matrices
#   - Merges them
#   - Annotates gene IDs with org.Mm.eg.db (gene_name)
#   - Runs DESeq2 interaction model
#   - Exports DEGs (annotated)
#   - Creates PCA, volcano, heatmap
###############################################################################

library(DESeq2)
library(tidyverse)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(dplyr)

setwd("C:/Users/Calcium_RN14N/Documents/Bulk_re")

###############################################################################
### 1. Load raw count tables (with gene_id as rownames)
###############################################################################

counts_drug <- read.csv("NR_Res_ketamine_readcount_genename.csv",
                        row.names = 1, check.names = FALSE)

counts_nodrug <- read.csv("NR_Res_no_ketamine.csv",
                          row.names = 1, check.names = FALSE)


# 3. Clean gene IDs (remove version numbers and handle invalid characters)


# Identify duplicated row names (gene IDs)
duplicated_rows <- rownames(counts_nodrug)[duplicated(rownames(counts_nodrug))]
duplicated_rows

# Check for empty gene IDs
empty_genes <- which(rownames(counts_nodrug) == "")
empty_genes

# Remove duplicated gene IDs and keep only the first occurrence
counts_nodrug <- counts_nodrug[!duplicated(rownames(counts_nodrug)), ]

# Verify that row names are now unique
stopifnot(length(rownames(counts_nodrug)) == length(unique(rownames(counts_nodrug))))

# Replace 'novel' gene ID with a unique identifier
rownames(counts_nodrug)[rownames(counts_nodrug) == "novel"] <- paste0("novel_", seq_along(which(rownames(counts_nodrug) == "novel")))

# Clean gene IDs by removing special characters
rownames(counts_nodrug) <- gsub("[^a-zA-Z0-9]", "_", rownames(counts_nodrug))

# Ensure row names are unique and non-empty
stopifnot(length(rownames(counts_nodrug)) == length(unique(rownames(counts_nodrug))))
stopifnot(all(rownames(counts_nodrug) != ""))

# Clean gene IDs in counts_nodrug
rownames(counts_nodrug) <- sub("\\.\\d+$", "", rownames(counts_nodrug))  # Remove version numbers
rownames(counts_nodrug) <- gsub("[^a-zA-Z0-9]", "_", rownames(counts_nodrug))  # Remove special characters

# Handle duplicates: Keep only the first occurrence of duplicated gene IDs
counts_nodrug <- counts_nodrug[!duplicated(rownames(counts_nodrug)), ]

# Replace any non-standard 'novel' gene ID with a unique identifier
rownames(counts_nodrug)[rownames(counts_nodrug) == "novel"] <- paste0("novel_", seq_along(which(rownames(counts_nodrug) == "novel")))

# Ensure row names are unique and non-empty
stopifnot(length(rownames(counts_nodrug)) == length(unique(rownames(counts_nodrug))))
stopifnot(all(rownames(counts_nodrug) != ""))

# Now you should be able to proceed with setting row names and the rest of the analysis


# Remove version numbers (e.g., ENSMUSG00000064351.1 becomes ENSMUSG00000064351)
rownames(counts_nodrug) <- sub("\\.\\d+$", "", rownames(counts_nodrug))
rownames(counts_drug) <- sub("\\.\\d+$", "", rownames(counts_drug))

# Clean gene IDs by removing non-alphanumeric characters
rownames(counts_nodrug) <- gsub("[^a-zA-Z0-9]", "_", rownames(counts_nodrug))
rownames(counts_drug) <- gsub("[^a-zA-Z0-9]", "_", rownames(counts_drug))

# Check for empty or missing gene IDs and remove any rows that have them
counts_nodrug <- counts_nodrug[rownames(counts_nodrug) != "", ]
counts_drug <- counts_drug[rownames(counts_drug) != "", ]
counts_nodrug <- counts_nodrug[!is.na(rownames(counts_nodrug)), ]
counts_drug <- counts_drug[!is.na(rownames(counts_drug)), ]

# Ensure row names are unique
stopifnot(length(rownames(counts_nodrug)) == length(unique(rownames(counts_nodrug))))
stopifnot(length(rownames(counts_drug)) == length(unique(rownames(counts_drug))))

# Identify the common genes between counts_nodrug and counts_drug
common_genes <- intersect(rownames(counts_nodrug), rownames(counts_drug))

# Subset both datasets to only include the common genes
counts_nodrug <- counts_nodrug[common_genes, , drop = FALSE]
counts_drug <- counts_drug[common_genes, , drop = FALSE]

# Check that both datasets now have the same number of rows
dim(counts_nodrug)  # Should have the same number of rows
dim(counts_drug)    # Should also have the same number of rows

# 4. Merge the count data for both conditions (no ketamine + ketamine)
counts_all <- cbind(counts_nodrug, counts_drug)

# 5. Create metadata for the samples
sample_names <- colnames(counts_all)

# Create metadata with conditions
metadata <- data.frame(
  sample = sample_names,
  genotype = ifelse(grepl("NR", sample_names), "CON", "MUT"),  # NR indicates control group
  drug = ifelse(grepl("k_", sample_names), "KET", "NO_KET"),  # k_ indicates ketamine group
  batch = ifelse(grepl("Res", sample_names), "exp2", "exp1")  # Assuming Res is in batch 2
)

# Set sample names as rownames for metadata
rownames(metadata) <- metadata$sample

# Ensure metadata rows correspond to the count columns
stopifnot(all(rownames(metadata) == colnames(counts_all)))

###############################################################################
### Annotate gene IDs using org.Mm.eg.db (Ensembl → Symbol)
###############################################################################

gene_ids <- rownames(counts_all)

# Remove version again to be safe
gene_ids_clean <- sub("\\.\\d+$", "", gene_ids)

# Map to SYMBOL
gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = gene_ids_clean,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Create a lookup table
gene_annot <- data.frame(
  gene_id = gene_ids_clean,
  gene_symbol = gene_symbols,
  stringsAsFactors = FALSE
)

# Replace NAs with gene_id
gene_annot$gene_symbol[is.na(gene_annot$gene_symbol)] <- gene_annot$gene_id[is.na(gene_annot$gene_symbol)]

# Store mapping in rownames
rownames(counts_all) <- gene_annot$gene_id

# 6. Create DESeq2 dataset
# Create DESeq2 dataset without the 'batch' term
# Remove batch term if not necessary (simplified model)
dds <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData = metadata,
  design = ~ genotype * drug  # Interaction term between genotype and drug
)

# Re-level factors to ensure proper reference levels
metadata$genotype <- factor(metadata$genotype, levels = c("CON", "MUT"))
metadata$drug <- factor(metadata$drug, levels = c("NO_KET", "KET"))

# Recreate DESeq2 dataset after setting levels explicitly
dds <- DESeqDataSetFromMatrix(
  countData = counts_all,
  colData = metadata,
  design = ~ genotype * drug  # Use the simplified model if necessary
)

# Check model matrix if needed
model_matrix <- model.matrix(~ genotype * drug, data = metadata)
head(model_matrix)


# 7. Run DESeq2 analysis
dds <- DESeq(dds)

# 8. Extract DEGs (Differentially Expressed Genes)
# Interaction term: Genotype * Drug interaction
res_interaction <- results(dds, name = "genotypeMUT.drugKET")
###############################################################################
### 8. Extract DEGs (with gene SYMBOL added)
###############################################################################

res_interaction <- results(dds, name = "genotypeMUT.drugKET")
res_interaction <- as.data.frame(res_interaction)

# Add gene_id column explicitly
res_interaction$gene_id <- rownames(res_interaction)

# Add gene_symbol using annotation table
res_interaction <- left_join(
  res_interaction,
  gene_annot,
  by = "gene_id"
)

# Reorder columns
res_interaction <- res_interaction %>%
  relocate(gene_id, gene_symbol)

write.csv(res_interaction,
          "DEG_genotype_drug_interaction_annotated.csv",
          row.names = FALSE)

# 9. Export results to CSV
write.csv(res_interaction, "DEG_genotype_drug_interaction.csv", row.names = TRUE)

# 10. PCA Plot
vsd <- vst(dds)
pdf("PCA_plot.pdf", width = 6, height = 5)
plotPCA(vsd, intgroup = c("genotype", "drug"))
dev.off()

###############################################################################
###############################################################################
### 11. Volcano Plot with gene name labels (robust type-safe version)
###############################################################################

# If not installed:
# install.packages("ggrepel")
library(ggplot2)
library(ggrepel)
library(dplyr)

# Start from plain base data.frame to avoid S4/Rle issues
res_plot <- as.data.frame(res_interaction, stringsAsFactors = FALSE)

# Ensure required columns are plain numeric
res_plot$log2FoldChange <- as.numeric(res_plot$log2FoldChange)
res_plot$pvalue         <- as.numeric(res_plot$pvalue)
res_plot$padj           <- as.numeric(res_plot$padj)

# If gene_symbol missing/NA, fall back to gene_id for labels
if (!"gene_symbol" %in% colnames(res_plot)) res_plot$gene_symbol <- res_plot$gene_id
res_plot$gene_symbol[is.na(res_plot$gene_symbol) | res_plot$gene_symbol == ""] <- res_plot$gene_id[is.na(res_plot$gene_symbol) | res_plot$gene_symbol == ""]

# Compute -log10 adjusted p
res_plot$log10padj <- -log10(res_plot$padj)

# Clean rows for plotting/filtering
res_plot <- res_plot %>%
  filter(!is.na(log2FoldChange), !is.na(padj), is.finite(log10padj))

# Label genes: padj < 0.05 and |log2FC| > 0.2 (top 20 by padj)
label_genes <- res_plot %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.2) %>%
  arrange(padj) %>%
  slice_head(n = 20)

# Optional significance class
res_plot <- res_plot %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange >  0.2 ~ "Up",
    padj < 0.05 & log2FoldChange < -0.2 ~ "Down",
    TRUE ~ "NS"
  ))

pdf("Volcano_interaction.pdf", width = 7, height = 6)
ggplot(res_plot, aes(x = log2FoldChange, y = log10padj)) +
  geom_point(aes(color = sig), alpha = 0.7, size = 1.2) +
  scale_color_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey70")) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = label_genes,
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = 200
  ) +
  labs(
    title = "Genotype × Ketamine Interaction",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = ""
  ) +
  theme_classic(base_size = 12)
dev.off()

