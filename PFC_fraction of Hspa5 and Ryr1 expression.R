library(Seurat)
library(ggplot2)

# -----------------------------
# 1️⃣ Setup Genes and Fetch Data
# -----------------------------
genes_to_check <- c("Hspa5", "Ryr1")

# Fetch all data at once (Safe and fast)
# This creates a table with: celltype, Hspa5 counts, Ryr1 counts
data_table <- FetchData(PFC1, vars = c("celltype", genes_to_check))

# -----------------------------
# 2️⃣ Define Calculation Function
# -----------------------------
# We make a simple function to calculate stats for any specific gene
calculate_stats <- function(df, gene_name) {
  # Sum the specific gene by celltype
  res <- aggregate(
    x = df[[gene_name]], 
    by = list(celltype = df$celltype), 
    FUN = sum
  )
  colnames(res) <- c("celltype", "total_counts")
  
  # Calculate fraction independently
  res$fraction <- res$total_counts / sum(res$total_counts) * 100
  
  # Add gene label so we know which rows correspond to which gene
  res$gene <- gene_name
  
  return(res)
}

# -----------------------------
# 3️⃣ Run Calculation for Both Genes
# -----------------------------
# Run the function for Hspa5
stats_hspa5 <- calculate_stats(data_table, "Hspa5")

# Run the function for Ryr1
stats_ryr1  <- calculate_stats(data_table, "Ryr1")

# Combine them into one big table
final_results <- rbind(stats_hspa5, stats_ryr1)

# -----------------------------
# 4️⃣ View Results
# -----------------------------
print(final_results)

# -----------------------------
# 5️⃣ Visualization (Optional)
# -----------------------------
# This plot separates the two genes into side-by-side panels
ggplot(final_results, aes(x = celltype, y = fraction, fill = celltype)) +
  geom_col() +
  facet_wrap(~gene, scales = "free_y") + # Creates separate panels per gene
  labs(y = "% of Total Expression", title = "Expression Contribution by Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))