library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
#library(umap)
library(cluster)
library(uwot)
#library(factoextra)  # for gap statistic plotting (optional)

df <- read_excel("File_Name.xlsx")

df <- df %>%
  mutate(
    behavior_group = case_when(
      sign(NSF) == 1 & sign(FST) == 1 ~ "both_pos",
      sign(NSF) == -1 & sign(FST) == -1 ~ "both_neg",
      TRUE ~ "opposite"
    ),
    NSF_abs = abs(NSF),
    FST_abs = abs(FST),
    NSF_sq = NSF^2,
    FST_sq = FST^2,
    NSF_FST_product = NSF * FST,
    abs_diff = abs(NSF - FST),
    same_sign = sign(NSF) == sign(FST),
    same_sign_num = as.integer(same_sign),
    has_zero_flag = ifelse(NSF == 0 | FST == 0, 1, 0)
  )

continuous_features <- c("NSF", "FST", "NSF_abs", "FST_abs", 
                         "NSF_sq", "FST_sq", "NSF_FST_product", "abs_diff")
df_scaled <- df %>%
  mutate(across(all_of(continuous_features), scale))

features_scaled <- df_scaled %>%
  select(NSF, FST, NSF_FST_product, same_sign_num)

# --- UMAP projection ---
set.seed(42)
umap_result <- umap(features_scaled, n_neighbors = 4, min_dist = 0.3)

cluster_df <- as.data.frame(umap_result)
colnames(cluster_df) <- c("PC1", "PC2")
cluster_df$behavior_group <- df$behavior_group

# --- ADDED: K-means clustering and statistical evaluation ---
set.seed(42)
kmeans_result <- kmeans(cluster_df[, c("PC1", "PC2")], centers = 3)
cluster_df$kmeans_cluster <- as.factor(kmeans_result$cluster)

# Silhouette score
sil <- silhouette(kmeans_result$cluster, dist(cluster_df[, c("PC1", "PC2")]))
cat("Average silhouette width:", mean(sil[, 3]), "\n")

# Chi-square test vs. known behavior groups
cat("Chi-square test:\n")
chisq_result <- chisq.test(table(cluster_df$kmeans_cluster, cluster_df$behavior_group))
print(chisq_result)


# --- END ADDED ---

# UMAP group label centers
group_labels <- cluster_df %>%
  group_by(behavior_group) %>%
  summarise(
    n = n(),
    center_x = median(PC1),
    center_y = median(PC2)
  ) %>%
  mutate(label = paste0(behavior_group, "\n(n=", n, ")"))

# --- UMAP PLOT by behavior_group ---
p <- ggplot(cluster_df, aes(x = PC1, y = PC2, color = behavior_group)) +
  stat_ellipse(type = "norm", level = 0.95, aes(fill = behavior_group), 
               geom = "polygon", alpha = 0.15, color = "black", size = 0.3) +  # <-- updated
  geom_point(
    shape = 21, size = 2, stroke = 1, fill = "white",
    aes(color = behavior_group), alpha = 0.9
  ) +
  geom_text_repel(data = group_labels,
                  aes(x = center_x, y = center_y, label = label),
                  color = "black", size = 2, fontface = "bold",
                  box.padding = 0.6, max.overlaps = Inf, inherit.aes = FALSE) +
  scale_color_manual(values = c("both_neg" = "blue", "both_pos" = "red", "opposite" = "gray50")) +
  scale_fill_manual(values = c("both_neg" = "blue", "both_pos" = "red", "opposite" = "gray50")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "UMAP with behavior group and sample counts",
       x = "UMAP-1", y = "UMAP-2")

print(p)

# --- OPTIONAL: Plot clusters found by k-means ---
p2 <- ggplot(cluster_df, aes(x = PC1, y = PC2, color = kmeans_cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(title = "K-means Clustering (k=3) on UMAP", x = "UMAP-1", y = "UMAP-2") +
  theme_minimal()
print(p2)
