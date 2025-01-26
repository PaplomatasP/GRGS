# Φόρτωση του αρχείου LR pairs
lr_data <- read.delim("human_lr_pair.txt", header = TRUE, sep = "\t")

# Μετατροπή σε κεφαλαία και αφαίρεση κενών
lr_genes_clean <- unique(toupper(trimws(c(lr_data$ligand_gene_symbol, lr_data$receptor_gene_symbol))))

# Δημιουργία λίστας με τα subsets και καθαρισμός τους
subsets <- list(
  LASSO = unique(toupper(trimws(results_GA_LASSO_Boruta$LASSO$genes))),
  GA = unique(toupper(trimws(results_GA_LASSO_Boruta$GA$genes))),
  Boruta = unique(toupper(trimws(results_GA_LASSO_Boruta$Boruta$genes))),
  GRGS = unique(toupper(trimws(results124.20.01.25$ga_genes)))
  #FilterData = unique(toupper(trimws(colnames(filtered_data)))),
  #FullDataset = unique(toupper(trimws(colnames(thca_gse_154763))))
)

# Υπολογισμός LR γονιδίων για κάθε subset
subset_counts <- sapply(subsets, function(genes) {
  length(intersect(genes, lr_genes_clean))
})

# Υπολογισμός συνολικού πλήθους γονιδίων σε κάθε subset
subset_sizes <- sapply(subsets, length)

# Υπολογισμός ποσοστού LR γονιδίων
lr_ratios <- subset_counts / subset_sizes * 100  # Ποσοστό (%)

# Δημιουργία dataframe για το πλοτ
lr_df <- data.frame(
  Subset = names(subsets),
  LR_Genes = subset_counts,
  Total_Genes = subset_sizes,
  LR_Ratio = lr_ratios
)

# Εμφάνιση του dataframe
print(lr_df)

#lr_df$LR_Ratio[1] = 19.6
# Δημιουργία πλοτ
# Δημιουργία πλοτ για paper
library(ggplot2)

ggplot(lr_df, aes(x = Subset, y = LR_Ratio, fill = Subset)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%%", LR_Ratio)),
    vjust = 1.5, size = 6, color = "white", fontface = "bold",
    position = position_dodge(width = 0.7)
  ) +
  scale_fill_brewer(palette = "Dark2") +  # Χρησιμοποιούμε μια επαγγελματική παλέτα
  labs(
    title = "Ligand-Receptor Gene Enrichment: GRGS vs Traditional Feature Selection Methods",
    #subtitle = "Evaluating the effectiveness of GRGS against Boruta, RFE, and GA in identifying ligand-receptor genes.",
    x = "Feature Selection Method",
    y = "Percentage of LR Genes (%)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.5, color = "grey85"),
    panel.grid.minor = element_blank()
  )

ggsave("lr_comparison.png", dpi = 300, width = 10, height = 8, bg = "white")






# Φόρτωση του αρχείου LR pairs
lr_data <- read.delim("human_lr_pair.txt", header = TRUE, sep = "\t")

# Δημιουργία λίστας με τα subsets
subsets <- list(
  RFE    = results_GA_RFE_Boruta$RFE$genes,
  GA     = results_GA_RFE_Boruta$GA$genes,
  Boruta = results_GA_RFE_Boruta$Boruta$genes,
  GRGS   = results$ga_genes#,
  #FullDataset = colnames(filtered_data)
)

# Υπολογισμός LR γονιδίων (απόλυτο πλήθος) για κάθε subset
subset_counts <- sapply(subsets, function(genes) {
  length(intersect(genes, unique(c(lr_data$ligand_gene_symbol, lr_data$receptor_gene_symbol))))
})

# Υπολογισμός συνολικού πλήθους γονιδίων σε κάθε subset (για επιπλέον πληροφορία αν χρειαστεί)
subset_sizes <- sapply(subsets, length)

# Δημιουργία dataframe για το πλοτ
lr_df <- data.frame(
  Subset      = names(subsets),
  LR_Genes    = subset_counts,
  Total_Genes = subset_sizes

)

# Δημιουργία πλοτ με απόλυτα νούμερα LR γονιδίων
library(ggplot2)

ggplot(lr_df, aes(x = Subset, y = LR_Genes, fill = Subset)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_text(
    aes(label = LR_Genes), 
    vjust = 1.5, size = 6, color = "white", fontface = "bold",
    position = position_dodge(width = 0.7)
  ) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Ligand-Receptor Gene Enrichment: GRGS vs Traditional Feature Selection Methods",
    x = "Feature Selection Method",
    y = "Number of LR Genes Found"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5, size = 14, face = "italic"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    axis.text.y      = element_text(size = 14),
    axis.title       = element_text(size = 16, face = "bold"),
    legend.position  = "none",
    panel.grid.major = element_line(size = 0.5, color = "grey85"),
    panel.grid.minor = element_blank()
  )

# Αποθήκευση πλοτ
ggsave("lr_comparison.png", dpi = 300, width = 10, height = 8, bg = "white")



filtered_data

##### ##### ##### ##### ##### Non-lr ##### ##### ##### ##### ##### ##### ##### 

# υργία λίστας με τα subsets
subsets <- list(
  RFE    = results_GA_RFE_Boruta$RFE$genes,
  GA     = results_GA_RFE_Boruta$GA$genes,
  Boruta = results_GA_RFE_Boruta$Boruta$genes,
  GRGS   = results$ga_genes
)

# Υπολογισμός Non-LR γονιδίων
non_lr_counts <- sapply(subsets, function(genes) {
  length(setdiff(genes, unique(c(lr_data$ligand_gene_symbol, lr_data$receptor_gene_symbol))))
})

# Υπολογισμός ποσοστού Non-LR γονιδίων
non_lr_ratios <- non_lr_counts / subset_sizes * 100  # Ποσοστό (%)

# Δημιουργία dataframe για το πλοτ
non_lr_df <- data.frame(
  Subset = names(subsets),
  Non_LR_Genes = non_lr_counts,
  Total_Genes = subset_sizes,
  Non_LR_Ratio = non_lr_ratios
)

# Δημιουργία πλοτ για το ποσοστό των Non-LR γονιδίων
library(ggplot2)

ggplot(non_lr_df, aes(x = Subset, y = Non_LR_Ratio, fill = Subset)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%%", Non_LR_Ratio)),
    vjust = 1.5, size = 6, color = "white", fontface = "bold",
    position = position_dodge(width = 0.7)
  ) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Percentage of Non-LR Genes Across Feature Selection Methods",
    x = "Feature Selection Method",
    y = "Percentage of Non-LR Genes (%)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    axis.text.y      = element_text(size = 14),
    axis.title       = element_text(size = 16, face = "bold"),
    legend.position  = "none",
    panel.grid.major = element_line(size = 0.5, color = "grey85"),
    panel.grid.minor = element_blank()
  )

# Αποθήκευση πλοτ
#ggsave("non_lr_comparison.png", dpi = 300, width = 10, height = 8, bg = "white")


##########################    #########################


# Φόρτωση αρχείου με τα LR pairs
lr_data <- read.delim("human_lr_pair.txt", header = TRUE, sep = "\t")

# Δημιουργία λίστας με τα subsets για κάθε αλγόριθμο
subsets <- list(
  RFE = results_GA_RFE_Boruta$RFE$genes,
  GA = results_GA_RFE_Boruta$GA$genes,
  Boruta = results_GA_RFE_Boruta$Boruta$genes,
  GRGS = results$ga_genes,
  FilterData = colnames(filtered_data),
  FullDataset = colnames(thca_gse_154763)
)

# Υπολογισμός συνολικών LR γονιδίων από τη βάση δεδομένων
total_lr_genes <- unique(c(lr_data$ligand_gene_symbol, lr_data$receptor_gene_symbol))
total_lr_count <- length(total_lr_genes)

# Υπολογισμός γονιδίων και LR γονιδίων για κάθε subset
algorithm_results <- lapply(names(subsets), function(algorithm) {
  genes <- subsets[[algorithm]]
  total_genes <- length(genes)
  lr_genes <- length(intersect(genes, total_lr_genes))
  list(
    Algorithm = algorithm,
    Total_Genes = total_genes,
    LR_Genes = lr_genes,
    LR_Ratio = (lr_genes / total_genes) * 100
  )
})

# Μετατροπή των αποτελεσμάτων σε dataframe
algorithm_df <- do.call(rbind, lapply(algorithm_results, as.data.frame))
rownames(algorithm_df) <- NULL

# Εμφάνιση αποτελεσμάτων
print(algorithm_df)

# Δημιουργία πλοτ για το paper
library(ggplot2)

ggplot(algorithm_df, aes(x = Algorithm, y = LR_Ratio, fill = Algorithm)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%%", LR_Ratio)),
    vjust = 1.5, size = 6, color = "white", fontface = "bold",
    position = position_dodge(width = 0.7)
  ) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Ligand-Receptor Gene Enrichment by Algorithm",
    x = "Algorithm",
    y = "Percentage of LR Genes (%)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.5, color = "grey85"),
    panel.grid.minor = element_blank()
  )
