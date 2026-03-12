# ============================================================
# 0) Load Libraries
# ============================================================
library(ggplot2)
library(DESeq2)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(grid)

# ============================================================
# 1-3) Data Loading & Alignment
# ============================================================
expression_data <- read.csv("expression_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)
meta_data <- read.csv("metadata.csv", header = TRUE, stringsAsFactors = FALSE)

# Robust metadata formatting
meta_formatted <- meta_data
rownames(meta_formatted) <- meta_formatted$Sample_id 
meta_formatted <- meta_formatted %>%
  mutate(
    subject_id = factor(subject_id),
    gender     = factor(gender),
    timepoint  = factor(timepoint, levels = c("baseline", "post intervention"))
  )

meta_formatted <- meta_formatted[colnames(expression_data), ]

# Filter low-count technical noise
keep <- rowSums(expression_data >= 10) >= 16
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(expression_data[keep,])),
                              colData = meta_formatted,
                              design = ~ subject_id + timepoint)

# ============================================================
# 5) Run DESeq2 & Map Gene Symbols
# ============================================================
dds <- DESeq(dds)
res <- results(dds, contrast = c("timepoint", "post intervention", "baseline"))

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  mutate(symbol = mapIds(org.Hs.eg.db, keys = gene_id, column = "SYMBOL", 
                         keytype = "ENSEMBL", multiVals = "first")) %>%
  filter(!is.na(symbol))

# ============================================================
# 6) MA Plot
# ============================================================
# Shows the relationship between mean expression and log fold change
plotMA(res, ylim = c(-3, 3), main = "MA Plot: Post vs Baseline")

# ============================================================
# 7) Volcano Plot (Significance Colored)
# ============================================================
res_df <- res_df %>%
  mutate(diffexpressed = case_when(
    log2FoldChange > 0.5 & pvalue < 0.05 ~ "Upregulated",
    log2FoldChange < -0.5 & pvalue < 0.05 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = diffexpressed)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "#de2d26", "Downregulated" = "#3182bd", "Not Significant" = "grey")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ylim(0, 5.8) + theme_classic() + labs(title = "Volcano Plot: Dietary Impact")

#interactive/ function for lable next to point

# ============================================================
# 8-9) Functional Enrichment (GO) & Dotplot
# ============================================================
sig_genes <- res_df %>% filter(pvalue < 0.05) %>% pull(gene_id)

ego_bp <- enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                   ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, readable = TRUE, minGSSize = 5,
                   maxGSSize = 500)

# Enrichment Visualization
dotplot(ego_bp, showCategory = 10) + ggtitle("GO Biological Process Enrichment")
barplot(ego_bp, showCategory = 10)

#more strict cutoff

# ============================================================
# 10) Heatmap Data Preparation
# ============================================================
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

# Create 16-person "Net Change" matrix
subjects <- levels(meta_formatted$subject_id)
lfc_matrix <- matrix(NA, nrow = nrow(vst_mat), ncol = length(subjects))
colnames(lfc_matrix) <- subjects
rownames(lfc_matrix) <- rownames(vst_mat)

for (i in seq_along(subjects)) {
  sub <- subjects[i]
  idx_post <- which(meta_formatted$subject_id == sub & meta_formatted$timepoint == "post intervention")
  idx_base <- which(meta_formatted$subject_id == sub & meta_formatted$timepoint == "baseline")
  lfc_matrix[,i] <- vst_mat[,idx_post] - vst_mat[,idx_base]
}

ann_colors = list(
  gender = c(F = "#f1a340", M = "#998ec3"),
  timepoint = c(baseline = "#7fbf7b", `post intervention` = "#af8dc3")
)

lfc_threshold <- 0.5
p_threshold <- 0.05

top_genes <- res_df %>% arrange(pvalue) %>% filter(abs(log2FoldChange) > lfc_threshold & pvalue < p_threshold) %>%
  arrange(pvalue)

# ============================================================
# 11) FINAL HEATMAPS
# ============================================================

# 1. Define your scale variables BEFORE the function call
limit <- max(abs(plot_16), na.rm = TRUE) 
my_breaks <- seq(-limit, limit, length.out = 101)
my_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# 2. Run the pheatmap function
pheatmap(plot_16, 
         main = "Heatmap 1: Net Change per Individual (16 Columns)",
         annotation_col = ann_16, 
         annotation_colors = ann_colors,
         cluster_cols = TRUE, 
         border_color = "white",
         
         # Use the variables we just created:
         color = my_colors,
         breaks = my_breaks
) # <--- This closing bracket is essential!

# 3. Add the legend label using the grid package
library(grid)
grid.text("log2 Fold Change", x=0.95, y=0.5, rot=270, gp=gpar(fontsize=10, fontface="bold"))
