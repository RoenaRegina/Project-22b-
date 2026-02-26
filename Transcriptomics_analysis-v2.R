# ============================================================
# RNA-seq Differential Expression Analysis
# Dietary Intervention Study
# Design: Paired (Repeated Measures) within subjects
# ============================================================


# ============================================================
# 0) Install and load required packages
# ============================================================

# Biological meaning:
# These packages allow us to:
# - Analyze count-based RNA-seq data (DESeq2)
# - Manipulate data (dplyr, tibble)
# - Visualize results (ggplot2)
# - Perform functional enrichment analysis (clusterProfiler)

pkgs_cran <- c("ggplot2", "dplyr", "tibble", "readr")
pkgs_bioc <- c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "enrichplot")

to_install_cran <- pkgs_cran[!pkgs_cran %in% rownames(installed.packages())]
if (length(to_install_cran) > 0) install.packages(to_install_cran)

if (!"BiocManager" %in% rownames(installed.packages())) install.packages("BiocManager")

to_install_bioc <- pkgs_bioc[!pkgs_bioc %in% rownames(installed.packages())]
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, ask = FALSE)

library(ggplot2)
library(DESeq2)
library(dplyr)
library(tibble)
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# ============================================================
# 1) Load expression matrix and metadata
# ============================================================

# Biological meaning:
# expression_data:
#   - rows = genes
#   - columns = samples
#   - values = raw RNA-seq counts (number of reads mapped to each gene)
#
# meta_data:
#   - contains information about each sample
#   - includes subject_id (individual)
#   - includes timepoint (baseline vs post intervention)


expression_data <- read.csv("expression_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
meta_data <- read.csv("metadata.csv", header = TRUE, stringsAsFactors = FALSE)
# ============================================================
# 2) Format metadata and define biological model
# ============================================================

# Biological meaning:
# We convert subject_id and timepoint into factors because:
# - subject_id accounts for individual-specific baseline expression differences
# - timepoint models the biological effect of the dietary intervention
#
# By setting baseline as reference level,
# the model estimates how gene expression changes AFTER intervention.

meta_formatted <- meta_data %>%
  mutate(
    subject_id = factor(subject_id),
    timepoint  = factor(timepoint,
                        levels = c("baseline", "post intervention"))
  ) %>%
  column_to_rownames("Sample_id")

# Align metadata to expression matrix
meta_formatted <- meta_formatted[colnames(expression_data), ]

stopifnot(all(rownames(meta_formatted) == colnames(expression_data)))

# ============================================================
# 3) Gene filtering
# ============================================================

# Biological meaning:
# Genes with extremely low counts across samples
# are often:
# - Not biologically active
# - Technical noise
#
# Removing them improves statistical power
# and reduces multiple testing burden.

keep <- rowSums(expression_data >= 10) >= 16
expression_data_filtered <- expression_data[keep, ]

# ============================================================
# 4) Construct DESeq2 model (paired design)
# ============================================================

# Statistical model:
# design = ~ subject_id + timepoint
#
# Biological meaning:
# This means:
# - First account for inter-individual variability
# - Then estimate the intervention effect
#
# So we are testing:
# "Within the same individual, did gene expression change after diet?"

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(expression_data_filtered)),
  colData   = meta_formatted,
  design    = ~ subject_id + timepoint
)

# ============================================================
# 5) Run differential expression analysis
# ============================================================

# Biological meaning:
# DESeq performs:
# 1) Normalization (adjust for sequencing depth)
# 2) Estimate dispersion (biological variability)
# 3) Fit negative binomial model
# 4) Hypothesis testing

dds <- DESeq(dds)

res <- results(dds,
               contrast = c("timepoint",
                            "post intervention",
                            "baseline"))
# Biological meaning:
# log2FoldChange:
#   > 0  → gene upregulated after intervention
#   < 0  → gene downregulated after intervention
#
# padj:
#   Adjusted p-value (multiple testing correction using BH method)
#   Controls false discovery rate (FDR)

# Which genes' adjust p-value are lower than 0.05? 

# ============================================================
# 6) MA plot
# ============================================================

# Biological meaning:
# MA plot shows:
# - Mean expression (A)
# - Log fold change (M)
#
# It helps identify:
# - Systematic bias
# - Strongly regulated genes

plotMA(res, ylim = c(-3, 3))

# ============================================================
# 7) Volcano plot
# ============================================================

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    neglog10_p = -log10(pvalue),
    sig = case_when(
      pvalue < 0.05 & log2FoldChange > 1  ~ "Up",
      pvalue < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Biological meaning:
# Volcano plot combines:
# - Effect size (log2FC)
# - Statistical significance (-log10 p-value)
#
# It helps identify biologically meaningful DE genes.

ggplot(res_df,
       aes(x = log2FoldChange,
           y = neglog10_p)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme_classic() +
  labs(
    x = "log2FC (post intervention vs baseline)",
    y = "-log10(p-value)"
  )

# ============================================================
# 8) Functional enrichment analysis (GO)
# ============================================================

# Biological meaning:
# Individual genes are hard to interpret.
# Therefore we ask:
#
# "Are specific biological processes collectively affected?"
#
# GO Biological Process (BP) tells us:
# - Immune response?
# - Lipid metabolism?
# - Inflammation?
# - Oxidative stress?
#
# This helps translate gene-level changes
# into pathway-level biological interpretation.

sig_genes <- res_df %>%
  filter(pvalue < 0.05) %>%
  pull(gene_id)


ego_bp <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

# Biological meaning:
# Enriched GO terms represent biological processes
# significantly overrepresented among DE genes.
#
# Example interpretation:
# If "inflammatory response" is enriched,
# the diet may modulate immune activation.


# ============================================================
# 9) Visualize GO results
# ============================================================

#Show first 10 pathways
#What are these pathways' function? How to intepret the intervention impact?

dotplot(ego_bp, showCategory = 10) +
  ggtitle("GO Biological Process Enrichment")

barplot(ego_bp, showCategory = 10)

dotplot(ego_bp, showCategory = 50) +
  ggtitle("GO Biological Process Enrichment")

barplot(ego_bp, showCategory = 50)

# ============================================================
# 10. Separate samples in males and females and run a PCA plot 
# ============================================================

#Run Differential Expression Analysis (Base Model)

# Paired analysis using DESeq2
# Baseline vs Post Intervention

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(expression_data_filtered)),
  colData   = meta_formatted,
  design    = ~ subject_id + timepoint
)

dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("timepoint", "post intervention", "baseline"))

summary(res)


#Transform Data for Visualization (PCA)

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA colored by intervention
plotPCA(vsd, intgroup = "timepoint") +
  theme_minimal() +
  ggtitle("PCA: Baseline vs Post Intervention")

# PCA colored by gender
plotPCA(vsd, intgroup = "gender") +
  theme_minimal() +
  ggtitle("PCA: Clustering by Gender")

#Test Whether Gender is a Confounder
#PERMANOVA test

if (!require("vegan")) install.packages("vegan")
library(vegan)

dist_matrix <- vegdist(t(assay(vsd)), method = "euclidean")

set.seed(123)
gender_test <- adonis2(dist_matrix ~ gender,
                       data = as.data.frame(colData(vsd)))

print(gender_test)

#Model Adjusting for Gender

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(expression_data_filtered)),
  colData   = meta_formatted,
  design    = ~ gender + subject_id + timepoint
)

dds <- DESeq(dds)

res_universal <- results(dds,
                         contrast=c("timepoint",
                                    "post intervention",
                                    "baseline"))

#Test Gender × Intervention Interaction

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(expression_data_filtered)),
  colData   = meta_formatted,
  design    = ~ gender + subject_id + timepoint + gender:timepoint
)

dds <- DESeq(dds)

res_interaction <- results(dds,
                           name="genderM.timepointpost.intervention")

#Female-Specific Analysis
#Subset data

female_samples <- rownames(meta_formatted)[meta_formatted$gender == "F"]

meta_female <- meta_formatted[female_samples, ]
expression_female <- expression_data_filtered[, female_samples]

stopifnot(all(rownames(meta_female) == colnames(expression_female)))

#Run paired DESeq2

dds_female <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(expression_female)),
  colData   = meta_female,
  design    = ~ subject_id + timepoint
)

dds_female <- DESeq(dds_female)

res_female <- results(dds_female,
                      contrast = c("timepoint",
                                   "post intervention",
                                   "baseline"))

summary(res_female)

#Male-Specific Analysis

male_samples <- rownames(meta_formatted)[meta_formatted$gender == "M"]

meta_male <- meta_formatted[male_samples, ]
expression_male <- expression_data_filtered[, male_samples]

dds_male <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(expression_male)),
  colData   = meta_male,
  design    = ~ subject_id + timepoint
)

dds_male <- DESeq(dds_male)

res_male <- results(dds_male,
                    contrast = c("timepoint",
                                 "post intervention",
                                 "baseline"))

summary(res_male)

#Identify Significant Female Genes

sig_female_genes <- subset(res_female, padj < 0.1)
print(sig_female_genes)

#PCA to Check Age Effect in Females

vsd_f <- vst(dds_female, blind = FALSE)

plotPCA(vsd_f, intgroup = "age") +
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("Female PCA: Age Effects")

#Likelihood Ratio Test (Does Gender Matter?)

dds_lrt <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(expression_data_filtered)),
  colData   = meta_formatted,
  design    = ~ gender + timepoint
)

dds_lrt <- DESeq(dds_lrt,
                 test = "LRT",
                 reduced = ~ timepoint)

res_lrt <- results(dds_lrt)
summary(res_lrt)

#Convert Ensembl IDs to Gene Names

if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

female_ens_ids <- rownames(sig_female_genes)

gene_map <- select(org.Hs.eg.db,
                   keys = female_ens_ids,
                   columns = c("SYMBOL", "GENENAME"),
                   keytype = "ENSEMBL")

final_female_table <- merge(as.data.frame(sig_female_genes),
                            gene_map,
                            by.x = 0,
                            by.y = "ENSEMBL")

print(final_female_table[, c("SYMBOL", "log2FoldChange", "padj", "GENENAME")])

#GO Enrichment Analysis

if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
library(clusterProfiler)

ego <- enrichGO(
  gene         = female_ens_ids,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENSEMBL",
  ont          = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

dotplot(ego, showCategory = 10) +
  ggtitle("Top Biological Pathways Affected by Diet in Females")