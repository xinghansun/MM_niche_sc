suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(stringr)
  library(data.table)
  library(harmony)
  library(ggplot2)
  library(ggrepel)
})


dir.create("../results/", recursive = TRUE, showWarnings = FALSE)
myeloid_labels <- c("Monocyte", "Macrophage", "DC", "Neutrophils")


# Define gene sets for state scoring
AP_GENES_REF <- c(
  "CIITA",
  "CD74",
  "HLA-DRA","HLA-DRB1","HLA-DRB5",
  "HLA-DPA1","HLA-DPB1",
  "HLA-DQA1","HLA-DQB1",
  "HLA-DMA","HLA-DMB",
  "HLA-DOA","HLA-DOB",
  "CTSS"
)

SUP_GENES_REF <- c(
  # checkpoint ligands / inhibitory receptors (myeloid side)
  "CD274","PDCD1LG2","VSIR","LILRB1","LILRB2",
  
  # suppressive cytokines / immune-metabolic axes
  "IL10","TGFB1","IDO1","ARG1","NOS2",
  
  # negative regulators / tolerance-associated nodes
  "IL1RN","SOCS3","HMOX1","LGALS9",
  
  # TAM/efferoctyosis-associated immunoregulatory receptors
  "MERTK","AXL"
)


# Load discovery scRNA data
obj_list <- readRDS("../data/seurat/GSE124310.seurat.list.rds")

# QC
obj_list <- lapply(obj_list, function(x) {
  # MT percentage
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  
  # remove low quality cells
  x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
  
  # remove low count genes
  counts <- GetAssayData(x, slot = "counts")
  keep_genes <- rowSums(counts > 0) >= 3
  x <- subset(x, features = keep_genes)
  
  return(x)
})

# Integrate samples using Harmony
combined <- merge(obj_list[[1]], 
                  y = obj_list[-1], 
                  add.cell.ids = names(obj_list))
rm(obj_list)
gc()

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)
ElbowPlot(combined, ndims = 50)

combined <- RunHarmony(combined, group.by.vars = "orig.ident")

# Clustering
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:50)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:50)
combined <- FindClusters(combined, resolution = 0.3)

DimPlot(combined, group.by = "orig.ident") # to check if batch effect is removed
UMAPPlot(combined, label = TRUE)


# Extract myeloid cells
library(SingleR)
library(celldex)
ref <- HumanPrimaryCellAtlasData() 
pred <- SingleR(test = as.SingleCellExperiment(combined), 
                ref = ref, labels = ref$label.main)
combined$SingleR_labels <- pred$labels

DimPlot(combined, group.by = "SingleR_labels") # based on SingleR annotation
tmp <- combined@meta.data
table(tmp[tmp$SingleR_labels %in% myeloid_labels, "seurat_clusters"], 
      tmp[tmp$SingleR_labels %in% myeloid_labels,"SingleR_labels"])


myeloid_markers <- c("CD68", "CD14", "FCGR3A", 
                     "LYZ", "S100A8", "S100A9", 
                     "CSF1R", "ITGAX")
combined <- AddModuleScore(
  object = combined,
  features = list(myeloid_markers),
  name = "Myeloid_Score"
) # based on Myeloid Score 
VlnPlot(combined, features = "Myeloid_Score1", group.by = "seurat_clusters")
FeaturePlot(combined, features = "Myeloid_Score1")

myeloid_obj <- subset(combined, idents = c(2, 9, 11))
myeloid_obj$disease_stage <- case_when(
  grepl("NBM", myeloid_obj$orig.ident) ~ "Normal",
  grepl("MGUS", myeloid_obj$orig.ident) ~ "MGUS",
  grepl("SMM", myeloid_obj$orig.ident) ~ "SMM",
  grepl("MM-", myeloid_obj$orig.ident) ~ "MM",
  TRUE ~ "Other"
)
myeloid_obj$disease_stage <- factor(myeloid_obj$disease_stage, 
                                    levels = c("Normal", "MGUS", "SMM", "MM"))
print(myeloid_obj)

# Recluster meyloid cells
myeloid_obj <- FindVariableFeatures(myeloid_obj, nfeatures = 2000)
myeloid_obj <- ScaleData(myeloid_obj)

myeloid_obj <- RunPCA(myeloid_obj, npcs = 30)
ElbowPlot(myeloid_obj, ndims = 30)
myeloid_obj <- RunHarmony(myeloid_obj, group.by.vars = "orig.ident")

myeloid_obj <- RunUMAP(myeloid_obj, reduction = "harmony", dims = 1:30)
myeloid_obj <- FindNeighbors(myeloid_obj, reduction = "harmony", dims = 1:30)
myeloid_obj <- FindClusters(myeloid_obj, resolution = 0.2)

DimPlot(myeloid_obj, group.by = "orig.ident")
DimPlot(myeloid_obj, group.by = "disease_stage")

UMAPPlot(myeloid_obj, label = TRUE)

# Check cluster identity
lineage_genes <- final_features <- c(
  # --- (Classical Monocytes: C0, C1, C4, C7) ---
  "CD14", "LYZ", "S100A8", "VCAN",
  
  # --- (IFN-Response) ---
  "ISG15", "IFITM3",
  
  # --- (Non-classical/Pathological) ---
  "FCGR3A", "MS4A7",
  "CD68", "C1QA",
  
  # --- (Dendritic Cells) ---
  "CD1C", "CLEC10A",
  "CLEC9A",
  "LILRA4"         
)
DotPlot(myeloid_obj, features = lineage_genes, group.by = "seurat_clusters") + coord_flip()

ref_data <- BlueprintEncodeData()
pred_results <- SingleR(test = myeloid_obj@assays$RNA@data, 
                        ref = ref_data, 
                        labels = ref_data$label.fine)
myeloid_obj$SingleR_fine <- pred_results$labels
table(myeloid_obj$seurat_clusters, myeloid_obj$SingleR_fine)


# Calculate module scores for antigen presentation and immunosuppressive signatures
AP_GENES  <- intersect(AP_GENES_REF, rownames(combined))
SUP_GENES <- intersect(SUP_GENES_REF, rownames(combined))

message(sprintf("Using %d antigen presentation genes: %s;\n\nDropped %d genes not found in dataset.", 
                length(AP_GENES), 
                paste(AP_GENES, collapse = ", "), 
                length(AP_GENES_REF) - length(AP_GENES)
)
)
message(sprintf("Using %d immunosuppressive genes: %s;\n\nDropped %d genes not found in dataset.", 
                length(SUP_GENES), 
                paste(SUP_GENES, collapse = ", "), 
                length(SUP_GENES_REF) - length(SUP_GENES)
)
)

define_myeloid_states <- function(obj, ap_genes, sup_genes) {
  obj <- AddModuleScore(obj, features = list(ap_genes), name = "AP_Score")
  obj <- AddModuleScore(obj, features = list(sup_genes), name = "SUP_Score")
  obj$AP_val <- obj$AP_Score1
  obj$SUP_val <- obj$SUP_Score1
  
  ap_thr  <- median(obj$AP_val, na.rm = TRUE)
  sup_thr <- median(obj$SUP_val, na.rm = TRUE)
  
  obj$functional_state <- "Intermediate"
  obj$functional_state[obj$AP_val <= ap_thr & obj$SUP_val >= sup_thr] <- "AP_low_Sup_high"
  obj$functional_state[obj$AP_val > ap_thr & obj$SUP_val < sup_thr] <- "AP_high_Sup_low"
  obj$functional_state[obj$AP_val <= ap_thr & obj$SUP_val < sup_thr] <- "AP_low_Sup_low"
  obj$functional_state[obj$AP_val > ap_thr & obj$SUP_val >= sup_thr] <- "AP_high_Sup_high"
  
  obj$functional_state <- factor(obj$functional_state, 
                                     levels = c("AP_high_Sup_high",
                                                "AP_high_Sup_low", 
                                                "Intermediate",
                                                "AP_low_Sup_high",
                                                "AP_low_Sup_low"))
  
  return(obj)
}

myeloid_obj <- define_myeloid_states(myeloid_obj, 
                                     ap_genes = AP_GENES, 
                                     sup_genes = SUP_GENES)

FeaturePlot(myeloid_obj, 
            features = c("AP_val", "SUP_val"), 
            ncol = 2, cols = c("lightgrey", "red"))

VlnPlot(myeloid_obj, features = c("AP_val", "SUP_val"), pt.size = 0)

DimPlot(myeloid_obj, 
        group.by = "functional_state")

DimPlot(myeloid_obj, 
        group.by = "functional_state", 
        split.by = "disease_stage")


# Check AP high SUP high distribution across disease stages
# A. pooled distribution
prop_data <- myeloid_obj@meta.data %>%
  group_by(disease_stage, functional_state) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

prop_data$disease_stage <- factor(prop_data$disease_stage, 
                                  levels = c("Normal", "MGUS", "SMM", "MM"))

ggplot(prop_data, aes(x = disease_stage, y = percentage, fill = functional_state)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "white") +
  scale_fill_manual(values = c("AP_high_Sup_high" = "#E41A1C",
                               "AP_high_Sup_low"  = "#4DAF4A",
                               "AP_low_Sup_high"  = "#377EB8",
                               "AP_low_Sup_low"   = "#984EA3")) +
  theme_classic() +
  labs(title = "Proportion of Functional States across Stages",
       x = "Disease Stage", y = "Proportion (%)", fill = "Functional State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# B. Proportion test by sample
sample_props <- myeloid_obj@meta.data %>%
  group_by(orig.ident, disease_stage, functional_state) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

library(ggpubr)

ggplot(filter(sample_props, functional_state == "AP_high_Sup_high"), 
       aes(x = disease_stage, y = prop, fill = disease_stage)) +
  geom_boxplot(outlier.shape = NA)  + 
  geom_jitter(width = 0.2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Normal", "MM"), 
                                        c("Normal", "SMM"),
                                        c("Normal", "MGUS")
                                        )) +
  theme_classic() +
  labs(title = "Proportion of AP_high-SUP_high per Sample", 
       y = "Proportion")

FeaturePlot(myeloid_obj, 
            features = c("AP_val", "SUP_val"), 
            blend = T)

# C. Compare expansion magnitude (Fold Change) between Red and Green states across disease stages
sample_props <- myeloid_obj@meta.data %>%
  group_by(orig.ident, disease_stage, functional_state) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

normal_baseline <- sample_props %>%
  filter(disease_stage == "Normal") %>%
  group_by(functional_state) %>%
  summarise(baseline_prop = mean(prop))

fc_data <- sample_props %>%
  left_join(normal_baseline, by = "functional_state") %>%
  mutate(fold_change = prop / baseline_prop) %>%
  filter(functional_state %in% c("AP_high_Sup_high", "AP_high_Sup_low"))
fc_data$disease_stage <- factor(fc_data$disease_stage, levels = c("Normal", "MGUS", "SMM", "MM"))

ggplot(fc_data, aes(x = disease_stage, y = fold_change, fill = functional_state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 1, alpha = 0.5) +
  stat_compare_means(aes(group = functional_state), 
                     label = "p.signif", 
                     method = "wilcox.test",
                     hide.ns = FALSE) +
  scale_fill_manual(values = c("AP_high_Sup_high" = "#E41A1C", "AP_high_Sup_low" = "#4DAF4A")) +
  theme_classic() +
  labs(title = "Expansion Magnitude Comparison: Red vs Green",
       y = "Fold Change (Relative to Normal Mean)",
       x = "Disease Stage") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")

# D. Check myeloid cell clusters' AP/SUP states distribution
cluster_dist <- myeloid_obj@meta.data %>%
  group_by(seurat_clusters, functional_state) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

cluster_order <- cluster_dist %>%
  filter(functional_state == "AP_high_Sup_high") %>%
  arrange(desc(proportion)) %>%
  pull(seurat_clusters)

cluster_dist$seurat_clusters <- factor(cluster_dist$seurat_clusters, levels = cluster_order)

ggplot(cluster_dist, aes(x = seurat_clusters, y = proportion, fill = functional_state)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, color = "white") +
  scale_fill_manual(values = c("AP_high_Sup_high" = "#E41A1C",
                               "AP_high_Sup_low"  = "#4DAF4A",
                               "AP_low_Sup_high"  = "#377EB8",
                               "AP_low_Sup_low"   = "#984EA3",
                               "Intermediate"     = "#D3D3D3")) +
  theme_classic() +
  labs(title = "Distribution of Functional States across Myeloid Clusters",
       x = "Seurat Cluster (Identity)", 
       y = "Proportion", 
       fill = "Functional State") +
  theme(axis.text.x = element_text(face = "bold"))


# Extract Ap_high_Sup_high cells markers for downstream analysis
# drop Cluster 6, 3 as they are likely DC, distinct from monocytes, and 6 has small cell num
# drop Cluster 7 as it's very distinct from the monocyte clusters and small cell num
# drop Cluster 4 as it may 
intrested_subset <- subset(myeloid_obj, idents = c("0", "1", "2", "4", "5"))
Idents(intrested_subset) <- "functional_state"

# DEG volcano
deg_ApSupHigh_posneg <- FindMarkers(intrested_subset, 
                             ident.1 = "AP_high_Sup_high", 
                             ident.2 = NULL,
                             test.use = "LR",
                             latent.vars = c("seurat_clusters"),
                             logfc.threshold = 0,
                             min.pct = 0.25)

plot_data <- deg_ApSupHigh_posneg
plot_data$gene_symbol <- rownames(plot_data)

pval_cutoff <- 0.05
logfc_cutoff <- 0.25 

plot_data <- plot_data %>%
  mutate(
    minus_log10_pval = -log10(p_val_adj + 1e-300),
    category = case_when(
      gene_symbol %in% AP_GENES ~ "AP Markers",
      gene_symbol %in% SUP_GENES ~ "SUP Markers",
      p_val_adj < pval_cutoff & avg_log2FC > logfc_cutoff ~ "Upregulated (Novel)",
      p_val_adj < pval_cutoff & avg_log2FC < -logfc_cutoff ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

my_colors <- c(
  "AP Markers" = "#FFD700",
  "SUP Markers" = "#8A2BE2",
  "Upregulated (Novel)" = "red3",
  "Downregulated" = "steelblue",
  "Not Significant" = "grey90"
)

top_novel_up <- plot_data %>% 
  filter(category == "Upregulated (Novel)") %>% 
  top_n(10, wt = avg_log2FC) %>% 
  pull(gene_symbol)
genes_to_label <- unique(c(AP_GENES, SUP_GENES, top_novel_up))
label_data <- plot_data %>% 
  filter(gene_symbol %in% genes_to_label & category != "Not Significant")

ggplot(plot_data, aes(x = avg_log2FC, y = minus_log10_pval, color = category)) +
  geom_point(data = subset(plot_data, category == "Not Significant"), 
             size = 1, alpha = 0.3) +
  geom_point(data = subset(plot_data, category != "Not Significant"), 
             size = 2.5, alpha = 0.9) +
  scale_color_manual(values = my_colors) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), lty = 2, col = "grey50", lwd = 0.5) +
  geom_hline(yintercept = -log10(pval_cutoff), lty = 2, col = "grey50", lwd = 0.5) +
  geom_text_repel(data = label_data,
                  aes(label = gene_symbol),
                  size = 3.5,
                  box.padding = 0.5,
                  max.overlaps = 30,
                  fontface = "bold",
                  show.legend = FALSE) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    axis.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "log2(Fold Change)", 
    y = "-log10(Adjusted P-value)",
    title = "Transcriptomic Signature: AP (Gold) vs SUP (Purple)",
    subtitle = "Distinguishing defined markers form novel upregulated genes"
  )


# Extract novel markers for bulk scoring
deg_ApSupHigh_pos <- FindMarkers(intrested_subset, 
                                 ident.1 = "AP_high_Sup_high", 
                                 ident.2 = NULL,
                                 test.use = "LR",
                                 only.pos = TRUE,
                                 latent.vars = c("seurat_clusters"),
                                 logfc.threshold = 0.25,
                                 min.pct = 0.25)

print(nrow(deg_ApSupHigh_pos))
deg_ApSupHigh_pos_novel <- deg_ApSupHigh_pos[!rownames(deg_ApSupHigh_pos) %in% c(AP_GENES, SUP_GENES), ]
print(nrow(deg_ApSupHigh_pos_novel))

write.table(deg_ApSupHigh_pos, 
            file = "../results/myeloid_APSupHigh_pos_markers.tsv", 
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = NA)

write.table(deg_ApSupHigh_pos_novel, 
            file = "../results/myeloid_APSupHigh_novel_markers.tsv", 
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = NA)
