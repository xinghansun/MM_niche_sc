suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(stringr)
  library(data.table)
  library(harmony)
  library(ggplot2)
})

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
myeloid_obj <- FindClusters(myeloid_obj, resolution = 0.4)

DimPlot(myeloid_obj, group.by = "orig.ident")
DimPlot(myeloid_obj, group.by = "disease_stage")

UMAPPlot(myeloid_obj, label = TRUE)


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

# B. 
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
