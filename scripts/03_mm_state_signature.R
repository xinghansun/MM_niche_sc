suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(stringr)
  library(data.table)
  library(harmony)
  library(ggplot2)
})

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

#SUP_GENES_REF <- c("CD274", "PDCD1LG2", "IDO1", "LGALS9", "LAIR1", "HAVCR2")
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
myeloid_markers <- c("CD68", "CD14", "FCGR3A", 
                     "LYZ", "S100A8", "S100A9", 
                     "CSF1R", "ITGAX")
combined <- AddModuleScore(
  object = combined,
  features = list(myeloid_markers),
  name = "Myeloid_Score"
)
VlnPlot(combined, features = "Myeloid_Score1", group.by = "seurat_clusters")
FeaturePlot(combined, features = "Myeloid_Score1")

myeloid_obj <- subset(combined, idents = c(1, 5, 8))
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
  
  # only use normal myeloid cells to define thresholds
  norm_cells <- which(obj$disease_stage == "Normal")
  ap_thr  <- median(obj$AP_val[norm_cells], na.rm = TRUE)
  sup_thr <- median(obj$SUP_val[norm_cells], na.rm = TRUE)
  
  obj$functional_state <- "Intermediate"
  obj$functional_state[obj$AP_val <= ap_thr & obj$SUP_val >= sup_thr] <- "Suppressive_like"
  obj$functional_state[obj$AP_val > ap_thr & obj$SUP_val < sup_thr] <- "Activated_like"
  
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
        group.by = "functional_state", 
        cols = c("red", "blue", "grey"))




