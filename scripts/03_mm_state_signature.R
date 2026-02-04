suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(DESeq2)
  library(limma)
  library(stringr)
  library(data.table)
})

dir.create("../results/scRNA", recursive = TRUE, showWarnings = FALSE)
dir.create("../results/signature", recursive = TRUE, showWarnings = FALSE)


# Define gene sets for state scoring
AP_GENES <- c(
  "CIITA",
  "CD74",
  "HLA-DRA","HLA-DRB1","HLA-DRB5",
  "HLA-DPA1","HLA-DPB1",
  "HLA-DQA1","HLA-DQB1",
  "HLA-DMA","HLA-DMB",
  "HLA-DOA","HLA-DOB",
  "CTSS"
)

SUP_GENES <- c(
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

# Integrate samples
obj_list <- lapply(obj_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x, verbose = FALSE)
  x <- RunPCA(x, features = VariableFeatures(x), npcs = 30, verbose = FALSE)
  return(x)
})

features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 2000)

anchors <- FindIntegrationAnchors(object.list = obj_list, 
                                  anchor.features = features,
                                  dims = 1:30)
rm(obj_list)
gc()
combined <- IntegrateData(anchorset = anchors, dims = 1:30)



# UMAP
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined, dims = 1:30)
ElbowPlot(combined, ndims = 30)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.2)

UMAPPlot(combined, label = TRUE)
DimPlot(combined, group.by = "orig.ident")

AP_GENES  <- intersect(AP_GENES, rownames(obj))
SUP_GENES <- intersect(SUP_GENES, rownames(obj))
