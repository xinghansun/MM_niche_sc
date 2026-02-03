suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(stringr)
})

dir.create("../data/untared_geo", recursive = TRUE, showWarnings = FALSE)
dir.create("../data/seurat", recursive = TRUE, showWarnings = FALSE)

# GSE124310
message("Processing GSE124310 ...")

untar("../data/geo/GSE124310/GSE124310_RAW.tar", 
      exdir = "../data/untared_geo/GSE124310")

files <- list.files(path = "../data/untared_geo/GSE124310/", pattern = "filtered_gene_bc_matrices.tar.gz$")
for (f in files) {
  tared_file = paste0("../data/untared_geo/GSE124310/", f)
  sample <- sub("\\.filtered_gene_bc_matrices.tar.gz$", "", f)
  sample_dir <- paste0("../data/untared_geo/GSE124310/", sample)
  
  dir.create(sample_dir, showWarnings = FALSE)
  untar(tared_file, 
        exdir = sample_dir)
  file.remove(tared_file)
}

sample_dirs <- list.dirs(
  "../data/untared_geo/GSE124310",
  recursive = FALSE
)

objs <- lapply(sample_dirs, function(d) {
  data_dir <- file.path(d, "filtered_gene_bc_matrices", "GRCh38")
  
  counts <- Read10X(data_dir)
  
  CreateSeuratObject(
    counts = counts,
    project = basename(d)
  )
})

names(objs) <- basename(sample_dirs)
saveRDS(objs, "../data/seurat/GSE124310.seurat.list.rds")


# GSE163278 
message("Processing GSE163278 ...")

untar("../data/geo/GSE163278/GSE163278_RAW.tar", 
      exdir = "../data/untared_geo/GSE163278")

gz_files <- list.files("../data/untared_geo/GSE163278", 
                       pattern = "\\.gz$", 
                       full.names = TRUE)
for (f in gz_files) {
  R.utils::gunzip(f, overwrite = FALSE, remove = TRUE)
}

files <- list.files("../data/untared_geo/GSE163278", 
                    pattern = "\\.tsv$|\\.mtx$", 
                    full.names = TRUE)
samples <- unique(
  sub("_(barcodes|genes|matrix)\\.tsv$|_matrix\\.mtx$", 
      "", 
      basename(files)
      )
  )

for (s in samples) {
  sample_dir <- file.path("../data/untared_geo/GSE163278", s)
  dir.create(sample_dir, showWarnings = FALSE)

  for (f in files[grepl(s, files)]) {
    
    file.rename(
      f, 
      file.path(sample_dir, 
                tail(strsplit(basename(f), "_")[[1]], 1) # to fit Read10x
                )
      )
  }
}

sample_dirs <- list.dirs(
  "../data/untared_geo/GSE163278",
  recursive = FALSE
)

objs <- lapply(sample_dirs, function(d) {
  counts <- Read10X(d)
  
  CreateSeuratObject(
    counts = counts,
    project = basename(d)
  )
})

names(objs) <- basename(sample_dirs)
saveRDS(objs, "../data/seurat/GSE163278.seurat.list.rds")




