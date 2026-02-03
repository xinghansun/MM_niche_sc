suppressPackageStartupMessages({
  library(GEOquery)
})

dir.create("../data/geo", recursive = TRUE, showWarnings = FALSE)

download_geo_supp <- function(gse, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("Downloading supplementary files for: ", gse)
  
  GEOquery::getGEOSuppFiles(
    GEO = gse,
    baseDir = outdir,
    makeDirectory = FALSE,
    filter_regex = NULL
  )
  
  message("Done: ", gse)
}

download_geo_matrix <- function(gse, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  message("Downloading series matrix for: ", gse)
  g <- GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = TRUE)
  saveRDS(g, file = file.path(outdir, paste0(gse, "_GSEMatrix.rds")))
  message("Done: ", gse)
}

# scRNA
download_geo_supp("GSE124310", "../data/geo/GSE124310") # discovery scRNA
download_geo_supp("GSE163278", "../data/geo/GSE163278") # validation scRNA

# bulk
download_geo_matrix("GSE136337", "../data/geo/GSE136337") # validation bulk RNA

message("All downloads finished.")
