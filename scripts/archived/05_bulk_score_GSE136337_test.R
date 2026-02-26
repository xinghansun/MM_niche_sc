suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(Biobase)
  library(GSVA)
  library(survival)
  library(ggplot2)
  library(lubridate)
  library(survminer)
  library(svglite)
})

# Load expression data
expr_data <- readRDS("../data/geo/GSE136337/GSE136337_GSEMatrix.rds")
expr_matrix <- exprs(expr_data$GSE136337_series_matrix.txt.gz)

clinical_meta <- pData(expr_data$GSE136337_series_matrix.txt.gz)
clinical_df <- clinical_meta %>%
  mutate(
    os_months = as.numeric(ifelse(`monthsos:ch1` == ".", NA, `monthsos:ch1`)),
    os_days = os_months * 30.44,
    os_status = ifelse(`datedeath:ch1` == ".", 0, 1),
    submitter_id = geo_accession,
    gender = `gender:ch1`,
    age_at_index = as.numeric(`ageatsampledate:ch1`)
  ) %>%
  filter(!is.na(os_days) & os_days > 0)


# Calculate gsva scores
model <- readRDS("../results/Myeloid_State_Risk_Model.rds")

percent_marker_included <- 100*length(intersect(model$gene_signature, rownames(expr_matrix)))/length(model$gene_signature)
message(sprintf("TCGA-LAML covers %.2f%% genes of the model (%d/%d)", 
                percent_marker_included, 
                length(intersect(model$gene_signature, rownames(expr_matrix))), 
                length(model$gene_signature)))

target_genes_list <- list(State_Score = model$gene_signature,
                          Myeloid_Abundance = c("CD14", "CSF1R", "CD68", "LYZ"))
gsva_res_tcga <- gsva(gsvaParam(
  expr = expr_matrix,
  geneSets = target_genes_list,
  kcdf = model$gsva_params$kcdf, 
  absRanking = model$gsva_params$absRanking
))

gsva_df_tcga <- as.data.frame(t(gsva_res_tcga))
gsva_df_tcga$submitter_id <- rownames(gsva_df_tcga)
validation_df <- left_join(clinical_df, gsva_df_tcga, by = "submitter_id")

# Test model
# multi continuous
cox_multi_continuous <- coxph(Surv(os_days, os_status) ~ State_Score + 
                                Myeloid_Abundance +
                                gender + age_at_index, 
                              data = validation_df)
summary(cox_multi_continuous)

zph <- cox.zph(cox_multi_continuous)
print(zph)
plot(zph)

