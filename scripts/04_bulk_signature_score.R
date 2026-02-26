suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(edgeR)
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
query_exp <- GDCquery(
  project = "MMRF-COMMPASS",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts" 
)
mmrf_data <- GDCprepare(query_exp, directory = "../data/MMRF_COMMPASS")
counts_matrix <- assay(mmrf_data)
keep <- rowSums(cpm(counts_matrix) > 1) >= 5
counts_matrix <- counts_matrix[keep, ]
expr_matrix <- cpm(counts_matrix, log = TRUE) # to log2-CPM

gene_metadata <- as.data.frame(rowData(mmrf_data)) 
id_map <- gene_metadata %>% 
  dplyr::select(gene_id, gene_name)
valid_rows <- rownames(expr_matrix) %in% id_map$gene_id
expr_matrix_mapped <- expr_matrix[valid_rows, ]
match_idx <- match(rownames(expr_matrix_mapped), id_map$gene_id)
new_rownames <- id_map$gene_name[match_idx]
expr_df <- as.data.frame(expr_matrix_mapped)
expr_df$Symbol <- new_rownames
expr_clean <- expr_df %>%
  filter(!is.na(Symbol) & Symbol != "") %>%
  group_by(Symbol) %>%
  summarise(across(everything(), max)) # Convert Ensemble IDs to Gene Symbols
final_mat <- as.matrix(expr_clean[, -1])
rownames(final_mat) <- expr_clean$Symbol
expr_matrix <- final_mat

clinical <- as.data.frame(colData(mmrf_data))
clinical_clean <- clinical %>%
  mutate(
    os_status = ifelse(vital_status == "Dead", 1, 0),
    os_days = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up)
  ) %>%
  filter(!is.na(os_days) & os_days > 0) %>%
  select(submitter_id, os_status, os_days, gender, age_at_index)

# Define signature
apsup_novel_pos <- read.table(
  "../results/myeloid_APSupHigh_novel_markers.tsv",
  header = TRUE)

refined_novel_genes <- apsup_novel_pos %>%
  filter(p_val_adj < 5e-8,
         avg_log2FC >= 0.5
         )
print(nrow(refined_novel_genes))
print(rownames(refined_novel_genes))

final_genes <- list(
  State_Score = rownames(refined_novel_genes),
  # to control for myeloid abundance later
  Myeloid_Abundance = c("CD14", "CSF1R", "CD68", "LYZ")
)

# Compute GSVA scores
params <- gsvaParam(
  expr = as.matrix(expr_matrix), 
  geneSets = final_genes,
  kcdf = "Gaussian",
  absRanking = FALSE
)

gsva_results <- gsva(params)

# survival analysis
gsva_df <- as.data.frame(t(gsva_results))
if(!all(rownames(clinical_clean) == rownames(gsva_df))) {
  gsva_df <- gsva_df[rownames(clinical_clean), ]
}

final_df <- cbind(clinical_clean, gsva_df)

res.cut <- surv_cutpoint(final_df, 
                         time = "os_days", 
                         event = "os_status", 
                         variables = "State_Score")

final_df$group <- ifelse(final_df$State_Score > res.cut$cutpoint$cutpoint, "High", "Low") # using res.cut
print(table(final_df$group))
fit <- survfit(Surv(os_days, os_status) ~ group, data = final_df)

svglite::svglite("../results/figures/04_MyeloidScore_MMRF-COMMPASS_SurvivalCurve.svg", 
                 width = 7, height = 6)
ggsurvplot(fit, data = final_df,
           pval = TRUE,
           risk.table = TRUE,
           palette = c("#2E9FDF", "#E7B800"), 
           #title = "Survival Analysis: 10-Gene Myeloid State",
           xlab = "Days")
dev.off()

png("../results/figures/04_MyeloidScore_MMRF-COMMPASS_SurvivalCurve.png", 
    width = 7, height = 6, units = "in", res = 600)
ggsurvplot(fit, data = final_df,
           pval = TRUE,
           risk.table = TRUE,
           palette = c("#2E9FDF", "#E7B800"), 
           #title = "Survival Analysis: 10-Gene Myeloid State",
           xlab = "Days")
dev.off()

# Cox regression to confirm independence from myeloid abundance
final_df$group <- factor(final_df$group, levels = c("Low", "High"))
# uni
cox_uni_group <- coxph(Surv(os_days, os_status) ~ group, data = final_df)
summary(cox_uni_group)

zph <- cox.zph(cox_uni_group)
print(zph)
plot(zph)

# multi
cox_multi_group <- coxph(Surv(os_days, os_status) ~ group + Myeloid_Abundance, 
                         data = final_df)
summary(cox_multi_group)

zph <- cox.zph(cox_multi_group)
print(zph)
plot(zph)


# Save model for validation
myeloid_model <- list(
  gene_signature = final_genes$State_Score, 
  risk_cutoff = res.cut$cutpoint$cutpoint,
  gsva_params = list(
    kcdf = "Gaussian",
    absRanking = FALSE
  ),
  source_note = "Derived from MMRF-COMMPASS using surv_cutpoint maxstat"
)
saveRDS(myeloid_model, "../results/Myeloid_State_Risk_Model.rds")



