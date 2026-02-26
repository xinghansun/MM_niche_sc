suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(GSVA)
  library(survival)
  library(survminer)
  library(edgeR)
})

# Download TCGA-LAML
query_exp <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts" 
)

GDCdownload(query_exp, directory = "../data/TCGA_LAML") 
tcga_data <- GDCprepare(query_exp, directory = "../data/TCGA_LAML")

# Prepare expression matrix
counts_matrix <- assay(tcga_data)
keep <- rowSums(cpm(counts_matrix) > 1) >= 5
counts_matrix <- counts_matrix[keep, ]
tcga_expr_matrix <- cpm(counts_matrix, log = TRUE) # to log2-CPM

# prepare clinical data
clinical <- as.data.frame(colData(tcga_data))

clinical_clean <- clinical %>%
  mutate(
    os_status = ifelse(vital_status == "Dead", 1, 0),
    os_days = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up)
  ) %>%
  filter(!is.na(os_days) & os_days > 0) %>%
  select(submitter_id, os_status, os_days, gender, age_at_index)

# pair expression and clinical
colnames(tcga_expr_matrix) <- substr(colnames(tcga_expr_matrix), 1, 12)
common_samples <- intersect(colnames(tcga_expr_matrix), clinical_clean$submitter_id)
tcga_expr_final <- tcga_expr_matrix[, common_samples]
clinical_final <- clinical_clean %>% filter(submitter_id %in% common_samples)
clinical_final <- clinical_final[match(common_samples, clinical_final$submitter_id), ]

# Convert Ensemble IDs to Gene Symbols
gene_metadata <- as.data.frame(rowData(tcga_data)) 
id_map <- gene_metadata %>% 
  dplyr::select(gene_id, gene_name)
valid_rows <- rownames(tcga_expr_final) %in% id_map$gene_id
tcga_expr_mapped <- tcga_expr_final[valid_rows, ]
match_idx <- match(rownames(tcga_expr_mapped), id_map$gene_id)
new_rownames <- id_map$gene_name[match_idx]
expr_df <- as.data.frame(tcga_expr_mapped)
expr_df$Symbol <- new_rownames
expr_clean <- expr_df %>%
  filter(!is.na(Symbol) & Symbol != "") %>%
  group_by(Symbol) %>%
  summarise(across(everything(), max))
final_mat <- as.matrix(expr_clean[, -1])
rownames(final_mat) <- expr_clean$Symbol


# Calculate gsva scores
model <- readRDS("../results/Myeloid_State_Risk_Model.rds")

percent_marker_included <- 100*length(intersect(model$gene_signature, rownames(final_mat)))/length(model$gene_signature)
message(sprintf("TCGA-LAML covers %.2f%% genes of the model (%d/%d)", 
                percent_marker_included, 
                length(intersect(model$gene_signature, rownames(final_mat))), 
                length(model$gene_signature)))

target_genes_list <- list(State_Score = model$gene_signature,
                          Myeloid_Abundance = c("CD14", "CSF1R", "CD68", "LYZ"))
gsva_res_tcga <- gsva(gsvaParam(
  expr = final_mat,
  geneSets = target_genes_list,
  kcdf = model$gsva_params$kcdf, 
  absRanking = model$gsva_params$absRanking
))

gsva_df_tcga <- as.data.frame(t(gsva_res_tcga))
gsva_df_tcga$submitter_id <- rownames(gsva_df_tcga)
validation_df <- left_join(clinical_final, gsva_df_tcga, by = "submitter_id")

# Test model
validation_df$group_fixed <- ifelse(
  validation_df$State_Score > model$risk_cutoff, 
  "High", "Low"
)

print(table(validation_df$group_fixed))

fit_fixed <- survfit(Surv(os_days, os_status) ~ group_fixed, data = validation_df)
svglite::svglite("../results/figures/05_MyeloidScore_TCGA-LAML_SurvivalCurve.svg", 
                 width = 7, height = 6)
ggsurvplot(fit_fixed, data = validation_df,
           pval = TRUE,
           risk.table = TRUE,
           palette = c("#2E9FDF", "#E7B800"), 
           #title = "Validation: Fixed Cutoff",
           xlab = "Days")
dev.off()

png("../results/figures/05_MyeloidScore_TCGA-LAML_SurvivalCurve.png", 
    width = 7, height = 6, units = "in", res = 600)
ggsurvplot(fit_fixed, data = validation_df,
           pval = TRUE,
           risk.table = TRUE,
           palette = c("#2E9FDF", "#E7B800"), 
           #title = "Validation: Fixed Cutoff",
           xlab = "Days")
dev.off()


# uni continuous
cox_fit_continuous <- coxph(Surv(os_days, os_status) ~ State_Score, data = validation_df)
summary(cox_fit_continuous)

zph <- cox.zph(cox_fit_continuous)
print(zph)
plot(zph)

# uni group
cox_uni_group <- coxph(Surv(os_days, os_status) ~ group_fixed, data = validation_df)
summary(cox_uni_group)

zph <- cox.zph(cox_uni_group)
print(zph)
plot(zph)

# multi
cox_multi_group <- coxph(Surv(os_days, os_status) ~ group_fixed + Myeloid_Abundance, 
                         data = validation_df)
summary(cox_multi_group)

zph <- cox.zph(cox_multi_group)
print(zph)
plot(zph)






