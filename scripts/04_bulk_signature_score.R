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

# Define signature
apsup_novel_pos <- read.table(
  "../results/myeloid_APSupHigh_novel_markers.tsv",
  header = TRUE)

refined_novel_genes <- apsup_novel_pos %>%
  filter(p_val_adj < 1e-5,
         avg_log2FC >= 0.5
         )

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
if(!all(rownames(clinical_meta) == rownames(gsva_df))) {
  gsva_df <- gsva_df[rownames(clinical_meta), ]
}

final_df <- cbind(clinical_meta, gsva_df)

df_surv_clean <- final_df %>%
  mutate(
    date_sample = ymd(`sampledatetime__1:ch1`), 
    clean_death_num = as.numeric(ifelse(`datedeath:ch1` == ".", NA, `datedeath:ch1`)),
    date_death = as.Date(clean_death_num, origin = "1899-12-30"),
    os_status = ifelse(is.na(date_death), 0, 1),
    days_to_death = as.numeric(date_death - date_sample)
  ) 

df_ready <- final_df %>%
  mutate(
    os_months = as.numeric(ifelse(`monthsos:ch1` == ".", NA, `monthsos:ch1`)),
    os_days = os_months * 30.4,
    
    os_status = ifelse(`datedeath:ch1` == ".", 0, 1)
  ) %>%
  filter(!is.na(os_days) & os_days > 0)

res.cut <- surv_cutpoint(df_ready, 
                         time = "os_days", 
                         event = "os_status", 
                         variables = "State_Score")

df_ready$group <- ifelse(df_ready$State_Score > res.cut$cutpoint$cutpoint, "High", "Low") # using res.cut
print(table(df_ready$group))
fit <- survfit(Surv(os_days, os_status) ~ group, data = df_ready)

svglite::svglite("../results/figures/04_MyeloidScore_GSE136337_SurvivalCurve.svg", 
                 width = 7, height = 6)
ggsurvplot(fit, data = df_ready,
           pval = TRUE,
           risk.table = TRUE,
           palette = c("#2E9FDF", "#E7B800"), 
           #title = "Survival Analysis: 10-Gene Myeloid State",
           xlab = "Days")
dev.off()

png("../results/figures/04_MyeloidScore_GSE136337_SurvivalCurve.png", 
    width = 7, height = 6, units = "in", res = 600)
ggsurvplot(fit, data = df_ready,
           pval = TRUE,
           risk.table = TRUE,
           palette = c("#2E9FDF", "#E7B800"), 
           #title = "Survival Analysis: 10-Gene Myeloid State",
           xlab = "Days")
dev.off()

# Cox regression to confirm independence from myeloid abundance
df_ready$group <- factor(df_ready$group, levels = c("Low", "High"))
# uni
cox_uni_group <- coxph(Surv(os_days, os_status) ~ group, data = df_ready)
summary(cox_uni_group)

zph <- cox.zph(cox_uni_group)
print(zph)
plot(zph)

# multi
cox_multi_group <- coxph(Surv(os_days, os_status) ~ group + Myeloid_Abundance, 
                         data = df_ready)
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
  source_note = "Derived from GSE136337 using surv_cutpoint maxstat"
)
saveRDS(myeloid_model, "../results/Myeloid_State_Risk_Model.rds")



