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
  library(broom)
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
           #title = "Survival Analysis: 7-Gene Myeloid State",
           xlab = "Days")
dev.off()

png("../results/figures/04_MyeloidScore_MMRF-COMMPASS_SurvivalCurve.png", 
    width = 7, height = 6, units = "in", res = 600)
ggsurvplot(fit, data = final_df,
           pval = TRUE,
           risk.table = TRUE,
           palette = c("#2E9FDF", "#E7B800"), 
           #title = "Survival Analysis: 7-Gene Myeloid State",
           xlab = "Days")
dev.off()

# Cox regression to confirm independence from myeloid abundance
final_df$group <- factor(final_df$group, levels = c("Low", "High"))
# uni binary
cox_uni_group <- coxph(Surv(os_days, os_status) ~ group, data = final_df)
summary(cox_uni_group)

zph <- cox.zph(cox_uni_group)
print(zph)
plot(zph)

# multi binary
cox_multi_group <- coxph(Surv(os_days, os_status) ~ group + 
                                                    Myeloid_Abundance +
                                                    gender + age_at_index, 
                         data = final_df)
summary(cox_multi_group)

zph <- cox.zph(cox_multi_group)
print(zph)
plot(zph)

# multi continuous
cox_multi_continuous <- coxph(Surv(os_days, os_status) ~ State_Score + 
                           Myeloid_Abundance +
                           gender + age_at_index, 
                         data = final_df)
summary(cox_multi_continuous)

zph <- cox.zph(cox_multi_continuous)
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


# Plot State Score Distribution
cutoff_value <- res.cut$cutpoint$cutpoint
col_density <- "#4C72B0"
col_cutoff  <- "#D62728"

p <- ggplot(final_df, aes(x = State_Score)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 fill = "grey80",
                 color = "black",
                 linewidth = 0.2,
                 alpha = 0.6) +
  geom_density(color = col_density,
               linewidth = 1.2,
               adjust = 1.2) +
  geom_vline(xintercept = cutoff_value,
             linetype = "dashed",
             linewidth = 1,
             color = col_cutoff) +
  annotate("text",
           x = cutoff_value,
           y = Inf,
           label = paste0("Cutoff = ", round(cutoff_value, 3)),
           vjust = 2,
           hjust = -0.05,
           size = 5,
           color = col_cutoff,
           fontface = "bold") +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  labs(
    x = "State Score",
    y = "Density"
  )

svglite::svglite("../results/figures/04_MyeloidScore_MMRF-COMMPASS_StateScoreDistribution.svg", 
                 width = 6, height = 5)
p
dev.off()
png("../results/figures/04_MyeloidScore_MMRF-COMMPASS_StateScoreDistribution.png", 
    width = 6, height = 5, units = "in", res = 600)
p
dev.off()


# Plot Forest Plot for Multi-var Cox Model
cox_summary <- tidy(cox_multi_continuous, 
                    exponentiate = TRUE, 
                    conf.int = TRUE)

forest_df <- cox_summary %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = case_when(
      term == "State_Score" ~ "AP/SUP-High State Score",
      term == "Myeloid_Abundance" ~ "Myeloid Abundance Score",
      term == "gendermale" ~ "Gender (Male vs Female)",
      term == "age_at_index" ~ "Age (Year)",
      TRUE ~ term
    ),
    term = factor(term, levels = rev(term)),
    p_label = ifelse(p.value < 0.001,
                     "P < 0.001",
                     paste0("P = ", sprintf("%.3f", p.value)))
  )

x_max <- max(forest_df$conf.high) * 1.6
x_min <- min(forest_df$conf.low) * 0.8

p_forest <- ggplot(forest_df, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             linewidth = 0.8,
             color = "grey40") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2,
                 linewidth = 1) +
  geom_point(size = 3) +
  geom_text(aes(x = x_max, label = p_label),
            hjust = 1,
            size = 5) +
  scale_x_log10(limits = c(x_min, x_max)) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  labs(
    x = "Hazard Ratio (log scale)",
    y = NULL
  )

svglite::svglite("../results/figures/04_MyeloidScore_MMRF-COMMPASS_CoxForestPlot.svg", 
                 width = 10, height = 5)
p_forest
dev.off()
png("../results/figures/04_MyeloidScore_MMRF-COMMPASS_CoxForestPlot.png", 
    width = 10, height = 5, units = "in", res = 600)
p_forest
dev.off()


# Export Multi-var Cox Model Summary
cox_table <- tidy(cox_multi_continuous,
                  exponentiate = TRUE,
                  conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Variable = case_when(
      term == "State_Score" ~ "AP/SUP-High State Score",
      term == "Myeloid_Abundance" ~ "Myeloid Abundance Score",
      term == "gendermale" ~ "Gender (Male vs Female)",
      term == "age_at_index" ~ "Age (Year)",
      TRUE ~ term
    ),
    `Hazard Ratio` = round(estimate, 3),
    `Lower 95% CI` = round(conf.low, 3),
    `Upper 95% CI` = round(conf.high, 3),
    `P value` = p.value
  ) %>%
  select(Variable,
         `Hazard Ratio`,
         `Lower 95% CI`,
         `Upper 95% CI`,
         `P value`)

write.csv(cox_table,
          "../results/04_MyeloidScore_MMRF-COMMPASS_MultiVarCoxModelSummary.csv",
          row.names = FALSE)

# Export AP/SUP-High State Score Signature Gene List
write.table(final_genes$State_Score,
            "../results/04_MyeloidScore_MMRF-COMMPASS_StateScoreSignatureGenes.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
