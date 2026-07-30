############################################################################################
#Packages
library(readxl)
library(ShortRead)
library(ggplot2)
library(dplyr)
library(tidyr)
library(brms)
library(matrixStats)

############################################################################################
# Datasets
# Public (NCBI)
# Metadata
exp1 <- read.csv2("C:/files/sra_runinfo_PRJNA513137.csv",sep=",")
exp2 <- read.csv2("C:/files/sra_runinfo_PRJNA1072695.csv",sep=",")
exp <- rbind(exp1,exp2)
rm(list = setdiff(ls(), "exp"))
exp <- exp[,colnames(exp) %in% c("Run","spots","BioSample")]
colnames(exp) <- c("NCBI_Run","NCBI_raw_reads","NCBI_BioSample")

# Downloaded
path <- "C:/files/merged_vsearch"
sample_names <- sort(list.files(path, pattern=".assembled.fastq", full.names = FALSE))
sample_names <- sub(".assembled.fastq$", "", sample_names)
merged <- sort(list.files(path, pattern=".assembled.fastq", full.names = TRUE))
not_merged_r1 <- sort(list.files(path, pattern=".notmerged_R1.fastq", full.names = TRUE))
not_merged_r2 <- sort(list.files(path, pattern=".notmerged_R2.fastq", full.names = TRUE))
merge_summary <- data.frame(
  Downloaded_Run = sample_names,
  merge = countFastq(merged)$records,
  nmerged1 = countFastq(not_merged_r1)$records,
  nmerged2 = countFastq(not_merged_r2)$records
)
merge_summary$Downloaded_raw_reads <- merge_summary$merge + merge_summary$nmerged1
merge_summary <- merge_summary[,c(1,5)]
metadata <- merge(exp,merge_summary,by.x="NCBI_Run",by.y="Downloaded_Run")

# Supplementary data
supplementar_data <- read_excel("C:/dePCR/Table_S5_Meta_Data.xlsx")
supplementar_data <- supplementar_data[,c(4,13,18)]
colnames(supplementar_data) <- c("Supplementary_BioSample","Supplementary_Name","Supplementary_raw_reads")
metadata <- merge(metadata,supplementar_data,by.x="NCBI_BioSample",by.y="Supplementary_BioSample")
metadata$comparison_raw_reads <- apply(
  metadata[, c("Downloaded_raw_reads", "Supplementary_raw_reads")],1,function(x) x[1] == x[2]
)

############################################################################################
# Count reads (Published data)
# count_raw_reads_paper <- read.csv2("C:/dePCR/Supplemental_Material_3/data_tables/raw_feature_table.csv",sep=",")
# count_raw_reads_paper <- merge(metadata,
#                                count_raw_reads_paper,
#                                by.x="Supplementary_Name",
#                                by.y="sample_name")
# count_raw_reads_paper$percent <- count_raw_reads_paper$total_reads/count_raw_reads_paper$Supplementary_raw_reads
# cores <- rep("lightblue", length(unique(count_raw_reads_paper$Exp)))
# names(cores) <- sort(unique(count_raw_reads_paper$Exp))
# cores[c("A10", "A64")] <- "red"
# par(mfrow = c(1, 2), mar = c(8, 4, 4, 1))
# boxplot(
#   percent ~ Exp,
#   data = count_raw_reads_paper,
#   xlab = "Experiment",
#   ylab = "%",
#   main = "Mapped reads by experiment (Primer-Template) - Their results",
#   col = cores
# )

# Our results
our_results <- read.csv2("C:/files/FinalCount_samples_rows2.tsv",sep="")
our_results$total_counts <- rowSums(our_results[, 2:641])
our_results <- merge(metadata,our_results,by.x="NCBI_Run",by.y="sample")
our_results$Exp <- sub("_.*", "", our_results$Supplementary_Name)
our_results$percent <- our_results$total_counts/our_results$Downloaded_raw_reads

# boxplot(
#   percent ~ Exp,
#   data = our_results,
#   xlab = "Experiment",
#   ylab = "%",
#   main = "Mapped reads by experiment (Primer-Template) - Our results",
#   col = cores
# )
rm(list = setdiff(ls(), "our_results"))

# write.table(metadata, "C:/files/comparison.txt",
#           row.names=FALSE,
#           col.names=TRUE,
#           quote=FALSE,
#           sep="\t")

############################################################################################
# Doing regression analysis
our_results <- our_results[,c(5,8:648)]
our_results <- our_results %>%
  separate(
    col = Supplementary_Name,
    into = c("Exp", "Method", "Enzyme", "Temp", "Replicate"),
    sep = "_",
    remove = FALSE
  )
our_results <- our_results %>%
  pivot_longer(
    cols = 7:646,
    names_to = "feature",
    values_to = "Counts"
  )
our_results$id_feature <- paste(our_results$Exp,our_results$feature,sep="_")
metadados <- read_excel("C:/dePCR/dePCR_feature_metadata.xlsx")
metadados_feature <- read.table("C:/dePCR/dePCR_experiment_features.txt", header = TRUE)
metadados_feature$id_feature <- paste(metadados_feature$CurrentName,metadados_feature$feature,sep="_")
our_results_filtered <- our_results[our_results$id_feature %in% unique(metadados_feature$id_feature),]
our_results_filtered <- merge(
  our_results_filtered,
  metadados,
  by = "feature"
)
colnames(our_results_filtered)[12] <- "p3_status"
colnames(our_results_filtered)[14] <- "p5_status"
colnames(our_results_filtered)[15] <- "p3_base"
colnames(our_results_filtered)[17] <- "p5_base"
our_results_filtered$Temp <- factor(our_results_filtered$Temp,levels = c("45C", "55C"))
our_results_filtered$p3_status <- factor(our_results_filtered$p3_status,levels = c("Match", "Mismatch"))
our_results_filtered$mid_status <- factor(our_results_filtered$mid_status,levels = c("Match", "Mismatch"))
our_results_filtered$p5_status <- factor(our_results_filtered$p5_status,levels = c("Match", "Mismatch"))
our_results_filtered$p3_base <- factor(our_results_filtered$p3_base,levels = c("A", "C", "G", "T"))
our_results_filtered$mid_base <- factor(our_results_filtered$mid_base,levels = c("T", "C", "G", "A"))
our_results_filtered$p5_base <- factor(our_results_filtered$p5_base,levels = c("C", "A", "G", "T"))
our_results_filtered$group <- paste(our_results_filtered$p3_base,
                                    our_results_filtered$mid_base,
                                    our_results_filtered$p5_base,
                                    sep="_")
our_results_filtered$group <- factor(our_results_filtered$group,levels = c("A_T_C",
                                                                           "A_G_C",
                                                                           "A_T_G",
                                                                           "A_A_C",
                                                                           "A_C_C",
                                                                           "A_T_A",
                                                                           "A_T_T",
                                                                           "C_T_C",
                                                                           "G_T_C",
                                                                           "T_T_C"))
rm(list = setdiff(ls(), "our_results_filtered"))

###############################################################################
# Bayesian modelling
# Data spliting
# PCR
set.seed(10231991)
exp_10 <- our_results_filtered[our_results_filtered$Exp=="A10",]
#exp_10 <- exp_10[exp_10$Temp=="55C",]
#exp_10_55 <- exp_10[exp_10$Temp=="55C",]
exp_10$log_total_counted <- log(exp_10$total_counts)
exp_10_pcr <- exp_10[exp_10$Method=="PCR",]
#exp_10_pcr_train <- exp_10_pcr[exp_10_pcr$Replicate %in% as.character(sample(1:8, 6)),]
#exp_10_pcr_test <- exp_10_pcr[!(exp_10_pcr$Replicate %in% exp_10_pcr_train$Replicate),]

# DePCR
exp_10_depcr <- exp_10[exp_10$Method=="DePCR",]
#exp_10_depcr_train <- exp_10_depcr[exp_10_depcr$Replicate %in% as.character(sample(1:8, 6)),]
#exp_10_depcr_test <- exp_10_depcr[!(exp_10_depcr$Replicate %in% exp_10_depcr_train$Replicate),]

# ################################################################################
# One by one
# Matrix model
# bf_disp_1 <- bf(
#   Counts ~ offset(log_total_counted),
#   shape ~ 1
# )
bf_disp_status <- bf(
  Counts ~ offset(log_total_counted) + p3_status + mid_status + p5_status,
  shape ~ p3_base + mid_base + p5_base + Temp
)
# fit_depcr <- brm(
#     formula = bf_disp_1,
#     data = exp_10_depcr_train,
#     family = negbinomial(link = "log", link_shape = "log"),
#     chains = 4,
#     cores = 8,
#     iter = 6000,
#     warmup = 1000,
#     control = list(adapt_delta = 0.99, max_treedepth = 12),
#     refresh = 50,
#     seed = 10231991
#   )
# summary(fit_depcr)
# loo_1 <- loo(fit_depcr)

fit_depcr_status <- brm(
  formula = bf_disp_status,
  data = exp_10_depcr,
  family = negbinomial(link = "log", link_shape = "log"),
  chains = 4,
  cores = 8,
  iter = 6000,
  warmup = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  refresh = 50,
  seed = 10231991
)
summary(fit_depcr_status)

loo_status <- loo(fit_depcr_status)

results_loo <- loo_compare(loo_1, loo_status)
summary(fit_depcr_status)
summary(fit_depcr)

################################################################################
# Using the trained model on test set
epred_test <- posterior_epred(
  fit_depcr_status,
  newdata = exp_10_depcr_test
)
pred_mean_test <- colMeans(epred_test)
pred_mean_ci <- apply(
  epred_test,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)
pred_counts_test <- posterior_predict(
  fit_depcr_status,
  newdata = exp_10_depcr_test
)
pred_counts_pi <- apply(
  pred_counts_test,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)
)
obs_test <- exp_10_depcr_test$Counts
rmse_test <- sqrt(mean((obs_test - pred_mean_test)^2))
mae_test <- mean(abs(obs_test - pred_mean_test))
cor_pearson_test <- cor(obs_test, pred_mean_test, method = "pearson")
rmse_test
mae_test
cor_pearson_test
coverage_95 <- mean(
  obs_test >= pred_counts_pi[1, ] &
    obs_test <= pred_counts_pi[3, ]
)
log_lik_test <- log_lik(
  fit_depcr_full,
  newdata = exp_10_depcr_test
)
lpd_i_test <- matrixStats::colLogSumExps(log_lik_test) -
  log(nrow(log_lik_test))
lpd_test <- sum(lpd_i_test)
mean_lpd_test <- mean(lpd_i_test)
lpd_test
mean_lpd_test

test_predictions <- data.frame(
  observed = obs_test,
  predicted_mean = pred_mean_test,
  predicted_median_count = pred_counts_pi[2, ],
  pred_lower_95 = pred_counts_pi[1, ],
  pred_upper_95 = pred_counts_pi[3, ],
  p3_status = exp_10_depcr_test$p3_status,
  mid_status = exp_10_depcr_test$mid_status,
  p5_status = exp_10_depcr_test$p5_status,
  p3_base = exp_10_depcr_test$p3_base,
  mid_base = exp_10_depcr_test$mid_base,
  p5_base = exp_10_depcr_test$p5_base
)
coverage_95 <- mean(
  obs_test >= pred_counts_pi[1, ] &
    obs_test <= pred_counts_pi[3, ]
)
coverage_95

ggplot(test_predictions, aes(x = observed, y = predicted_mean)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  labs(
    x = "Observed counts",
    y = "Predicted expected counts",
    title = "Observed vs predicted counts - test set, log scale"
  )

test_predictions$status_combo <- paste(
  test_predictions$p3_base,
  test_predictions$mid_base,
  test_predictions$p5_base,
  sep = "_"
)

ggplot(test_predictions, aes(x = observed, y = predicted_mean, color = status_combo)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme_bw() +
  labs(
    x = "Observed counts",
    y = "Predicted expected counts",
    color = "Status combo",
    title = "Observed vs predicted counts by mismatch status"
  )

################################################################################
# Looping
# out_dir_pcr <- "C:/files/brms_A10_PCR_Group"
# dir.create(out_dir_pcr, recursive = TRUE, showWarnings = FALSE)
# out_dir_depcr <- "C:/files/brms_A10_DePCR_Group"
# dir.create(out_dir_depcr, recursive = TRUE, showWarnings = FALSE)
# mean_models <- list(
#   base = "Counts ~ offset(log_total_counted) + group"
# )
# shape_models <- list(
#   shape_Temp = "Temp",
#   shape_group = "group",
#   shape_p3_status = "p3_status",
#   shape_mid_status = "mid_status",
#   shape_p5_status = "p5_status",
#   shape_status = "p3_status + mid_status + p5_status",
#   shape_p3_base = "p3_base",
#   shape_mid_base = "mid_base",
#   shape_p5_base = "p5_base",
#   shape_base = "p3_base + mid_base + p5_base",
#   # shape_Temp_p3_status = "Temp + p3_status",
#   # shape_Temp_mid_status = "Temp + mid_status",
#   # shape_Temp_p5_status = "Temp + p5_status",
#   # shape_Temp_status = "Temp + p3_status + mid_status + p5_status",
#   # shape_Temp_p3_base = "Temp + p3_base",
#   # shape_Temp_mid_base = "Temp + mid_base",
#   # shape_Temp_p5_base = "Temp + p5_base",
#   # shape_Temp_base = "Temp + p3_base + mid_base + p5_base",
#   # shape_Temp_group = "Temp + group",
#   shape_1 = "1"
# )
# model_grid <- data.frame(
#   mean_name = character(),
#   shape_name = character(),
#   stringsAsFactors = FALSE
# )
# model_grid <- rbind(
#   model_grid,
#   data.frame(
#     mean_name = "base",
#     shape_name = "shape_1",
#     stringsAsFactors = FALSE
#   )
# )
# model_grid <- rbind(
#   model_grid,
#   data.frame(
#     mean_name = "base",
#     shape_name = names(shape_models),
#     stringsAsFactors = FALSE
#   )
# )
# fits <- list()
# for (i in seq_len(nrow(model_grid))) {
#   mean_name <- model_grid$mean_name[i]
#   shape_name <- model_grid$shape_name[i]
#   model_id_pcr <- paste(mean_name, shape_name, "PCR", sep = "_")
#   model_id_depcr <- paste(mean_name, shape_name, "DePCR", sep = "_")
#   bf_disp <- bf(
#     as.formula(mean_models[[mean_name]]),
#     as.formula(paste("shape ~", shape_models[[shape_name]]))
#   )
#   cat("\n====================================\n")
#   cat("Running model:", model_id_pcr, "\n")
#   cat("====================================\n")
#   fit_pcr <- tryCatch(
#     brm(
#       formula = bf_disp,
#       data = exp_10_pcr_train,
#       family = negbinomial(link = "log", link_shape = "log"),
#       chains = 4,
#       cores = 4,
#       iter = 6000,
#       warmup = 1000,
#       control = list(adapt_delta = 0.99, max_treedepth = 12),
#       refresh = 50,
#       seed = 10231991
#     ),
#     error = function(e) {
#       cat("\nERROR in model:", model_id_pcr, "\n")
#       cat(e$message, "\n")
#       return(NULL)
#     }
#   )
#   fits[[model_id_pcr]] <- fit_pcr
#   if (!is.null(fit_pcr)) {
#     # saveRDS(
#     #   fit_pcr,
#     #   file = file.path(out_dir_pcr, paste0(model_id_pcr, "base.rds"))
#     # )
#     writeLines(
#       capture.output(summary(fit_pcr)),
#       con = file.path(out_dir_pcr, paste0(model_id_pcr, "_summary_base.txt"))
#     )
#   }
#   cat("\n====================================\n")
#   cat("Running model:", model_id_depcr, "\n")
#   cat("====================================\n")
#   fit_depcr <- tryCatch(
#     brm(
#       formula = bf_disp,
#       data = exp_10_depcr_train,
#       family = negbinomial(link = "log", link_shape = "log"),
#       chains = 4,
#       cores = 8,
#       iter = 6000,
#       warmup = 1000,
#       control = list(adapt_delta = 0.99, max_treedepth = 12),
#       refresh = 50,
#       seed = 10231991
#     ),
#     error = function(e) {
#       cat("\nERROR in model:", model_id_depcr, "\n")
#       cat(e$message, "\n")
#       return(NULL)
#     }
#   )
#   fits[[model_id_depcr]] <- fit_depcr
#   if (!is.null(fit_depcr)) {
#     saveRDS(
#       fit_depcr,
#       file = file.path(out_dir_depcr, paste0(model_id_depcr, "base.rds"))
#     )
#     writeLines(
#       capture.output(summary(fit_depcr)),
#       con = file.path(out_dir_depcr, paste0(model_id_depcr, "_summary_base.txt"))
#     )
#   }
# }
