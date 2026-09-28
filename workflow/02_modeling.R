################################################################################
# PACKAGES
library(brms)
library(rstan)
library(loo)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)
set.seed(10231991)
options(mc.cores = parallel::detectCores())

################################################################################
# INPUT DATA
data_model <- data %>% filter(Exp %in% c("A10", "A64"), total_counted > 13000)
data_PCR <- data_model %>% filter(Protocol == "PCR")
data_DePCR <- data_model %>% filter(Protocol == "DePCR")

################################################################################
# DEFINE REFERENCE LEVELS
prepare_factors <- function(df) {
  df %>%
    mutate(
      Exp = relevel(factor(Exp), ref = "A10"),
      Temp = relevel(factor(Temp), ref = "45C"),
      group = relevel(factor(group), ref = "A_T_C"),
      p3_status = relevel(factor(p3_status), ref = "Match"),
      mid_status = relevel(factor(mid_status), ref = "Match"),
      p5_status = relevel(factor(p5_status), ref = "Match")
    )
}
data_PCR   <- prepare_factors(data_PCR)
data_DePCR <- prepare_factors(data_DePCR)

################################################################################
# MODEL SETTINGS
control_settings <- list(adapt_delta = 0.99, max_treedepth = 15)
final_formula <- bf(Counts ~ offset(log_total_counted) + Exp + Temp + group, shape ~ Exp + Temp + p3_status + mid_status + p5_status)

################################################################################
# FIT FINAL PCR MODEL
fit_PCR <- brm(
  formula = final_formula,
  data = data_PCR,
  family = negbinomial(link = "log"),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  seed = 10231991,
  backend = "rstan",
  control = control_settings,
  save_pars = save_pars(all = TRUE)
)

################################################################################
# FIT FINAL DePCR MODEL
fit_DePCR <- brm(
  formula = final_formula,
  data = data_DePCR,
  family = negbinomial(link = "log"),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  seed = 10231991,
  backend = "rstan",
  control = control_settings,
  save_pars = save_pars(all = TRUE)
)

################################################################################
# MODEL SUMMARIES
summary(fit_PCR)
summary(fit_DePCR)

################################################################################
# CONSTANT-SHAPE NEGATIVE BINOMIAL MODELS
constant_shape_formula <- bf(Counts ~ offset(log_total_counted) + Exp + Temp + group, shape ~ 1)
fit_PCR_NB <- brm(
  formula = constant_shape_formula,
  data = data_PCR,
  family = negbinomial(link = "log"),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  seed = 10231991,
  backend = "rstan",
  control = control_settings,
  save_pars = save_pars(all = TRUE)
)
fit_DePCR_NB <- brm(
  formula = constant_shape_formula,
  data = data_DePCR,
  family = negbinomial(link = "log"),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  seed = 10231991,
  backend = "rstan",
  control = control_settings,
  save_pars = save_pars(all = TRUE)
)

################################################################################
# PSIS-LOO
loo_PCR_distNB <- loo(fit_PCR)
loo_PCR_NB     <- loo(fit_PCR_NB)
loo_DePCR_distNB <- loo(fit_DePCR)
loo_DePCR_NB     <- loo(fit_DePCR_NB)

################################################################################
# DISTRIBUTIONAL VS CONSTANT-SHAPE MODEL
comp_PCR <- loo_compare(list(NB_distributional = loo_PCR_distNB, NB = loo_PCR_NB))
comp_DePCR <- loo_compare(list(NB_distributional = loo_DePCR_distNB, NB = loo_DePCR_NB))
print(comp_PCR)
print(comp_DePCR)

################################################################################
# SAVE
dir.create("results/model_fits", recursive = TRUE, showWarnings = FALSE)
saveRDS(fit_PCR, "results/model_fits/final_model_PCR.rds")
saveRDS(fit_DePCR, "results/model_fits/final_model_DePCR.rds")
saveRDS(fit_PCR_NB,"results/model_fits/constant_shape_NB_PCR.rds")
saveRDS(fit_DePCR_NB, "results/model_fits/constant_shape_NB_DePCR.rds")

dir.create("results/model_comparison", recursive = TRUE, showWarnings = FALSE)
write.table(as.data.frame(comp_PCR), "results/model_comparison/distributional_vs_NB_PCR.tsv", sep = "\t", quote = FALSE, row.names = TRUE)
write.table(as.data.frame(comp_DePCR), "results/model_comparison/distributional_vs_NB_DePCR.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

################################################################################
# POSTERIOR COEFFICIENT TABLES
PCR_parameters <- as.data.frame(posterior_summary(fit_PCR, probs = c(0.025, 0.975)))
DePCR_parameters <- as.data.frame(posterior_summary(fit_DePCR, probs = c(0.025, 0.975)))
PCR_parameters$parameter <- rownames(PCR_parameters)
DePCR_parameters$parameter <- rownames(DePCR_parameters)
rownames(PCR_parameters) <- NULL
rownames(DePCR_parameters) <- NULL
write.table(PCR_parameters, "results/model_comparison/final_PCR_parameters.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(DePCR_parameters, "results/model_comparison/final_DePCR_parameters.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

################################################################################
# EXPLORATORY TRANSFER DATA (ARRUMAR)
################################################################################
# 0. LIBRARY PATHS + PACKAGES + BUILD TOOLS
################################################################################
new_lib <- "/home/luangc/.ondemand/luangc/rstudio/libs/4.5.1-clean"
old_lib <- "/home/luangc/.ondemand/luangc/rstudio/libs/4.4.1"
if (dir.exists(new_lib)) {
  .libPaths(c(new_lib,setdiff(.libPaths(), old_lib)))
}
options(buildtools.check = NULL)
build_tools_ok <- pkgbuild::has_build_tools(debug = TRUE)
options(buildtools.check = function(action) TRUE)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 10)

library(brms)
library(rstan)
library(loo)
library(posterior)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
packageVersion("brms")
packageVersion("rstan")
packageVersion("StanHeaders")
################################################################################
# 1. DATA
################################################################################
# Importing data
data <- readRDS("./files/data.rds")
data <- data[data$total_counts > 13000,]
data$log_total_counted <- log(data$total_counts)
data$group <- paste(data$p3_base, data$mid_base, data$p5_base, sep="_")
data$ST <- stringr::str_extract(data$feature, "ST\\d{2}")

data_PCR <- data[data$Method=="PCR",]
data_PCR <- data_PCR[data_PCR$Exp %in% c("A10","A64"),]
 
data_DePCR <- data[data$Method=="DePCR",]
data_DePCR <- data_DePCR[data_DePCR$Exp %in% c("A10","A64"),]

################################################################################
# 2. DESIGN MATRIX
################################################################################
fixed_mean <- "offset(log_total_counted) + Exp + Temp"
fixed_shape <- "Exp + Temp"
var_status <- c("p3_status","mid_status","p5_status")
var_mut <- c("p3_mut","mid_mut","p5_mut")
var_base <- c("p3_base","mid_base","p5_base")

# Function to create the combinations
make_pairs <- function(v, name, fixed_mean, fixed_shape) {
  out <- do.call(
    rbind,
    lapply(seq_along(v), function(i) {
      value <- paste(v[1:i], collapse = " + ")
      data.frame(
        vector = name,
        col1 = c(
          paste(fixed_mean, value, sep = " + "),
          fixed_mean,
          paste(fixed_mean, value, sep = " + ")
        ),
        col2 = c(
          fixed_shape,
          paste(fixed_shape, value, sep = " + "),
          paste(fixed_shape, value, sep = " + ")
        ),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  return(out)
}

# Create the final data frame
design_matrix <- rbind(
  make_pairs(var_base, "base", fixed_mean, fixed_shape),
  make_pairs(var_mut, "mut", fixed_mean, fixed_shape),
  make_pairs(var_status, "status", fixed_mean, fixed_shape)
)
design_matrix[nrow(design_matrix) + 1, ] <- c("group", "offset(log_total_counted) + Exp + Temp + group",
                                              "Exp + Temp + p3_mut + mid_mut + p5_mut")
design_matrix[nrow(design_matrix) + 1, ] <- c("group", "offset(log_total_counted) + Exp + Temp + group",
                                              "Exp + Temp + p3_status + mid_status + p5_status")
# design_matrix[nrow(design_matrix) + 1, ] <- c("group", "offset(log_total_counted) + Exp + Temp + group",
#                                               "Exp + Temp + p3_base + mid_base + p5_base")
design_matrix[nrow(design_matrix) + 1, ] <- c("group", "offset(log_total_counted) + Exp + Temp + group",
                                              "Exp + Temp")
design_matrix$loo_code <- paste(design_matrix$vector, rownames(design_matrix), sep="")
replacements <- c(
  "offset(log_total_counted) + Exp + Temp + " = "",
  "Exp + Temp + " = "",
  "Exp + Temp" = "",
  "_base" = "Base",
  "_status" = "Status",
  "_mut" = "Mut",
  "p3" = "P3",
  "mid" = "Mid",
  "p5" = "P5",
  "offset(log_total_counted) + " = "",
  "group" = "Group"
)
design_matrix$loo_name <- paste("#Mean=", design_matrix$col1, "#Shape=", design_matrix$col2, sep="")
design_matrix$loo_name <- Reduce(
  function(text, pattern) {
    str_replace_all(
      text,
      fixed(pattern),
      replacements[[pattern]]
    )
  },
  names(replacements),
  init = design_matrix$loo_name
)
design_matrix$loo_name <- gsub("\\s+", "", design_matrix$loo_name)
design_matrix$group_loo <- NA
design_matrix$group_loo[c(1,4,7,10,13,16,19,22,25)] <- "Only mean"
design_matrix$group_loo[c(2,5,8,11,14,17,20,23,26)] <- "Only shape"
design_matrix$group_loo[c(3,6,9,12,15,18,21,24,27)] <- "Both"
design_matrix$group_loo[28] <- "Mean Group - Shape Mut"
design_matrix$group_loo[29] <- "Mean Group - Shape Status"
design_matrix$group_loo[30] <- "Mean Group - Shape"

################################################################################
# 3. MODELING
################################################################################
#PCR
# for (rows in 1:nrow(design_matrix)){
#   aux <- file.exists(paste("./results_final_PCR_DePCR/PCR/rds/fit/",
#                            design_matrix$vector[rows],
#                            rownames(design_matrix)[rows],
#                            ".rds",
#                            sep=""))
#   if (aux) {
#     fit <- readRDS(paste("./results_final_PCR_DePCR/PCR/rds/fit/",
#                          design_matrix$vector[rows],
#                          rownames(design_matrix)[rows],
#                          ".rds",
#                          sep=""))
#     print(paste("Loo for",
#                 design_matrix$vector[rows],
#                 rownames(design_matrix)[rows],
#                 sep=" "))
#     loo_fit <- loo(fit, save_psis = TRUE)
#     tiff(
#       filename = paste(
#         "./results_final_PCR_DePCR/PCR/rds/loo_psis/",
#         design_matrix$vector[rows],
#         rownames(design_matrix)[rows],
#         ".tiff",
#         sep = ""
#       ),
#       width = 7,
#       height = 5,
#       units = "in",
#       res = 400,
#       compression = "lzw"
#     )
#     plot(
#       loo_fit,
#       main = design_matrix$loo_name_clean[rows]
#     )
#     dev.off()
#   }
#   else {
#     print(paste("Creating design formula for",
#                 design_matrix$vector[rows],
#                 rownames(design_matrix)[rows],
#                 sep=" "))
#     design_formula <- bf(
#       as.formula(paste("Counts ~", design_matrix[rows, 2])),
#       as.formula(paste("shape ~", design_matrix[rows, 3]))
#     )
#     
#     print(paste("Adjusting model for",
#                 design_matrix$vector[rows],
#                 rownames(design_matrix)[rows],
#                 sep=" "))
#     fit <- brm(
#       formula = design_formula,
#       data = data_PCR,
#       family = negbinomial(link = "log", link_shape = "log"),
#       chains = 4,
#       cores = 4,
#       iter = 4000,
#       warmup = 1000,
#       control = list(adapt_delta = 0.99, max_treedepth = 15),
#       refresh = 100,
#       seed = 10231991
#     )
#     saveRDS(fit,
#             paste("./results_final_PCR_DePCR/PCR/rds/fit/",
#                   design_matrix$vector[rows],
#                   rownames(design_matrix)[rows],
#                   ".rds",
#                   sep=""))
#     
#     print(paste("Loo for",
#                 design_matrix$vector[rows],
#                 rownames(design_matrix)[rows],
#                 sep=" "))
#     loo_fit <- loo(fit)
#     saveRDS(loo_fit,
#             paste("./results_final_PCR_DePCR/PCR/rds/loo/",
#                   design_matrix$vector[rows],
#                   rownames(design_matrix)[rows],
#                   ".rds",
#                   sep=""))
#     loo_fit <- loo(fit, save_psis = TRUE)
#     tiff(
#       filename = paste(
#         "./results_final_PCR_DePCR/PCR/rds/loo_psis/",
#         design_matrix$vector[rows],
#         rownames(design_matrix)[rows],
#         ".tiff",
#         sep = ""
#       ),
#       width = 7,
#       height = 5,
#       units = "in",
#       res = 400,
#       compression = "lzw"
#     )
#     plot(
#       loo_fit,
#       main = design_matrix$loo_name_clean[rows]
#     )
#     dev.off()
#   }
# }
# 
# #DePCR
# for (rows in 1:nrow(design_matrix)){
#   aux <- file.exists(paste("./results_final_PCR_DePCR/DePCR/rds/fit/",
#                            design_matrix$vector[rows],
#                            rownames(design_matrix)[rows],
#                            ".rds",
#                            sep=""))
#   if (aux) {
#     fit <- readRDS(paste("./results_final_PCR_DePCR/DePCR/rds/fit/",
#                          design_matrix$vector[rows],
#                          rownames(design_matrix)[rows],
#                          ".rds",
#                          sep=""))
#     print(paste("Loo for",
#                 design_matrix$vector[rows],
#                 rownames(design_matrix)[rows],
#                 sep=" "))
#     loo_fit <- loo(fit, save_psis = TRUE)
#     tiff(
#       filename = paste(
#         "./results_final_PCR_DePCR/DePCR/rds/loo_psis/",
#         design_matrix$vector[rows],
#         rownames(design_matrix)[rows],
#         ".tiff",
#         sep = ""
#       ),
#       width = 7,
#       height = 5,
#       units = "in",
#       res = 400,
#       compression = "lzw"
#     )
#     
#     plot(
#       loo_fit,
#       main = design_matrix$loo_name[rows]
#     )
#     
#     dev.off()
#   }
#   else {
#     print(paste("Creating design formula for",
#                 design_matrix$vector[rows],
#                 rownames(design_matrix)[rows],
#                 sep=" "))
#     design_formula <- bf(
#       as.formula(paste("Counts ~", design_matrix[rows, 2])),
#       as.formula(paste("shape ~", design_matrix[rows, 3]))
#     )
#     
#     print(paste("Adjusting model for",
#                 design_matrix$vector[rows],
#                 rownames(design_matrix)[rows],
#                 sep=" "))
#     fit <- brm(
#       formula = design_formula,
#       data = data_DePCR,
#       family = negbinomial(link = "log", link_shape = "log"),
#       chains = 4,
#       cores = 4,
#       iter = 4000,
#       warmup = 1000,
#       control = list(adapt_delta = 0.99, max_treedepth = 15),
#       refresh = 100,
#       seed = 10231991
#     )
#     saveRDS(fit,
#             paste("./results_final_PCR_DePCR/DePCR/rds/fit/",
#                   design_matrix$vector[rows],
#                   rownames(design_matrix)[rows],
#                   ".rds",
#                   sep=""))
#     
#     print(paste("Loo for",
#                 design_matrix$vector[rows],
#                 rownames(design_matrix)[rows],
#                 sep=" "))
#     loo_fit <- loo(fit)
#     saveRDS(loo_fit,
#             paste("./results_final_PCR_DePCR/DePCR/rds/loo/",
#                   design_matrix$vector[rows],
#                   rownames(design_matrix)[rows],
#                   ".rds",
#                   sep=""))
#     loo_fit <- loo(fit, save_psis = TRUE)
#     tiff(
#       filename = paste(
#         "./results_final_PCR_DePCR/DePCR/rds/loo_psis/",
#         design_matrix$vector[rows],
#         rownames(design_matrix)[rows],
#         ".tiff",
#         sep = ""
#       ),
#       width = 7,
#       height = 5,
#       units = "in",
#       res = 400,
#       compression = "lzw"
#     )
#     
#     plot(
#       loo_fit,
#       main = design_matrix$loo_name[rows]
#     )
#     
#     dev.off()
#   }
# }

################################################################################
# 3. LOO COMPARE AND GENERAL GRAPHICS
################################################################################
#PCR
loo_files_PCR <- list.files("./results_final_PCR_DePCR/PCR/rds/loo/", pattern = "\\.rds$", full.names = TRUE)
loos_PCR <- setNames(lapply(loo_files_PCR, readRDS), tools::file_path_sans_ext(basename(loo_files_PCR)))
comp_loo_PCR <- loo_compare(loos_PCR)
comp_loo_PCR <- comp_loo_PCR %>% left_join(design_matrix, by = c("model" = "loo_code"))
comp_loo_PCR$group_plot <- comp_loo_PCR$group_loo
comp_loo_PCR$group_plot[grepl(comp_loo_PCR$group_plot, pattern="Group")] <- "Group"
comp_loo_PCR$group_plot <- factor(comp_loo_PCR$group_plot, levels = c("Group", "Both", "Only mean", "Only shape"))
comp_loo_PCR$loo_name2 <- comp_loo_PCR$loo_name
comp_loo_PCR$vector[1] <- "status"
comp_loo_PCR$vector[2] <- "mut"
comp_loo_PCR$vector[3] <- "none"
for (rows in 1:nrow(comp_loo_PCR)){
  if (grepl(comp_loo_PCR$loo_name2[rows],pattern="P5")){
    comp_loo_PCR$loo_name2[rows] <- "P3_Mid_P5"
  }
    else {
      if (grepl(comp_loo_PCR$loo_name2[rows],pattern="Mid")){
        comp_loo_PCR$loo_name2[rows] <- "P3_Mid"
      }
      else {
        if (grepl(comp_loo_PCR$loo_name2[rows],pattern="P3")) {
          comp_loo_PCR$loo_name2[rows] <- "P3"
        }
        else {
          comp_loo_PCR$loo_name2[rows] <- "P3_Mid_P5"
        }
      }
    }
  }

#GENERAL PLOT
ggplot(comp_loo_PCR, aes(x = elpd_diff, y = reorder(model, elpd_diff), shape = vector, color = loo_name2)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(xmin = elpd_diff - se_diff, xmax = elpd_diff + se_diff), orientation = "y", width = 0.15, linewidth = 0.5) +
  geom_point(size = 5) +
  facet_grid(rows = vars(group_plot), scales = "free_y", space = "free_y") +
  labs(x = expression(Delta*"ELPD vs. reference model"),y = NULL,
       shape = "Variable",title = "Predictive gain by model component - PCR dataset") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey95"),
        strip.text.y = element_text(face = "bold"),
        panel.spacing.y = unit(0.7, "lines"))

###############################
#DePCR
loo_files_DePCR <- list.files("./results_final_PCR_DePCR/DePCR/rds/loo/", pattern = "\\.rds$", full.names = TRUE)
loos_DePCR <- setNames(lapply(loo_files_DePCR, readRDS), tools::file_path_sans_ext(basename(loo_files_DePCR)))
comp_loo_DePCR <- loo_compare(loos_DePCR)
comp_loo_DePCR <- comp_loo_DePCR %>% left_join(design_matrix, by = c("model" = "loo_code"))
comp_loo_DePCR$group_plot <- comp_loo_DePCR$group_loo
comp_loo_DePCR$group_plot[grepl(comp_loo_DePCR$group_plot, pattern="Group")] <- "Group"
comp_loo_DePCR$group_plot <- factor(comp_loo_DePCR$group_plot, levels = c("Group", "Both", "Only mean", "Only shape"))
comp_loo_DePCR$loo_name2 <- comp_loo_DePCR$loo_name
comp_loo_DePCR$vector[1] <- "status"
comp_loo_DePCR$vector[2] <- "mut"
comp_loo_DePCR$vector[3] <- "none"
for (rows in 1:nrow(comp_loo_DePCR)){
  if (grepl(comp_loo_DePCR$loo_name2[rows],pattern="P5")){
    comp_loo_DePCR$loo_name2[rows] <- "P3_Mid_P5"
  }
  else {
    if (grepl(comp_loo_DePCR$loo_name2[rows],pattern="Mid")){
      comp_loo_DePCR$loo_name2[rows] <- "P3_Mid"
    }
    else {
      if (grepl(comp_loo_DePCR$loo_name2[rows],pattern="P3")) {
        comp_loo_DePCR$loo_name2[rows] <- "P3"
      }
      else {
        comp_loo_DePCR$loo_name2[rows] <- "P3_Mid_P5"
      }
    }
  }
}

#GENERAL PLOT
ggplot(comp_loo_DePCR, aes(x = elpd_diff, y = reorder(model, elpd_diff), shape = vector, color = loo_name2)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(xmin = elpd_diff - se_diff, xmax = elpd_diff + se_diff), orientation = "y", width = 0.15, linewidth = 0.5) +
  geom_point(size = 5) +
  facet_grid(rows = vars(group_plot), scales = "free_y", space = "free_y") +
  labs(x = expression(Delta*"ELPD vs. reference model"),y = NULL,
       shape = "Variable",title = "Predictive gain by model component - DePCR dataset") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey95"),
        strip.text.y = element_text(face = "bold"),
        panel.spacing.y = unit(0.7, "lines"))

# ################################################################################
# # 4. SPECIFIC PLOT
# ################################################################################
#PCR
cmp <- function(a, b) {
  d <- a$pointwise[, "elpd_loo"] - b$pointwise[, "elpd_loo"]
  tibble(elpd_diff = sum(d), se_diff = sqrt(length(d) * var(d)))
}

pairs_PCR <- comp_loo_PCR %>%
  filter(group_plot %in% c("Both", "Only mean", "Only shape")) %>%
  select(model, vector, loo_name2, group_plot) %>%
  pivot_wider(names_from = group_plot, values_from = model)

gain_results_PCR <- bind_rows(
  pairs_PCR %>% rowwise() %>% mutate(x = list(cmp(loos_PCR[[Both]], loos_PCR[[`Only mean`]]))) %>% ungroup() %>% unnest(x) %>% mutate(component_added = "Add to shape"),
  pairs_PCR %>% rowwise() %>% mutate(x = list(cmp(loos_PCR[[Both]], loos_PCR[[`Only shape`]]))) %>% ungroup() %>% unnest(x) %>% mutate(component_added = "Add to mean")
)

ggplot(gain_results_PCR, aes(elpd_diff, loo_name2, color = component_added)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbar(aes(xmin = elpd_diff - se_diff, xmax = elpd_diff + se_diff),
                orientation = "y", width = .15, position = position_dodge(.5)) +
  geom_point(size = 3, position = position_dodge(.5)) +
  facet_wrap(~vector, ncol = 1) +
  labs(
    x = expression(Delta*"ELPD gained - PCR dataset"),
    y = NULL,
    color = "Component added"
  ) +
  theme_bw()

###############################
#DePCR
pairs_DePCR <- comp_loo_DePCR %>%
  filter(group_plot %in% c("Both", "Only mean", "Only shape")) %>%
  select(model, vector, loo_name2, group_plot) %>%
  pivot_wider(names_from = group_plot, values_from = model)

gain_results_DePCR <- bind_rows(
  pairs_DePCR %>% rowwise() %>% mutate(x = list(cmp(loos_DePCR[[Both]], loos_DePCR[[`Only mean`]]))) %>% ungroup() %>% unnest(x) %>% mutate(component_added = "Add to shape"),
  pairs_DePCR %>% rowwise() %>% mutate(x = list(cmp(loos_DePCR[[Both]], loos_DePCR[[`Only shape`]]))) %>% ungroup() %>% unnest(x) %>% mutate(component_added = "Add to mean")
)

ggplot(gain_results_DePCR, aes(elpd_diff, loo_name2, color = component_added)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbar(aes(xmin = elpd_diff - se_diff, xmax = elpd_diff + se_diff),
                orientation = "y", width = .15, position = position_dodge(.5)) +
  geom_point(size = 3, position = position_dodge(.5)) +
  facet_wrap(~vector, ncol = 1) +
  labs(
    x = expression(Delta*"ELPD gained - DePCR dataset"),
    y = NULL,
    color = "Component added"
  ) +
  theme_bw()

###################################################
# BEST MODELS
# #PCR
design_formula <- bf(
  as.formula(paste("Counts ~", design_matrix[29, 2])),
  as.formula(paste("shape ~", design_matrix[29, 3]))
)
data_PCR$group <- relevel(factor(data_PCR$group), ref = "A_T_C")
data_PCR$group_status <- relevel(factor(data_PCR$group_status), ref = "Match_Match_Match")
fit_PCR <- brm(
  formula = design_formula,
  data = data_PCR,
  family = negbinomial(link = "log", link_shape = "log"),
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  refresh = 100,
  seed = 10231991
)
fit_PCR

#DePCR
data_DePCR$group <- relevel(factor(data_DePCR$group), ref = "A_T_C")
data_DePCR$group_status <- relevel(factor(data_DePCR$group_status), ref = "Match_Match_Match")

bf_disp_status <- bf(
  Counts ~ offset(log_total_counted) + group + Exp + Temp,
  shape ~ p3_status + mid_status + p5_status + Exp + Temp
)

fit_DePCR <- brm(
  formula = bf_disp_status,
  data = data_DePCR,
  family = negbinomial(link = "log", link_shape = "log"),
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  refresh = 100,
  seed = 10231991
)

fit_DePCR

################################################################################
# 4. DEFAULT PRIORS USED BY brms
################################################################################
prior_summary(fit_PCR)
prior_summary(fit_DePCR)

################################################################################
# 5. EXTERNAL VALIDATION: B10/B64
################################################################################
#PCR
data2 <- data
data2 <- data2 %>% filter(!Exp %in% c("A10", "A64"))

data2$Exp <- dplyr::recode(data2$Exp, "B10" = "A10")
data2$Exp <- dplyr::recode(data2$Exp, "A27" = "A64")
data2$Exp <- dplyr::recode(data2$Exp, "B27" = "A64")

# data2_PCR <- data2[data2$Method=="PCR",]
# data2_PCR <- data2_PCR %>% filter(Exp %in% c("A10"))
# 
# ep_PCR <- posterior_epred(fit_PCR, newdata = data2_PCR, ndraws = 2000)
# pp_PCR <- posterior_predict(fit_PCR, newdata = data2_PCR, ndraws = 2000)
# 
# pred_PCR <- data2_PCR %>%
#   mutate(
#     predicted = colMeans(ep_PCR),
#     ep_lower = apply(ep_PCR, 2, quantile, 0.025),
#     ep_upper = apply(ep_PCR, 2, quantile, 0.975),
#     pp_lower = apply(pp_PCR, 2, quantile, 0.025),
#     pp_upper = apply(pp_PCR, 2, quantile, 0.975)
#   )
# 
# pred_PCR <- pred_PCR %>%
#   mutate(template = sub("V[0-9]+$", "", id_feature)) %>%
#   group_by(template) %>%
#   mutate(raw_counts_percent  = 100 * Counts / sum(Counts), pred_counts_percent = 100 * predicted / sum(predicted)) %>%
#   ungroup()
# 
# ggplot(pred_PCR, aes(x = raw_counts_percent, y = pred_counts_percent, color=factor(group_status))) +
#   geom_abline(slope = 1, intercept = 0, linetype = 2) +
#   geom_point(size = 2, alpha = 0.7) +
#   facet_wrap(~template) +
#   labs(
#     x = "Observed proportion within template (%)",
#     y = "Predicted proportion within template (%)",
#     title = "PCR - Prediction of within-template primer distribution"
#   ) + theme_bw()

#Prediction for DePCR
data2_DePCR <- data2[data2$Method=="DePCR",]
data2_DePCR <- data2_DePCR %>% filter(Exp %in% c("A64"))
fit_PCR
fit_DePCR
ep_DePCR <- posterior_epred(fit_DePCR, newdata = data2_DePCR, ndraws = 2000)
pp_DePCR <- posterior_predict(fit_DePCR, newdata = data2_DePCR, ndraws = 2000)

pred_DePCR <- data2_DePCR %>%
  mutate(
    predicted = colMeans(ep_DePCR),
    ep_lower = apply(ep_DePCR, 2, quantile, 0.025),
    ep_upper = apply(ep_DePCR, 2, quantile, 0.975),
    pp_lower = apply(pp_DePCR, 2, quantile, 0.025),
    pp_upper = apply(pp_DePCR, 2, quantile, 0.975)
  )
pred_DePCR <- pred_DePCR %>%
  mutate(template = sub("V[0-9]+$", "", id_feature)) %>%
  group_by(template) %>%
  mutate(raw_counts_percent  = 100 * Counts / sum(Counts), pred_counts_percent = 100 * predicted / sum(predicted)) %>%
  ungroup()
pred_DePCR$inter_temp_status <- paste(pred_DePCR$Temp, pred_DePCR$group_status)

ggplot(pred_DePCR, aes(x = raw_counts_percent , y = pred_counts_percent, color=inter_temp_status)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~template) +
  labs(
    x = "Observed proportion within template (%)",
    y = "Predicted proportion within template (%)",
    title = "DePCR dataset with PCR model - Prediction of within-template primer distribution"
  ) +
  theme_bw()

################################################################################
# INVESTIGACAO
################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)

temp_group <- data_DePCR %>%
  mutate(rate = Counts / total_counts) %>%
  group_by(group, Temp) %>%
  summarise(
    mean_rate = mean(rate),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Temp,
    values_from = mean_rate
  ) %>%
  mutate(
    log_ratio = log(`55C` / `45C`)
  )
ggplot(temp_group,
       aes(x = reorder(group, log_ratio),
           y = log_ratio)) +
  geom_point() +
  geom_hline(
    yintercept = -0.94,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted"
  ) +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Primer group",
    y = "Observed log rate ratio: 55C / 45C - DePCR"
  )

data_test <- data[data$Exp %in% c("A10", "A27"),]
temp_effect <- data_test %>%
  mutate(rate = Counts / total_counts) %>%
  group_by(Method, group, Temp) %>%
  summarise(rate = mean(rate), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = Temp,
    values_from = rate
  ) %>%
  mutate(
    log_ratio = log(`55C` / `45C`)
  )
temp_effect %>%
  group_by(Method) %>%
  summarise(
    mean_log_ratio = mean(log_ratio),
    sd_log_ratio = sd(log_ratio),
    min_log_ratio = min(log_ratio),
    max_log_ratio = max(log_ratio)
  )

ggplot(temp_effect, aes(x = log_ratio, fill = Method)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw()

library(dplyr)
library(ggplot2)
library(tibble)
library(stringr)

################################################################################
# 1) INPUT DATA
################################################################################

# ---------------------------
# PCR
# ---------------------------
coef_pcr <- tribble(
  ~Block, ~Term, ~Estimate, ~Lower, ~Upper,
  
  # Mean - global
  "Mean - global effects", "ExpA64",    -1.15, -1.17, -1.14,
  "Mean - global effects", "Temp55C",   -0.01, -0.02,  0.00,
  
  # Shape - global
  "Shape - global effects", "shape_ExpA64",   -2.55, -3.70, -1.65,
  "Shape - global effects", "shape_Temp55C",   0.20, -0.63,  1.05,
  
  # Shape - mismatch-position
  "Shape - mismatch-position effects", "shape_p3_statusMismatch",   0.22, -0.43,  0.89,
  "Shape - mismatch-position effects", "shape_mid_statusMismatch",  1.53,  0.88,  2.26,
  "Shape - mismatch-position effects", "shape_p5_statusMismatch",   1.75,  1.03,  2.62,
  
  # Mean - group
  "Mean - primer-group effects", "groupA_A_A", -1.21, -1.25, -1.16,
  "Mean - primer-group effects", "groupA_A_C", -0.30, -0.33, -0.28,
  "Mean - primer-group effects", "groupA_A_G", -0.95, -1.00, -0.91,
  "Mean - primer-group effects", "groupA_A_T", -1.10, -1.14, -1.05,
  "Mean - primer-group effects", "groupA_C_A", -1.15, -1.19, -1.11,
  "Mean - primer-group effects", "groupA_C_C", -0.13, -0.15, -0.11,
  "Mean - primer-group effects", "groupA_C_G", -1.20, -1.24, -1.16,
  "Mean - primer-group effects", "groupA_C_T", -1.19, -1.23, -1.14,
  "Mean - primer-group effects", "groupA_G_A", -0.70, -0.74, -0.66,
  "Mean - primer-group effects", "groupA_G_C", -0.11, -0.13, -0.09,
  "Mean - primer-group effects", "groupA_G_G", -0.99, -1.03, -0.95,
  "Mean - primer-group effects", "groupA_G_T", -1.27, -1.32, -1.23,
  "Mean - primer-group effects", "groupA_T_A", -0.15, -0.17, -0.13,
  "Mean - primer-group effects", "groupA_T_G", -0.03, -0.05, -0.01,
  "Mean - primer-group effects", "groupA_T_T", -0.28, -0.30, -0.26,
  "Mean - primer-group effects", "groupC_A_A", -0.65, -0.69, -0.61,
  "Mean - primer-group effects", "groupC_A_C", -0.78, -0.83, -0.74,
  "Mean - primer-group effects", "groupC_A_G", -1.18, -1.22, -1.14,
  "Mean - primer-group effects", "groupC_A_T", -1.38, -1.43, -1.33,
  "Mean - primer-group effects", "groupC_C_A", -0.89, -0.93, -0.85,
  "Mean - primer-group effects", "groupC_C_C", -1.06, -1.11, -1.01,
  "Mean - primer-group effects", "groupC_C_G", -1.15, -1.19, -1.11,
  "Mean - primer-group effects", "groupC_C_T", -1.69, -1.74, -1.64,
  "Mean - primer-group effects", "groupC_G_A", -0.62, -0.65, -0.58,
  "Mean - primer-group effects", "groupC_G_C", -1.14, -1.18, -1.09,
  "Mean - primer-group effects", "groupC_G_G", -0.81, -0.85, -0.78,
  "Mean - primer-group effects", "groupC_G_T", -1.22, -1.27, -1.18,
  "Mean - primer-group effects", "groupC_T_A", -1.05, -1.09, -1.00,
  "Mean - primer-group effects", "groupC_T_C", -0.34, -0.37, -0.32,
  "Mean - primer-group effects", "groupC_T_G", -1.33, -1.38, -1.28,
  "Mean - primer-group effects", "groupC_T_T", -1.21, -1.26, -1.17,
  "Mean - primer-group effects", "groupG_A_A", -0.79, -0.83, -0.75,
  "Mean - primer-group effects", "groupG_A_C", -1.42, -1.47, -1.37,
  "Mean - primer-group effects", "groupG_A_G", -1.66, -1.70, -1.61,
  "Mean - primer-group effects", "groupG_A_T", -1.85, -1.90, -1.79,
  "Mean - primer-group effects", "groupG_C_A", -1.61, -1.66, -1.56,
  "Mean - primer-group effects", "groupG_C_C", -1.41, -1.46, -1.37,
  "Mean - primer-group effects", "groupG_C_G", -1.69, -1.74, -1.64,
  "Mean - primer-group effects", "groupG_C_T", -1.95, -2.01, -1.90,
  "Mean - primer-group effects", "groupG_G_A", -0.75, -0.79, -0.71,
  "Mean - primer-group effects", "groupG_G_C", -1.13, -1.17, -1.08,
  "Mean - primer-group effects", "groupG_G_G", -1.05, -1.09, -1.01,
  "Mean - primer-group effects", "groupG_G_T", -0.84, -0.88, -0.80,
  "Mean - primer-group effects", "groupG_T_A", -1.03, -1.08, -0.99,
  "Mean - primer-group effects", "groupG_T_C", -0.87, -0.90, -0.85,
  "Mean - primer-group effects", "groupG_T_G", -1.69, -1.74, -1.64,
  "Mean - primer-group effects", "groupG_T_T", -1.90, -1.96, -1.84,
  "Mean - primer-group effects", "groupT_A_A", -0.86, -0.90, -0.82,
  "Mean - primer-group effects", "groupT_A_C", -1.26, -1.31, -1.21,
  "Mean - primer-group effects", "groupT_A_G", -1.34, -1.39, -1.30,
  "Mean - primer-group effects", "groupT_A_T", -1.46, -1.51, -1.41,
  "Mean - primer-group effects", "groupT_C_A", -1.47, -1.52, -1.42,
  "Mean - primer-group effects", "groupT_C_C", -1.41, -1.46, -1.36,
  "Mean - primer-group effects", "groupT_C_G", -1.16, -1.21, -1.12,
  "Mean - primer-group effects", "groupT_C_T", -1.76, -1.82, -1.71,
  "Mean - primer-group effects", "groupT_G_A", -0.62, -0.65, -0.58,
  "Mean - primer-group effects", "groupT_G_C", -0.81, -0.85, -0.76,
  "Mean - primer-group effects", "groupT_G_G", -0.75, -0.79, -0.71,
  "Mean - primer-group effects", "groupT_G_T", -0.99, -1.03, -0.95,
  "Mean - primer-group effects", "groupT_T_A", -1.17, -1.22, -1.13,
  "Mean - primer-group effects", "groupT_T_C", -0.50, -0.53, -0.48,
  "Mean - primer-group effects", "groupT_T_G", -1.55, -1.60, -1.49,
  "Mean - primer-group effects", "groupT_T_T", -1.74, -1.79, -1.68
)

# ---------------------------
# DePCR
# ---------------------------
fit_DePCR
coef_depcr <- tribble(
  ~Block, ~Term, ~Estimate, ~Lower, ~Upper,
  
  # Mean - global
  "Mean - global effects", "ExpA64",  -0.15, -0.16, -0.14,
  "Mean - global effects", "Temp55C", -0.96, -1.03, -0.88,
  
  # Shape - global
  "Shape - global effects", "shape_ExpA64",   0.16, -0.27,  0.60,
  "Shape - global effects", "shape_Temp55C",  6.75,  6.23,  7.32,
  
  # Shape - mismatch-position effects
  "Shape - mismatch-position effects", "shape_p3_statusMismatch",   0.46,  0.11,  0.82,
  "Shape - mismatch-position effects", "shape_mid_statusMismatch",  0.54,  0.20,  0.88,
  "Shape - mismatch-position effects", "shape_p5_statusMismatch",   0.71,  0.37,  1.05,
  
  # Mean - primer-group effects
  "Mean - primer-group effects", "groupA_A_A", -3.89, -3.96, -3.81,
  "Mean - primer-group effects", "groupA_A_C", -1.25, -1.28, -1.23,
  "Mean - primer-group effects", "groupA_A_G", -3.82, -3.89, -3.75,
  "Mean - primer-group effects", "groupA_A_T", -4.25, -4.33, -4.16,
  "Mean - primer-group effects", "groupA_C_A", -4.05, -4.13, -3.97,
  "Mean - primer-group effects", "groupA_C_C", -1.41, -1.43, -1.38,
  "Mean - primer-group effects", "groupA_C_G", -4.37, -4.46, -4.28,
  "Mean - primer-group effects", "groupA_C_T", -4.66, -4.76, -4.55,
  "Mean - primer-group effects", "groupA_G_A", -2.01, -2.04, -1.97,
  "Mean - primer-group effects", "groupA_G_C", -0.33, -0.35, -0.30,
  "Mean - primer-group effects", "groupA_G_G", -2.37, -2.41, -2.32,
  "Mean - primer-group effects", "groupA_G_T", -3.10, -3.15, -3.05,
  "Mean - primer-group effects", "groupA_T_A", -0.52, -0.55, -0.50,
  "Mean - primer-group effects", "groupA_T_G", -0.61, -0.63, -0.58,
  "Mean - primer-group effects", "groupA_T_T", -1.07, -1.10, -1.05,
  "Mean - primer-group effects", "groupC_A_A", -5.02, -5.15, -4.90,
  "Mean - primer-group effects", "groupC_A_C", -4.10, -4.18, -4.02,
  "Mean - primer-group effects", "groupC_A_G", -5.64, -5.81, -5.48,
  "Mean - primer-group effects", "groupC_A_T", -6.30, -6.53, -6.09,
  "Mean - primer-group effects", "groupC_C_A", -5.42, -5.57, -5.27,
  "Mean - primer-group effects", "groupC_C_C", -4.69, -4.80, -4.59,
  "Mean - primer-group effects", "groupC_C_G", -5.72, -5.90, -5.55,
  "Mean - primer-group effects", "groupC_C_T", -6.50, -6.75, -6.25,
  "Mean - primer-group effects", "groupC_G_A", -4.01, -4.09, -3.93,
  "Mean - primer-group effects", "groupC_G_C", -3.27, -3.33, -3.21,
  "Mean - primer-group effects", "groupC_G_G", -4.44, -4.53, -4.35,
  "Mean - primer-group effects", "groupC_G_T", -5.05, -5.17, -4.93,
  "Mean - primer-group effects", "groupC_T_A", -3.78, -3.85, -3.71,
  "Mean - primer-group effects", "groupC_T_C", -1.20, -1.23, -1.17,
  "Mean - primer-group effects", "groupC_T_G", -4.09, -4.17, -4.01,
  "Mean - primer-group effects", "groupC_T_T", -4.49, -4.59, -4.40,
  "Mean - primer-group effects", "groupG_A_A", -5.15, -5.29, -5.03,
  "Mean - primer-group effects", "groupG_A_C", -4.77, -4.88, -4.66,
  "Mean - primer-group effects", "groupG_A_G", -6.49, -6.73, -6.26,
  "Mean - primer-group effects", "groupG_A_T", -6.50, -6.74, -6.27,
  "Mean - primer-group effects", "groupG_C_A", -6.04, -6.24, -5.85,
  "Mean - primer-group effects", "groupG_C_C", -5.15, -5.29, -5.02,
  "Mean - primer-group effects", "groupG_C_G", -6.56, -6.82, -6.32,
  "Mean - primer-group effects", "groupG_C_T", -6.88, -7.17, -6.60,
  "Mean - primer-group effects", "groupG_G_A", -4.18, -4.27, -4.10,
  "Mean - primer-group effects", "groupG_G_C", -2.94, -2.99, -2.89,
  "Mean - primer-group effects", "groupG_G_G", -4.74, -4.85, -4.63,
  "Mean - primer-group effects", "groupG_G_T", -4.77, -4.88, -4.66,
  "Mean - primer-group effects", "groupG_T_A", -3.64, -3.71, -3.58,
  "Mean - primer-group effects", "groupG_T_C", -1.62, -1.65, -1.59,
  "Mean - primer-group effects", "groupG_T_G", -4.48, -4.57, -4.38,
  "Mean - primer-group effects", "groupG_T_T", -5.14, -5.28, -5.01,
  "Mean - primer-group effects", "groupT_A_A", -5.37, -5.52, -5.23,
  "Mean - primer-group effects", "groupT_A_C", -4.80, -4.91, -4.69,
  "Mean - primer-group effects", "groupT_A_G", -5.77, -5.95, -5.59,
  "Mean - primer-group effects", "groupT_A_T", -6.22, -6.44, -6.01,
  "Mean - primer-group effects", "groupT_C_A", -6.01, -6.21, -5.82,
  "Mean - primer-group effects", "groupT_C_C", -5.16, -5.30, -5.03,
  "Mean - primer-group effects", "groupT_C_G", -5.89, -6.08, -5.71,
  "Mean - primer-group effects", "groupT_C_T", -6.64, -6.91, -6.39,
  "Mean - primer-group effects", "groupT_G_A", -4.25, -4.34, -4.17,
  "Mean - primer-group effects", "groupT_G_C", -3.38, -3.44, -3.32,
  "Mean - primer-group effects", "groupT_G_G", -4.60, -4.71, -4.50,
  "Mean - primer-group effects", "groupT_G_T", -5.08, -5.21, -4.96,
  "Mean - primer-group effects", "groupT_T_A", -4.07, -4.15, -3.99,
  "Mean - primer-group effects", "groupT_T_C", -1.69, -1.72, -1.67,
  "Mean - primer-group effects", "groupT_T_G", -4.46, -4.56, -4.37,
  "Mean - primer-group effects", "groupT_T_T", -5.10, -5.23, -4.97
)
################################################################################
# 2) FUNCTION TO CLEAN LABELS
################################################################################

clean_term <- function(x) {
  x %>%
    str_replace("^group", "") %>%
    str_replace_all("_", " ") %>%
    str_replace("^shape_", "") %>%
    str_replace("ExpA64", "Experiment: A64") %>%
    str_replace("Temp55C", "Temperature: 55C") %>%
    str_replace("p3_statusMismatch", "3' mismatch") %>%
    str_replace("mid_statusMismatch", "Middle mismatch") %>%
    str_replace("p5_statusMismatch", "5' mismatch")
}

################################################################################
# 3) PLOTTING FUNCTION
################################################################################

plot_coef_blocks <- function(df, title_text, file_name, height = 18) {
  
  block_levels <- c(
    "Mean - global effects",
    "Mean - primer-group effects",
    "Shape - global effects",
    "Shape - mismatch-position effects"
  )
  
  df2 <- df %>%
    mutate(
      Block = factor(Block, levels = block_levels),
      Term_clean = clean_term(Term),
      sig = ifelse(Lower <= 0 & Upper >= 0, "CI includes 0", "CI excludes 0")
    ) %>%
    group_by(Block) %>%
    arrange(Estimate, .by_group = TRUE) %>%
    mutate(Term_clean = factor(Term_clean, levels = unique(Term_clean))) %>%
    ungroup()
  
  p <- ggplot(df2, aes(x = Estimate, y = Term_clean, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_errorbar(
      aes(xmin = Lower, xmax = Upper),
      orientation = "y",
      width = 0.15,
      linewidth = 0.65
    ) +
    geom_point(size = 2.2) +
    facet_grid(Block ~ ., scales = "free_y", space = "free_y") +
    scale_color_manual(
      values = c("CI excludes 0" = "red3", "CI includes 0" = "grey60")
    ) +
    labs(
      x = "Posterior coefficient (log scale)",
      y = NULL,
      color = NULL,
      title = title_text
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      strip.text.y = element_text(face = "bold", angle = 0),
      legend.position = "top",
      panel.spacing.y = unit(0.8, "lines")
    )
  
  ggsave(file_name, p, width = 11, height = height, dpi = 400)
  print(p)
}

################################################################################
# 4) SAVE FIGURES
################################################################################

plot_coef_blocks(
  coef_pcr,
  title_text = "Posterior coefficient estimates - final PCR model",
  file_name = "coef_blocks_PCR.png",
  height = 20
)

plot_coef_blocks(
  coef_depcr,
  title_text = "Posterior coefficient estimates - final DePCR model",
  file_name = "coef_blocks_DePCR.png",
  height = 20
)

names(data)


################################################################################
library(dplyr)
library(ggplot2)

feature_counts_plot <- data %>%
  filter(total_counts > 13000) %>%
  mutate(
    Exp = factor(Exp),
    Temp = factor(Temp),
    Method = factor(Method)
  )
feature_counts_plot <- feature_counts_plot[feature_counts_plot$Exp %in% c("A10","A27","A64","B10","B27"),]
p_feature_counts <- ggplot(
  feature_counts_plot,
  aes(
    x = reorder(feature, Counts, FUN = median),
    y = Counts
  )
) +
  geom_boxplot(
    outlier.alpha = 0.25,
    linewidth = 0.35
  ) +
  scale_y_log10() +
  facet_grid(
    Method + Temp ~ Exp,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    x = "Primer-template feature",
    y = "Observed count (log10 scale)",
    title = "Distribution of primer-template counts across experimental conditions"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6
    ),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.7, "lines")
  )

ggsave(
  "feature_counts_by_experiment.png",
  p_feature_counts,
  width = 14,
  height = 9,
  dpi = 400
)

p_feature_counts


p_feature_counts_main <- data %>%
  filter(
    total_counts > 13000,
    Exp %in% c("A10", "A64")
  ) %>%
  ggplot(
    aes(
      x = reorder(feature, Counts, FUN = median),
      y = Counts
    )
  ) +
  geom_boxplot(
    outlier.alpha = 0.25,
    linewidth = 0.35
  ) +
  scale_y_log10() +
  facet_grid(
    Method + Temp ~ Exp,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    x = "Primer-template feature",
    y = "Observed count (log10 scale)",
    title = "Primer-template count heterogeneity in the modeling datasets"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6
    ),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  "feature_counts_A10_A64.png",
  p_feature_counts_main,
  width = 12,
  height = 8,
  dpi = 400
)

template_comp <- data %>%
  filter(
    total_counts > 13000,
    grepl("^B", Exp),
    !is.na(ST)
  ) %>%
  group_by(
    Exp,
    Method,
    Temp,
    Replicate,
    ST
  ) %>%
  summarise(
    template_counts = sum(Counts),
    .groups = "drop"
  ) %>%
  group_by(
    Exp,
    Method,
    Temp,
    Replicate
  ) %>%
  mutate(
    template_percent =
      100 * template_counts / sum(template_counts)
  ) %>%
  ungroup()

p_template_comp <- ggplot(
  template_comp,
  aes(
    x = ST,
    y = template_percent
  )
) +
  geom_hline(
    yintercept = 10,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_boxplot(
    outlier.alpha = 0.3,
    linewidth = 0.4
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.45,
    size = 1.2
  ) +
  facet_grid(
    Method + Temp ~ Exp
  ) +
  labs(
    x = "Synthetic template",
    y = "Observed template proportion (%)",
    title = "Observed template composition in equimolar B-series mixtures",
    subtitle = "Dashed line indicates the 10% proportion expected from an equimolar ten-template input"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.7, "lines")
  )

ggsave(
  "template_composition_equimolar.png",
  p_template_comp,
  width = 13,
  height = 9,
  dpi = 400
)

p_template_comp

template_distortion <- template_comp %>%
  group_by(
    Exp,
    Method,
    Temp,
    Replicate
  ) %>%
  summarise(
    mean_percent = mean(template_percent),
    sd_percent = sd(template_percent),
    min_percent = min(template_percent),
    max_percent = max(template_percent),
    range_percent = max(template_percent) - min(template_percent),
    deviation_from_equimolar =
      sum(abs(template_percent - 10)),
    .groups = "drop"
  )

template_distortion_summary <- template_distortion %>%
  group_by(
    Exp,
    Method,
    Temp
  ) %>%
  summarise(
    mean_sd = mean(sd_percent),
    mean_range = mean(range_percent),
    mean_deviation = mean(deviation_from_equimolar),
    .groups = "drop"
  )

template_distortion_summary
fit_DePCR
fit_PCR
