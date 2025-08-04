library(dplyr)
library(tidyr)
source("functions.R")

set.seed(123)
source_sample_sizes <- c(100, 1000, 5000, 10000)
target_sample_sizes <- c(100, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000)
num_simulations <- 100
degree <- 2
d <- 5
alpha <- 2.5
target_sample_size_fixed <- 1000

for (x0_label in c("targetonly", "interior")) {
  x0 <- if (x0_label == "targetonly") runif(d, 0.7, 0.8) else runif(d, 0.2, 0.3)
  print(x0)
  # Part 1: Line plot data
  all_results <- list()
  for (n_s in source_sample_sizes) {
    results <- calculate_mse_fixed_nP(
      target_sample_sizes, n_s, x0, alpha,
      num_simulations, degree,
      kernel_type = "gaussian"
    )
    results$source_sample_size <- n_s
    all_results[[as.character(n_s)]] <- results
  }
  final_results_line <- bind_rows(all_results)
  saveRDS(final_results_line, paste0("mse_simulation_", x0_label, ".rds"))
  
  # Part 2: Boxplot data
  all_results <- list()
  for (n_s in source_sample_sizes) {
    for (sim in 1:num_simulations) {
      x_s <- matrix(runif(n_s * d, min = -0.5, max = 0.5), ncol = d)
      y_s <- f(x_s, x0, alpha) + rnorm(n_s)
      x_t <- matrix(runif(target_sample_size_fixed * d, min = 0, max = 1), ncol = d)
      y_t <- f(x_t, x0, alpha) + rnorm(target_sample_size_fixed)
      
      n_phi_combined <- n_s * ifelse(all(x0 >= -0.5 & x0 <= 0.5), 1, 0) +
        target_sample_size_fixed * ifelse(all(x0 >= 0 & x0 <= 1), 1, 0)
      n_phi_target <- target_sample_size_fixed * ifelse(all(x0 >= 0 & x0 <= 1), 1, 0)
      
      h_combined <- n_phi_combined^(-1 / (2 * alpha + d))
      h_target <- n_phi_target^(-1 / (2 * alpha + d))
      
      x_combined <- rbind(x_s, x_t)
      y_combined <- c(y_s, y_t)
      
      m_hat_combined <- local_polynomial_estimator(x_combined, y_combined, x0, degree, h_combined)
      m_hat_target <- local_polynomial_estimator(x_t, y_t, x0, degree, h_target)
      f_true <- f(matrix(x0, ncol = d, byrow = TRUE), x0, alpha)
      
      all_results[[length(all_results) + 1]] <- data.frame(
        target_sample_size = target_sample_size_fixed,
        source_sample_size = n_s,
        estimator = "Proposed",
        mse = (m_hat_combined - f_true)^2
      )
      all_results[[length(all_results) + 1]] <- data.frame(
        target_sample_size = target_sample_size_fixed,
        source_sample_size = n_s,
        estimator = "Target-Only",
        mse = (m_hat_target - f_true)^2
      )
    }
  }
  final_results_box <- bind_rows(all_results)
  saveRDS(final_results_box, paste0("box_simulation_", x0_label, ".rds"))
}



library(dplyr)
library(tidyr)
library(ggplot2)

source_sample_sizes <- c(100, 1000, 5000, 10000)

for (x0_label in c("targetonly", "interior")) {
  # Line plot
  final_results <- readRDS(paste0("mse_simulation_", x0_label, ".rds")) %>%
    distinct(target_sample_size, estimator, source_sample_size, .keep_all = TRUE)
  
  final_results$estimator <- factor(tolower(trimws(final_results$estimator)),
                                    levels = c("combined", "target-only"))
  final_results <- final_results %>%
    filter(mse > 0 & lower_bound > 0 & upper_bound > 0) %>%
    mutate(
      mse = ifelse(mse <= 0, 1e-8, mse),
      lower_bound = ifelse(lower_bound <= 0, 1e-8, lower_bound),
      upper_bound = ifelse(upper_bound <= 0, 1e-8, upper_bound)
    )
  final_results$estimator <- recode(final_results$estimator,
                                    "combined" = "Proposed",
                                    "target-only" = "Target only")
  
  ggplot(final_results, aes(x = target_sample_size, y = mse, color = estimator, linetype = estimator, group = estimator)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.1, alpha = 0.7) +
    geom_ribbon(aes(ymin = mse - 1.96 * se, ymax = mse + 1.96 * se, fill = estimator),
                alpha = 0.2, linetype = 0, show.legend = TRUE) +
    scale_x_log10(labels = function(x) paste0(x / 1000)) +
    scale_y_log10() +
    geom_text(
      aes(label = target_sample_size),
      vjust = 0.1, hjust = 0.5,
      size = 5, check_overlap = TRUE
    ) +
    labs(
      title = "MSE Comparison: Proposed vs Target only Estimators",
      x = expression("Target Sample Size" ~ (10^3)),
      y = "Mean Squared Error (log scale)",
      color = "Estimator",
      linetype = "Estimator",
      fill = "95% Confidence Interval"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 18, angle = 30, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 18),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      strip.text = element_text(size = 20, face = "bold"),
      panel.spacing = unit(0.8, "lines"),
      legend.position = "top"
    ) +
    scale_color_manual(values = c("Proposed" = "#5f97d2", "Target only" = "#d76364")) +
    scale_linetype_manual(values = c("Proposed" = "solid", "Target only" = "dashed")) +
    scale_fill_manual(values = c("Proposed" = "#3f74a5", "Target only" = "#a84343")) +
    facet_wrap(~ source_sample_size, ncol = 2, labeller = labeller(source_sample_size = label_both))
  
  ggsave(paste0("mse_plot_", x0_label, ".png"), width = 10, height = 8, dpi = 300, bg = "white")
  
  # Boxplot
  final_results <- readRDS(paste0("box_simulation_", x0_label, ".rds"))
  final_results$estimator <- factor(tolower(trimws(final_results$estimator)),
                                    levels = c("proposed", "target-only"))
  final_results$estimator <- recode(final_results$estimator,
                                    "proposed" = "Proposed",
                                    "target-only" = "Target only")
  
  ggplot(final_results, aes(x = as.factor(source_sample_size), y = mse, fill = estimator)) +
    geom_boxplot(outlier.shape = NA) +
    scale_y_log10() +
    scale_fill_manual(values = c("Proposed" = "#5f97d2", "Target only" = "#d76364")) +
    labs(
      title = paste0("Boxplot of MSE (n_Q = 1000)"),
      x = "Source Sample Size",
      y = "Mean Squared Error (log scale)",
      fill = "Estimator"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 18, angle = 30, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 18),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    )
  
  ggsave(paste0("box_plot_", x0_label, ".png"), width = 10, height = 7, dpi = 300, bg = "white")
}




