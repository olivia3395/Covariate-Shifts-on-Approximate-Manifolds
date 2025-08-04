library(dplyr)
library(tidyr)
library(ggplot2)
source("functions.R")  # assumes you have defined local_polynomial_estimator, f(), calculate_mse_fixed_nP, and calculate_mse_fixed_nP_manifold

# Dimensions and manifold mapping
D <- 5     # Ambient space dimension
d <- 2     # Intrinsic manifold dimension
alpha <- 2.5
degree <- 2

# Smooth nonlinear mapping from manifold to ambient space
phi <- function(z) {
  cbind(
    (z + 1) / 2,
    z[, 1]^2,
    z[, 2]^2,
    (z[, 1] + 1) * (z[, 2] + 1) / 4
  )
}

# Estimation point on the manifold
set.seed(123)
# z_manifold <- runif(d, -1, 1)
# x0_manifold <- phi(matrix(z_manifold, nrow = 1))


# -----------------------------------------------
# Target Function f(x)
# -----------------------------------------------
f <- function(x, x0, alpha) {
  rowSums(abs(x - matrix(x0, nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^alpha)
}


# Target point
z0 <- runif(d, -1, 1)
x0_manifold <- phi(matrix(z0, nrow = 1))


# Simulation setup
source_sample_sizes <- c(100, 1000, 5000, 10000)
target_sample_sizes <- c(100, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000)
num_simulations <- 100
target_sample_size_fixed <- 1000

# -------------------------------
# Part 1: Line Plot - MSE vs n_Q
# -------------------------------
all_results <- list()
for (n_s in source_sample_sizes) {
  results <- calculate_mse_fixed_nP_manifold(
    target_sample_sizes, n_s, x0_manifold, alpha,
    num_simulations, degree, kernel_type = "gaussian"
  )
  results$source_sample_size <- n_s
  all_results[[as.character(n_s)]] <- results
}
final_results <- bind_rows(all_results) %>%
  distinct(target_sample_size, estimator, source_sample_size, .keep_all = TRUE)

final_results$estimator <- factor(tolower(trimws(final_results$estimator)),
                                  levels = c("combined", "target-only"))
final_results$estimator <- recode(final_results$estimator,
                                  "combined" = "Combined",
                                  "target-only" = "Target only")

final_results <- final_results %>%
  filter(mse > 0 & lower_bound > 0 & upper_bound > 0) %>%
  mutate(
    mse = ifelse(mse <= 0, 1e-8, mse),
    lower_bound = ifelse(lower_bound <= 0, 1e-8, lower_bound),
    upper_bound = ifelse(upper_bound <= 0, 1e-8, upper_bound)
  )

# Plot MSE curve
ggplot(final_results, aes(x = target_sample_size, y = mse, color = estimator, linetype = estimator, group = estimator)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.1, alpha = 0.7) +
  geom_ribbon(aes(ymin = mse - 1.96 * se, ymax = mse + 1.96 * se, fill = estimator),
              alpha = 0.2, linetype = 0, show.legend = TRUE) +
  scale_x_log10(labels = function(x) paste0(x / 1000)) +
  scale_y_log10() +
  geom_text(aes(label = target_sample_size), vjust = -0.8, size = 3, check_overlap = TRUE) +
  labs(
    title = "MSE Comparison on Manifold: Combined vs Target only",
    x = expression("Target Sample Size" ~ (10^3)),
    y = "Mean Squared Error (log scale)",
    color = "Estimator", linetype = "Estimator", fill = "95% Confidence Interval"
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
  )+
  scale_color_manual(values = c("Combined" = "#5f97d2", "Target only" = "#d76364")) +
  scale_linetype_manual(values = c("Combined" = "solid", "Target only" = "dashed")) +
  scale_fill_manual(values = c("Combined" = "#3f74a5", "Target only" = "#a84343")) +
  facet_wrap(~ source_sample_size, ncol = 2, labeller = labeller(source_sample_size = label_both))

ggsave("mse_plot_manifold.png", width = 10, height = 8, dpi = 300, bg = "white")


# -------------------------------
# Part 2: Boxplot - Fixed n_Q, vary n_P
# -------------------------------
all_results <- list()
for (n_s in source_sample_sizes) {
  for (sim in 1:num_simulations) {
   
    x_s <- matrix(runif(n_s * D, 0, 1), ncol = D)
    
    y_s <- f(x_s, x0_manifold, alpha) + rnorm(n_s)
      
    # Target covariates
    z_t <- matrix(runif(target_sample_size_fixed * d, -1, 1), ncol = d)
    rho_n <- 0
    u_t <- matrix(runif(target_sample_size_fixed * D, -1, 1), ncol = D)
    x_t <- phi(z_t) + rho_n * u_t
    
    
    y_t <- f(x_t, x0_manifold, alpha) + rnorm(target_sample_size_fixed)
    
    
    
    
    h_combined <- (target_sample_size_fixed + n_s^((2 * alpha + d) / (2 * alpha + D)))^(-1 / (2 * alpha + d))
   
    
    
    
   
    h_target <- target_sample_size_fixed^(-1 / (2 * alpha + d))
    
    x_combined <- rbind(x_s, x_t)
    y_combined <- c(y_s, y_t)
    
    m_hat_combined <- local_polynomial_estimator(x_combined, y_combined, x0_manifold, degree, h_combined)
    m_hat_target <- local_polynomial_estimator(x_t, y_t, x0_manifold, degree, h_target)
    f_true <- f(matrix(x0_manifold, ncol = D, byrow = TRUE), x0_manifold, alpha)
    
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

final_results <- bind_rows(all_results)
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
    title = paste0("Boxplot of MSE on Manifold (n_Q = ", target_sample_size_fixed, ")"),
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

ggsave("box_plot_manifold.png", width = 10, height = 7, dpi = 300, bg = "white")





