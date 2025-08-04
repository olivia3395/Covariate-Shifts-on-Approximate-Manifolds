# simulation_adaptive.R
# -----------------------------------------------------------
# For each target sample size, estimate intrinsic dimension d
# from target data, then perform Lepskii-based MSE simulation.
# -----------------------------------------------------------



library(dplyr)
library(tidyr)
library(ggplot2)
library(FNN)
library(scales)

library(parallel)


rm(list = ls())
setwd("/Users/yuyao/Desktop/JASA-CODE")

# Load estimator function
source("functions.R") 







# ----------------------------------------------------
# Intrinsic Dimension Estimation (Levina-Bickel Method)
# ----------------------------------------------------

estimate_intrinsic_dimension_manifold <- function(X, k = 10) {
  # Compute k-nearest neighbor distances
  nn_distances <- FNN::get.knn(X, k = k)$nn.dist
  r_k <- nn_distances[, k]
  r_half_k <- nn_distances[, ceiling(k / 2)]
  
  # Local dimension estimate
  local_d <- log(2) / log(r_k / r_half_k)
  local_d <- local_d[is.finite(local_d) & local_d > 0]
  
  # Average-based estimator
  d_avg <- round(mean(pmin(local_d, D)))
  
  # Voting-based estimator
  d_freq_table <- table(round(local_d))
  d_vote <- as.numeric(names(d_freq_table)[which.max(d_freq_table)])
  
  return(as.numeric(d_avg))  # or return both d_avg and d_vote if needed
}



# --- Fixed Source Size on Manifold with Lepskii's Adaptive Procedure ---

# This function performs a simulation study to compare adaptive bandwidth selection using Lepskii's method
# for both combined and target-only estimators, against oracle choices, under manifold setting.
# The function returns MSEs and selected smoothness levels (beta) for each estimator.
cores = 10



# ------------------------------------------------------------------------
# Function: calculate_mse_fixed_nP_lepski_oracle_split
# Purpose:  Compare the MSE performance of local polynomial estimators
#           under covariate shift using (1) Lepskii-adaptive bandwidth 
#           and (2) Oracle bandwidth, for both Combined and Target-only data.
#           This function simulates data under a manifold setting.
#
# Inputs:
#   - target_sample_sizes: vector of target domain sample sizes
#   - n_s: fixed number of source domain samples
#   - x0_interior: evaluation point on the manifold
#   - alpha_true: true smoothness of the regression function
#   - num_simulations: number of repetitions per setting
#   - degree: polynomial degree for local regression
#   - kernel_type: kernel type for local polynomial regression
#
# Outputs:
#   - results_df: data.frame with MSEs, SEs, and estimator labels
#   - beta_df_combined: selected beta values for Combined estimator
#   - beta_df_target: selected beta values for Target-only estimator
# ------------------------------------------------------------------------

calculate_mse_fixed_nP_lepski_oracle_split <- function(
    target_sample_sizes, n_s, x0_interior, alpha_true, d,
    num_simulations, degree, kernel_type = "gaussian"
) {
  simulation_results <- list()
  beta_record_combined <- list()
  beta_record_target <- list()
  
  for (i in seq_along(target_sample_sizes)) {
    n_t <- target_sample_sizes[i]
    
    # Storage for simulations
    beta_values_combined <- numeric(num_simulations)
    beta_values_target <- numeric(num_simulations)
    h_combined_record <- numeric(num_simulations)
    h_target_record <- numeric(num_simulations)
    mse_c_sim <- oracle_c_sim <- oracle_t_sim <- numeric(num_simulations)
    
    sim_results <- mclapply(1:num_simulations, function(j) {
      set.seed(123 + j)
      # Generate data
      x_s <- matrix(runif(n_s * D, 0, 1), ncol = D)
      y_s <- f(x_s, x0_interior, alpha_true) + rnorm(n_s)
      
      z_t <- matrix(runif(n_t * d, -1, 1), ncol = d)
      x_t <- phi(z_t)
      y_t <- f(x_t, x0_interior, alpha_true) + rnorm(n_t)
      
      x_combined <- rbind(x_s, x_t)
      y_combined <- c(y_s, y_t)
      
      # -----------------------------
      # Lepskii Bandwidth: Combined
      # -----------------------------
      n_eff <- n_t + n_s
      beta_candidates <- seq(1, 5, by = 1 / log(n_eff))
      C_combined <- 0.1
      best_beta_c <- NULL
      
      for (beta in beta_candidates) {
        # h_beta <- (n_t + n_s^((2 * beta + d) / (2 * beta + D)))^(-1 / (2 * beta + d))
        h_beta <-  1.2*((n_t + n_s^((2 * beta + d) / (2 * beta + D)))/log(n_eff))^(-1 / (2 * beta + d))
        stable <- TRUE
        for (eta in beta_candidates) {
          if (eta <= beta) {
            ## h_eta <- (n_t + n_s^((2 * eta + d) / (2 * eta + D)))^(-1 / (2 * eta + d))
            h_eta <- 1.2*((n_t + n_s^((2 * eta + d) / (2 * eta + D)))/log(n_eff))^(-1 / (2 * eta + d))
            diff <- abs(local_polynomial_estimator(x_combined, y_combined, x0_interior, degree, h_beta, kernel_type) -
                          local_polynomial_estimator(x_combined, y_combined, x0_interior, degree, h_eta, kernel_type))
            if (diff > C_combined * h_eta^eta) {
              stable <- FALSE
              break
            }
          }
        }
        if (stable) best_beta_c <- beta
      }
      
      if (is.null(best_beta_c)) best_beta_c <- min(beta_candidates)
      beta_values_combined[j] <- best_beta_c
      h_combined <- 1.2*((n_t + n_s^((2 * best_beta_c + d) / (2 * best_beta_c + D)))/log(n_eff))^(-1 / (2 * best_beta_c + d))
      h_combined_record[j] <- h_combined
      
      
      # -----------------------------
      # Oracle Bandwidth
      # -----------------------------
      oracle_bandwidth_combined <- 1.2 * (n_t + n_s^((2 * alpha_true + d) / (2 * alpha_true + D)))^(-1 / (2 * alpha_true + d))
      oracle_bandwidth_target <- 0.8* n_t^(-1 / (2 * alpha_true + d))
      
      # Estimation
      m_hat_c <- local_polynomial_estimator(x_combined, y_combined, x0_interior, degree, h_combined, kernel_type)
      m_hat_c_oracle <- local_polynomial_estimator(x_combined, y_combined, x0_interior, degree, oracle_bandwidth_combined, kernel_type)
      m_hat_t_oracle <- local_polynomial_estimator(x_t, y_t, x0_interior, degree, oracle_bandwidth_target, kernel_type)
      
      # MSE
      truth <- f(matrix(x0_interior, ncol = d, byrow = TRUE), x0_interior, alpha_true)
      
      return(list(
        mse_c = (m_hat_c - truth)^2,
        mse_c_oracle = (m_hat_c_oracle - truth)^2,
        mse_t_oracle = (m_hat_t_oracle - truth)^2,
        beta_c = best_beta_c,
        h_c = h_combined
      ))
    }, mc.cores = cores)
    
    
    mse_c_sim         <- sapply(sim_results, function(res) res$mse_c)
    oracle_c_sim      <- sapply(sim_results, function(res) res$mse_c_oracle)
    oracle_t_sim      <- sapply(sim_results, function(res) res$mse_t_oracle)
    beta_values_combined <- sapply(sim_results, function(res) res$beta_c)
    h_combined_record <- sapply(sim_results, function(res) res$h_c)
    
    cat(sprintf("Target Size: %d | Mean β_combined: %.4f (SD: %.4f) | Mean h_comb: %.6f | Mean MSE_c: %.6f\n",
                n_t, mean(beta_values_combined), sd(beta_values_combined), mean(h_combined_record), mean(mse_c_sim)))
    cat(sprintf("              | Oracle MSE_c: %.6f | Oracle MSE_t: %.6f\n", 
                mean(oracle_c_sim), mean(oracle_t_sim)))
    beta_record_combined[[i]] <- data.frame(target_sample_size = n_t, beta_selected = beta_values_combined)
    beta_record_target[[i]] <- data.frame(target_sample_size = n_t, beta_selected = beta_values_target)
    
    simulation_results[[i]] <- data.frame(
      target_sample_size = n_t,
      mse_combined = mean(mse_c_sim),
      ## mse_target = mean(mse_t_sim),
      se_combined = sd(mse_c_sim) / sqrt(num_simulations),
      ## se_target = sd(mse_t_sim) / sqrt(num_simulations),
      oracle_mse_combined = mean(oracle_c_sim),
      oracle_mse_target = mean(oracle_t_sim),
      oracle_se_combined = sd(oracle_c_sim) / sqrt(num_simulations),
      oracle_se_target = sd(oracle_t_sim) / sqrt(num_simulations)
    )
  }
  
  # Format results for plotting
  combined_results <- bind_rows(simulation_results) %>%
    pivot_longer(
      cols = c(mse_combined, oracle_mse_combined, oracle_mse_target),
      names_to = "estimator", values_to = "mse"
    ) %>%
    mutate(
      se = case_when(
        estimator == "mse_combined" ~ se_combined,
        ## estimator == "mse_target" ~ se_target,
        estimator == "oracle_mse_combined" ~ oracle_se_combined,
        estimator == "oracle_mse_target" ~ oracle_se_target
      ),
      estimator = case_when(
        estimator == "mse_combined" ~ "Lepskii-Combined",
        ## estimator == "mse_target" ~ "Lepskii-TargetOnly",
        estimator == "oracle_mse_combined" ~ "Oracle-Combined",
        estimator == "oracle_mse_target" ~ "Oracle-TargetOnly"
      ),
      upper_bound = mse + 1.96 * se,
      lower_bound = mse - 1.96 * se
    )
  
  beta_df_combined <- do.call(rbind, beta_record_combined)
  beta_df_target <- do.call(rbind, beta_record_target)
  
  return(list(
    results_df = combined_results,
    beta_df_combined = beta_df_combined,
    beta_df_target = beta_df_target
  ))
}






## #This includes calculate_mse_fixed_nP_lepski_oracle_split

# ----------------------------- #
#      General Parameters       #
# ----------------------------- #
num_simulations <- 50
local_poly_degree <- 2
ambient_dim <- 5
true_smoothness <- 2.5
source_sample_size <- 5000
set.seed(123)

phi <- function(z) {
  cbind((z + 1) / 2,
        z[, 1]^2,
        z[, 2]^2,
        (z[, 1] + 1) * (z[, 2] + 1) / 4)
}

regression_function <- function(x, x0, alpha) {
  rowSums(abs(x - matrix(x0, nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^alpha)
}

z0 <- runif(2, -1, 1)
print(z0)
x0 <- phi(matrix(z0, nrow = 1))

# Target sample sizes to test
target_sample_sizes <- c(1000, 2000, 3000, 4000, 5000,
                         10000, 15000, 20000, 25000, 30000)


# Storage for all simulation results
all_results <- list()
all_beta_combined <- list()
all_beta_target <- list()

D = 5
# ----------------------------- #
#  Loop Over Sample Size Cases  #
# ----------------------------- #
for (n_t in target_sample_sizes) {
  cat("\nRunning for target sample size =", n_t, "\n")
  
  # Step 1: Generate target data for estimating d
  z_t <- matrix(runif(n_t * 2, -1, 1), ncol = 2)
  x_t <- phi(z_t)
  
  # Step 2: Estimate intrinsic dimension from x_t
  estimated_d <- estimate_intrinsic_dimension_manifold(x_t, k = 10)
  cat(sprintf("→ Estimated intrinsic dimension d̂ = %d\n", estimated_d))
  
  # Step 3: Run simulation using this estimated d
  res <- calculate_mse_fixed_nP_lepski_oracle_split(
    target_sample_sizes = c(n_t),
    n_s = source_sample_size,
    x0_interior = x0,
    alpha_true = true_smoothness,
    d = estimated_d,
    num_simulations = num_simulations,
    degree = local_poly_degree,
    kernel_type = "gaussian"
  )
  
  all_results[[as.character(n_t)]] <- res$results_df
  all_beta_combined[[as.character(n_t)]] <- res$beta_df_combined
  all_beta_target[[as.character(n_t)]] <- res$beta_df_target
}

# ----------------------------- #
#     Combine All Results       #
# ----------------------------- #
results_adaptive_all <- bind_rows(all_results)
beta_combined_all <- bind_rows(all_beta_combined)
beta_target_all <- bind_rows(all_beta_target)


save(results_adaptive_all,
     beta_combined_all,
     beta_target_all,
     all_results,
     all_beta_combined,
     all_beta_target,
     file = "adaptive_lepskii_oracle_results.RData")



timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
save(results_adaptive_all,
     beta_combined_all,
     beta_target_all,
     all_results,
     all_beta_combined,
     all_beta_target,
     file = paste0("lepskii_oracle_results_", timestamp, ".RData"))



# ----------------------------------------- #
#     Visualization: MSE under Covariate    #
#     Shift — Oracle vs Data-Driven Models  #
# ----------------------------------------- #

load("adaptive_lepskii_oracle_results.RData")




# Preprocess — Scale target sample size
results_adaptive_all <- results_adaptive_all %>%
  mutate(target_sample_size_scaled = target_sample_size / 1000)

target_breaks <- unique(results_adaptive_all$target_sample_size_scaled)

# Rename estimators to standardized form
results_adaptive_all <- results_adaptive_all %>%
  mutate(estimator = case_when(
    estimator == "Lepskii-Combined"     ~ "Adaptive (Pooled LPR)",
    estimator == "Oracle-Combined"      ~ "Oracle (Pooled LPR)",
    estimator == "Oracle-TargetOnly"    ~ "Oracle (Target only)",
    estimator == "Lepskii-TargetOnly"   ~ "Adaptive (Target only)",
    TRUE ~ estimator
  ))

# Filter only selected estimators for plotting
results_to_plot <- results_adaptive_all %>%
  filter(estimator %in% c("Adaptive (Pooled LPR)", "Oracle (Pooled LPR)", "Oracle (Target only)"))

mse_plot <- ggplot(results_to_plot, aes(x = target_sample_size_scaled, y = mse, color = estimator)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound), width = 0.5) +
  scale_x_continuous(breaks = target_breaks) +
  scale_y_continuous(labels = scales::scientific) +
  labs(
    title = "MSE Comparison: Adaptive vs Oracle Estimators",
    x = expression("Target Sample Size (" * 10^3 * ")"),
    y = "MSE (linear scale)",
    color = "Estimator"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.title.position = "plot",  
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12),
    ## legend.title = element_blank(),  
    legend.position = "right",                # 
    legend.justification = "center",          # 
    legend.box = "vertical",                  # 
    legend.margin = margin(0, 15, 0, 15),     # 
    legend.text = element_text(size = 11)
  )



# 
print(mse_plot)


ggsave("mse_comparison_adaptive_vs_oracle.png", 
       plot = mse_plot + theme(legend.position = "right",
                               legend.box = "vertical",
                               legend.justification = "center"),
       width = 8, height = 4.5, dpi = 300)



