# =============================================
# load packages
# =============================================
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(ggplot2)

# =============================================
# parameters
# =============================================
D <- 5
d <- 2
alpha <- 2.5
degree <- 2
num_simulations <- 100

source_sample_sizes <- c(1e4, 2e4, 3e4, 4e4, 5e4, 6e4, 7e4, 8e4, 9e4, 1e5, 2e5, 3e5, 4e5, 5e5, 6e5, 7e5,8e5, 9e5, 1e6)

set.seed(123)


# Smooth nonlinear mapping from manifold to ambient space
phi <- function(z) {
  cbind(
    (z + 1) / 2,
    z[, 1]^2,
    z[, 2]^2,
    (z[, 1] + 1) * (z[, 2] + 1) / 4
  )
}


# Target point
z0 <- runif(d, -1, 1)
u_0<- runif(D, -1, 1)
rho_n <- 0.5
x0_manifold <- phi(matrix(z0, nrow = 1)) + matrix(rho_n * u_0, nrow = 1)



f <- function(x, x0, alpha) {
  rowSums(abs(x - matrix(x0, nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^alpha)
}


local_polynomial_estimator <- function(x, y, x0_targetonly, degree, bandwidth, kernel_type = "box") {
  if (!is.matrix(x)) x <- as.matrix(x)
  d <- ncol(x)
  n <- nrow(x)
  
  
  kernel_function <- function(u, kernel_type) {
    if (kernel_type == "gaussian") {
      return(exp(-0.5 * rowSums(u^2)) / ((2 * pi)^(d / 2)))
    } else if (kernel_type == "epanechnikov") {
      return(pmax(0, 1 - rowSums(u^2)))
    } else if (kernel_type == "box") {
      return(ifelse(rowSums(abs(u) > 1) == 0, 1, 0)) 
    } else {
      stop("Unsupported kernel type")
    }
  }
  

  u <- sweep(x, 2, x0_targetonly, "-") / bandwidth
  
  weights <- kernel_function(u, kernel_type)
  

  terms_list <- lapply(1:d, function(i) 0:degree)
  all_combinations <- expand.grid(terms_list)
  valid_combinations <- all_combinations[rowSums(all_combinations) <= degree, , drop = FALSE]
  
  data_frame <- data.frame(y = y)
  term_names <- character(nrow(valid_combinations))
  
  for (i in 1:nrow(valid_combinations)) {
    k <- valid_combinations[i, ]
    term <- rep(1, n)
    term_name_parts <- character(0)
    
    for (j in 1:d) {
      kj <- k[[j]]
      term <- term * (x[, j] - x0_targetonly[j])^kj
      if (kj > 0) {
        term_name_parts <- c(term_name_parts, paste0("X", j, "_deg", kj))
      }
    }
    
    term_name <- if (length(term_name_parts) == 0) "Intercept" else paste(term_name_parts, collapse = "_")
    term_names[i] <- make.names(term_name, unique = TRUE)
    data_frame[[term_names[i]]] <- term
  }
  
  formula_str <- if ("Intercept" %in% term_names) {
    paste("y ~", paste(term_names, collapse = " + "), "-1")
  } else {
    paste("y ~", paste(term_names, collapse = " + "))
  }
  
  if (all(weights == 0)) {
    warning("No points within kernel support for estimation.")
    return(NA_real_)
  }
  
  lm_fit <- lm(as.formula(formula_str), data = data_frame, weights = weights)
  
  new_data <- as.data.frame(matrix(0, nrow = 1, ncol = length(term_names)))
  colnames(new_data) <- term_names
  if ("Intercept" %in% term_names) new_data$Intercept <- 1
  
  return(predict(lm_fit, newdata = new_data)[1])
}


# =============================================
## parallel
# =============================================
n_cores <- 18 
cl <- makeCluster(n_cores)
registerDoParallel(cl)


library(doRNG)
registerDoRNG(123)
# =============================================

simulation_results <- list()

for (n_s in source_sample_sizes) {
  n_t <- n_s * 0.2  
  h_n_star_combined <- (n_s^((2 * alpha + d) / (2 * alpha + D)) + n_t)^(-1 / (2 * alpha + d))
  cat(sprintf("[INFO] h_n_star (combined estimator): %.5f\n", h_n_star_combined))
  
  if (rho_n >= h_n_star_combined) {
    bandwidth_combined <- 1.7 * (n_s + n_t * rho_n^(d - D))^(-1 / (2 * alpha + D))
  } else {
    bandwidth_combined <- 1.7 * (n_s^((2 * alpha + d) / (2 * alpha + D)) + n_t)^(-1 / (2 * alpha + d))
  }
  
  sim_data <- foreach(j = 1:num_simulations, .combine = rbind, .packages = "stats") %dopar% {
    # Source covariates
    x_s <- matrix(runif(n_s * D, 0, 1), ncol = D)
    y_s <- f(x_s, x0_manifold, alpha) + rnorm(n_s)
    
    # Target covariates
    z_t <- matrix(runif(n_t * d, -1, 1), ncol = d)
    u_t <- matrix(runif(n_t * D, -1, 1), ncol = D)
    x_t <- phi(z_t) + rho_n * u_t
    y_t <- f(x_t, x0_manifold, alpha) + rnorm(n_t)
    
    x_combined <- rbind(x_s, x_t)
    y_combined <- c(y_s, y_t)
    
    m_hat_c <- local_polynomial_estimator(x_combined, y_combined, x0_manifold, degree, bandwidth_combined)
    
    true_val <- f(matrix(x0_manifold, ncol = D, byrow = TRUE), x0_manifold, alpha)
    
    data.frame(
      mse_combined = (m_hat_c - true_val)^2
    )
  }
  
  simulation_results[[as.character(n_s)]] <- data.frame(
    source_sample_size = n_s,
    target_sample_size = n_t,
    rho_n = rho_n,
    h_n_star_combined = h_n_star_combined,
    mse_combined = mean(sim_data$mse_combined),
    se_combined = sd(sim_data$mse_combined) / sqrt(num_simulations)
  )
}



stopCluster(cl)

results_loglog <- bind_rows(simulation_results) %>%
  pivot_longer(cols = c(mse_combined), names_to = "estimator", values_to = "mse") %>%
  pivot_longer(cols = c(se_combined), names_to = "se_estimator", values_to = "se") %>%
  mutate(
    estimator = case_when(
      estimator == "mse_combined" ~ "Proposed"
    ),
    eff_sample_size = case_when(
      rho_n >= h_n_star_combined ~ log10(source_sample_size + target_sample_size * rho_n^(d - D)),
      rho_n < h_n_star_combined ~ log10(source_sample_size^((2 * alpha + d) / (2 * alpha + D)) + target_sample_size)
    ),
    log_mse = log10(mse),
    upper_bound = log10(mse + 1.96 * se),
    lower_bound = log10(pmax(mse - 1.96 * se, 1e-10))
  )




saveRDS(results_loglog, "results_loglog_varyrho_proposedonly.rds")


results_loglog <- readRDS("results_loglog_varyrho_proposedonly.rds")



theory_slope <- -2 * alpha / (2 * alpha + D)

theory_intercept <- 0.08

slope_label <- sprintf("Theoretical slope (d = %d): %.3f", d, theory_slope)
slope_color <- "#e55b00" 

fit_proposed <- lm(log_mse ~ eff_sample_size, data = subset(results_loglog, estimator == "Proposed"))
slope_proposed <- coef(fit_proposed)[2]



x_vals <- seq(min(results_loglog$eff_sample_size), max(results_loglog$eff_sample_size), length.out = 100)
y_vals <- theory_intercept + theory_slope * x_vals
theory_line_df <- data.frame(x = x_vals, y = y_vals, type = "Theoretical slope")



mse_loglog_plot <- ggplot(results_loglog, aes(x = eff_sample_size, y = log_mse, color = estimator)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, size = 1.5) +
  
  # Theoretical slope line with legend
  geom_line(data = theory_line_df,
            aes(x = x, y = y, linetype = type),
            color = slope_color,
            size = 1.3) +
  
  annotate("text",
           x = max(results_loglog$eff_sample_size) - 1,
           y = quantile(results_loglog$log_mse, 0.25) + 0.4,
           label = sprintf("Proposed fit slope ≈ %.3f", slope_proposed),
           color = "#5f97d2",
           hjust = 1,
           size = 6) +
  annotate("text",
           x = max(results_loglog$eff_sample_size) - 1,
           y = quantile(results_loglog$log_mse, 0.25) + 0.25,
           label = sprintf("Theoretical slope ≈ %.3f", theory_slope),
           color = slope_color,
           hjust = 1,
           size = 6) +
  
  scale_color_manual(values = c("Proposed" = "#5f97d2")) +
  scale_linetype_manual(
    name = "Reference Lines",
    values = c("Theoretical slope" = "dashed")
  ) +
  
  labs(
    title = bquote("Log-Log MSE vs Effective Sample Size  (" ~ rho[n] == 0.05 ~ ")"),
    x = expression(log[10]("Effective Sample Size")),
    y = expression(log[10](MSE)),
    color = "Estimator"
  ) +
  ## coord_cartesian(ylim = c(-4.3, -2), clip = "off") +
  coord_cartesian(ylim = c(-3, -1.2), clip = "off") + 
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.position = "top"
  )

ggsave("loglog_mse_plot_varyrho_adaptive_slope.png", plot = mse_loglog_plot,
       width = 8, height = 6, dpi = 300, bg = "white")


# #