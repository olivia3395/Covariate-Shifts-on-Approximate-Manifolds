# functions.R
  
# -----------------------------------------------
# Local Polynomial Estimator
# -----------------------------------------------
local_polynomial_estimator <- function(x, y, x0_targetonly, degree, bandwidth, kernel_type = "gaussian") {
  if (!is.matrix(x)) x <- as.matrix(x)
  d <- ncol(x)
  n <- nrow(x)
  
  kernel_function <- function(u, kernel_type) {
    if (kernel_type == "gaussian") {
      return(exp(-0.5 * rowSums(u^2)) / ((2 * pi)^(d / 2)))
    } else if (kernel_type == "epanechnikov") {
      return(pmax(0, 1 - rowSums(u^2)))
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
    paste("y ~", paste(term_names, collapse = " + "), "- 1")
  } else {
    paste("y ~", paste(term_names, collapse = " + "))
  }
  
  lm_fit <- lm(as.formula(formula_str), data = data_frame, weights = weights)
  
  new_data <- as.data.frame(matrix(0, nrow = 1, ncol = length(term_names)))
  colnames(new_data) <- term_names
  if ("Intercept" %in% term_names) new_data$Intercept <- 1
  
  return(predict(lm_fit, newdata = new_data)[1])
}

# -----------------------------------------------
# Target Function f(x)
# -----------------------------------------------
f <- function(x, x0, alpha) {
  rowSums(abs(x - matrix(x0, nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^alpha)
}

# -----------------------------------------------
# MSE Calculation under Fixed Source Size
# -----------------------------------------------

calculate_mse_fixed_nP <- function(target_sample_sizes, n_s, x0_targetonly , alpha, num_simulations, degree, kernel_type = "gaussian") {
  
  simulation_results <- list()
  
  for (i in seq_along(target_sample_sizes)) {
    n_t <- target_sample_sizes[i]  
    
    
    n_phi_combined <- n_s * ifelse(all(x0_targetonly  >= -0.5 & x0_targetonly  <= 0.5), 1 / (1^d), 0) + 
      n_t * ifelse(all(x0_targetonly  >= 0 & x0_targetonly  <= 1), 1 / (1^d), 0)
    
    n_phi_target <- n_t * ifelse(all(x0_targetonly  >= 0 & x0_targetonly  <= 1), 1 / (1^d), 0)
    
   
    bandwidth_combined <- (n_phi_combined)^(-1 / (2 * alpha + d))
    bandwidth_target <- (n_phi_target)^(-1 / (2 * alpha + d))
    
   
    mse_c_simulations <- numeric(num_simulations)
    mse_t_simulations <- numeric(num_simulations)
    
    for (j in 1:num_simulations) {
     
      x_s <- matrix(runif(n_s * d, min = -1/2, max = 1/2), ncol = d)
      y_s <- f(x_s, x0_targetonly , alpha) + rnorm(n_s)  
      
      x_t <- matrix(runif(n_t * d, min = 0, max = 1), ncol = d)
      y_t <- f(x_t, x0_targetonly , alpha) + rnorm(n_t)
      
      
      x_combined <- rbind(x_s, x_t)
      y_combined <- c(y_s, y_t)
      
     
      m_hat_combined <- local_polynomial_estimator(x_combined, y_combined, x0_targetonly , degree, bandwidth_combined, kernel_type)
      m_hat_target <- local_polynomial_estimator(x_t, y_t, x0_targetonly , degree, bandwidth_target, kernel_type)
      
     
      mse_c_simulations[j] <- (m_hat_combined - f(matrix(x0_targetonly , ncol = d, byrow = TRUE), x0_targetonly , alpha))^2
      mse_t_simulations[j] <- (m_hat_target - f(matrix(x0_targetonly , ncol = d, byrow = TRUE), x0_targetonly , alpha))^2
    }
    
  
    mse_c_simulations[is.na(mse_c_simulations)] <- 0
    mse_t_simulations[is.na(mse_t_simulations)] <- 0
    
   
    simulation_results[[i]] <- data.frame(
      target_sample_size = n_t,
      mse_combined = mean(mse_c_simulations, na.rm = TRUE),
      mse_target = mean(mse_t_simulations, na.rm = TRUE),
      se_combined = sd(mse_c_simulations, na.rm = TRUE) / sqrt(num_simulations),
      se_target = sd(mse_t_simulations, na.rm = TRUE) / sqrt(num_simulations)
    )
  }
  
  
  combined_results <- bind_rows(simulation_results) %>%
  
    pivot_longer(cols = c(mse_combined, mse_target), 
                 names_to = "estimator", 
                 values_to = "mse") %>%
    pivot_longer(cols = c(se_combined, se_target), 
                 names_to = "se_estimator", 
                 values_to = "se") %>%
  
    mutate(
      estimator = case_when(
        estimator == "mse_combined" ~ "Combined",
        estimator == "mse_target" ~ "Target-Only"
      ),
      upper_bound = mse + 1.96 * se,
      lower_bound = mse - 1.96 * se
    )
  
  return(combined_results)
}


# Smooth nonlinear mapping from manifold to ambient space
phi <- function(z) {
  cbind(
    (z + 1) / 2,
    z[, 1]^2,
    z[, 2]^2,
    (z[, 1] + 1) * (z[, 2] + 1) / 4
  )
}




calculate_mse_fixed_nP_manifold <- function(target_sample_sizes, n_s, x0_manifold, alpha, num_simulations, degree, kernel_type = "gaussian") {
  simulation_results <- list()
  
  for (i in seq_along(target_sample_sizes)) {
    n_t <- target_sample_sizes[i]
    
    # Bandwidth definitions
    bandwidth_combined <- (n_t + n_s^((2 * alpha + d) / (2 * alpha + D)))^(-1 / (2 * alpha + d))
    bandwidth_target <- n_t^(-1 / (2 * alpha + d))
    
    mse_c_simulations <- numeric(num_simulations)
    mse_t_simulations <- numeric(num_simulations)
    
    for (j in 1:num_simulations) {
      # Generate source data in high dimensions (Uniform in [0,1]^D)
      ## x_s <- matrix(runif(n_s * D, 0, 1), ncol = D)
      x_s <- matrix(runif(n_s * D, 0, 1), ncol = D)
      
      y_s <- f(x_s, x0_manifold, alpha) + rnorm(n_s)
      
      # Generate target data on the manifold
      ## z_t <- matrix(runif(n_t * d, -1, 1), ncol = d)
      z_t <- matrix(runif(n_t * d, -1, 1), ncol = d)
      rho_n <- 0
      u_t <- matrix(runif(n_t * D, -1, 1), ncol = D)
      x_t <- phi(z_t) + rho_n * u_t
      
      ## x_t <- phi(z_t)
      y_t <- f(x_t, x0_manifold, alpha) + rnorm(n_t)
      
      # Combine data
      x_combined <- rbind(x_s, x_t)
      y_combined <- c(y_s, y_t)
      
      # Estimation
      m_hat_combined <- local_polynomial_estimator(x_combined, y_combined, x0_manifold, degree, bandwidth_combined, kernel_type)
      m_hat_target <- local_polynomial_estimator(x_t, y_t, x0_manifold, degree, bandwidth_target, kernel_type)
      true_val <- f(matrix(x0_manifold, ncol = D, byrow = TRUE), x0_manifold, alpha)
      
      mse_c_simulations[j] <- (m_hat_combined - true_val)^2
      mse_t_simulations[j] <- (m_hat_target - true_val)^2
    }
    
    # Handle NA
    mse_c_simulations[is.na(mse_c_simulations)] <- 0
    mse_t_simulations[is.na(mse_t_simulations)] <- 0
    
    simulation_results[[i]] <- data.frame(
      target_sample_size = n_t,
      mse_combined = mean(mse_c_simulations, na.rm = TRUE),
      mse_target = mean(mse_t_simulations, na.rm = TRUE),
      se_combined = sd(mse_c_simulations, na.rm = TRUE) / sqrt(num_simulations),
      se_target = sd(mse_t_simulations, na.rm = TRUE) / sqrt(num_simulations)
    )
  }
  
  combined_results <- bind_rows(simulation_results) %>%
    pivot_longer(cols = c(mse_combined, mse_target), names_to = "estimator", values_to = "mse") %>%
    pivot_longer(cols = c(se_combined, se_target), names_to = "se_estimator", values_to = "se") %>%
    mutate(
      estimator = case_when(
        estimator == "mse_combined" ~ "Combined",
        estimator == "mse_target" ~ "Target-Only"
      ),
      upper_bound = mse + 1.96 * se,
      lower_bound = mse - 1.96 * se
    )
  
  return(combined_results)
}
