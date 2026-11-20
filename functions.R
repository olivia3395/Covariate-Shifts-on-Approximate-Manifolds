# functions.R
  
# -----------------------------------------------
# Penalized Local Polynomial Estimator  (Section 2.2, Eq. (2.7))
# -----------------------------------------------
#   f_hat(x0) = e_1' (Z'WZ + lambda_h D_h)^{-1} Z'WY
#
#   Z      : i-th row is z((X_i - x0)/h), the scaled monomials
#            prod_j u_j^{a_j} / a_j!  with |a| <= degree
#   W      : diag(K((X_i - x0)/h))
#   D_h    : diag(h^{-2|a|} : |a| <= degree)
#   M_h    : sum_i W_ii
#   lambda : 1 + M_h * h^{2 * beta}
#
# `beta` is the smoothness index that enters the ridge parameter lambda_h
# (the true beta for the oracle estimator, the candidate beta for Lepski).
# Setting penalized = FALSE gives back the unregularized LPR of Eq. (2.5).
# The estimator is always well defined when penalized = TRUE (it returns 0
# if no observation falls in the kernel window).
# -----------------------------------------------

# All multi-indices a in N^D with |a| <= degree (first row is a = 0)
poly_multi_indices <- function(D, degree) {
  idx <- matrix(0L, nrow = 1, ncol = D)
  if (degree >= 1) {
    for (total in 1:degree) {
      prev <- idx[rowSums(idx) == total - 1, , drop = FALSE]
      cand <- do.call(rbind, lapply(1:D, function(j) {
        out <- prev
        out[, j] <- out[, j] + 1L
        out
      }))
      idx <- rbind(idx, unique(cand))
    }
  }
  idx
}

kernel_weights <- function(u, kernel_type = "gaussian") {
  D <- ncol(u)
  if (kernel_type == "gaussian") {
    exp(-0.5 * rowSums(u^2)) / ((2 * pi)^(D / 2))
  } else if (kernel_type == "epanechnikov") {
    pmax(0, 1 - rowSums(u^2))
  } else if (kernel_type == "box") {
    as.numeric(rowSums(abs(u) > 1) == 0)
  } else {
    stop("Unsupported kernel type")
  }
}

local_polynomial_estimator <- function(x, y, x0_targetonly, degree, bandwidth,
                                       kernel_type = "gaussian", beta = NULL,
                                       penalized = TRUE) {
  if (!is.matrix(x)) x <- as.matrix(x)
  if (penalized && is.null(beta)) {
    stop("`beta` is required for the penalized estimator: lambda_h = 1 + M_h * h^(2 * beta)")
  }
  D <- ncol(x)
  x0 <- as.numeric(x0_targetonly)

  u <- sweep(x, 2, x0, "-") / bandwidth
  w <- kernel_weights(u, kernel_type)

  # Only observations with positive weight contribute to Z'WZ and Z'WY
  keep <- which(w > 0)
  if (length(keep) == 0) {
    if (penalized) return(0)
    warning("No points within kernel support for estimation.")
    return(NA_real_)
  }
  u <- u[keep, , drop = FALSE]
  w <- w[keep]
  yk <- y[keep]

  # Design matrix Z with rows z(u_i)
  idx <- poly_multi_indices(D, degree)
  Z <- matrix(1, nrow = length(keep), ncol = nrow(idx))
  for (k in seq_len(nrow(idx))) {
    for (j in which(idx[k, ] > 0)) {
      a <- idx[k, j]
      Z[, k] <- Z[, k] * u[, j]^a / factorial(a)
    }
  }

  ZtWZ <- crossprod(Z, Z * w)
  ZtWY <- crossprod(Z, w * yk)

  if (penalized) {
    M_h <- sum(w)
    lambda_h <- 1 + M_h * bandwidth^(2 * beta)
    D_h <- bandwidth^(-2 * rowSums(idx))
    A <- ZtWZ + diag(lambda_h * D_h, nrow = nrow(idx))
    theta <- solve(A, ZtWY)
  } else {
    theta <- tryCatch(solve(ZtWZ, ZtWY), error = function(e) rep(NA_real_, nrow(idx)))
  }
  as.numeric(theta[1])
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
      
     
      m_hat_combined <- local_polynomial_estimator(x_combined, y_combined, x0_targetonly , degree, bandwidth_combined, kernel_type, beta = alpha)
      m_hat_target <- local_polynomial_estimator(x_t, y_t, x0_targetonly , degree, bandwidth_target, kernel_type, beta = alpha)
      
     
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
      m_hat_combined <- local_polynomial_estimator(x_combined, y_combined, x0_manifold, degree, bandwidth_combined, kernel_type, beta = alpha)
      m_hat_target <- local_polynomial_estimator(x_t, y_t, x0_manifold, degree, bandwidth_target, kernel_type, beta = alpha)
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
