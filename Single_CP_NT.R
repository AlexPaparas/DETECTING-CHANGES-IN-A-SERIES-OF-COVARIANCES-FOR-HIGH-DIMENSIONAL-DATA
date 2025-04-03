# Core Nagao Method Functions ----------------------------------------------
compute_nagao_components <- function(data_before, data_after, p, n, t) {
  # Regularized covariance matrices
  Sigma_1 <- cov(data_before) + 1e-8 * diag(p)
  Sigma_2 <- cov(data_after) + 1e-8 * diag(p)
  
  # Eigen decomposition
  F_t <- Sigma_1 %*% solve(Sigma_2)
  values <- eigen(F_t, symmetric = FALSE, only.values = TRUE)$values
  
  # Dimension ratios
  y1 <- p/(t-1)
  y2 <- p/(n-t-1)
  eta <- (t-1)/(n-2)
  
  # Asymptotic parameters
  mu <- y1 + y2 + 2*y1*y2
  sigma_sq <- 8*((y1 + y2)^2 + 
                   2*(y1 + y2)*(y1^2 + y2^2 - y1*y2) + 
                   y1*y2*(2*y1 - y2)*(2*y2 - y1))
  u <- y1 + y2
  
  list(values = values, y1 = y1, y2 = y2, eta = eta,
       mu = mu, sigma_sq = sigma_sq, u = u)
}

compute_nagao_statistic <- function(components, p) {
  with(components, {
    statistic <- sum(((values - 1)^2)/((eta*values + (1 - eta))^2))
    (statistic - p*u - mu)/sqrt(sigma_sq)
  })
}

# Main Analysis Function --------------------------------------------------
analyze_nagao <- function(Y, minseglen = NULL, alpha = 0.05) {
  # Validate input
  if (!is.matrix(Y) && !is.data.frame(Y)) {
    stop("Input data must be a matrix or dataframe")
  }
  
  # Set dimensions and parameters
  p <- ncol(Y)
  n <- nrow(Y)
  if (is.null(minseglen)) minseglen <- max(4*p, 30)
  
  # Candidate change points
  t_0 <- seq(p + minseglen + 1, n - p - minseglen)
  
  # Pre-allocate results
  results <- matrix(nrow = length(t_0), ncol = 3)
  colnames(results) <- c("candidate_cp", "standardized_stat", "raw_stat")
  
  # Compute statistics for each candidate point
  for (i in seq_along(t_0)) {
    t <- t_0[i]
    components <- compute_nagao_components(
      data_before = Y[1:t, , drop = FALSE],
      data_after = Y[(t+1):n, , drop = FALSE],
      p = p, n = n, t = t
    )
    
    results[i, ] <- c(
      t,
      compute_nagao_statistic(components, p),
      sum(((components$values - 1)^2)/((components$eta*components$values + 
                                          (1-components$eta))^2))
    )
  }
  
  # Find maximum statistic
  max_idx <- which.max(results[, "standardized_stat"])
  
  # Return comprehensive results
  list(
    test_statistics = results,
    change_point = results[max_idx, "candidate_cp"],
    max_statistic = results[max_idx, "standardized_stat"],
    critical_value = qnorm(1 - alpha/n),
    is_significant = results[max_idx, "standardized_stat"] > qnorm(1 - alpha/n),
    minseglen = minseglen,
    dimensions = c(n = n, p = p)
  )
}

# # Example Usage -----------------------------------------------------------
# 
# # Using your data
# data <- read.csv("Monthly_Prices.csv")[,4:22]
# 
# # Run analysis
# nagao_results <- analyze_nagao(data, minseglen = 20)
# 
# # View results
# print(nagao_results[c("change_point", "max_statistic", "is_significant")])