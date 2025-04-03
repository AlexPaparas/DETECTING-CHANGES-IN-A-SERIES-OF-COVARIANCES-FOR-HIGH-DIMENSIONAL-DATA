# Function to compute covariance matrices and eigenvalues
compute_eigenvalues <- function(Y, t, p) {
  data_before <- Y[1:t, , drop = FALSE]
  data_after <- Y[(t + 1):nrow(Y), , drop = FALSE]
  
  Sigma_1 <- cov(data_before) + 1e-8 * diag(p)
  Sigma_2 <- cov(data_after) + 1e-8 * diag(p)
  
  F_t <- Sigma_1 %*% solve(Sigma_2)
  eigen(F_t, symmetric = FALSE, only.values = TRUE)$values
}

# Function to compute the u term
compute_u_term <- function(y1, y2) {
  ((y1 + y2 - y1*y2)/(y1*y2))*log((y1 + y2)/(y1 + y2 - y1*y2)) +
    (y1*(1-y2)/(y2*(y1+y2)))*log(1-y2) +
    (y2*(1-y1)/(y1*(y1+y2)))*log(1-y1)
}

# Function to compute asymptotic parameters
compute_asymptotic_params <- function(y1, y2) {
  h <- sqrt(y1 + y2 - y1*y2)
  
  a1 <- (1 + h^2)/((1-y2)^2)
  b1 <- 2*h/((1-y2)^2)
  a2 <- (y1*(1-y2)^2 + y2*(1+h^2))/((y1+y2)*(1-y2)^2)
  b2 <- 2*y2*h/((y1+y2)*(1-y2)^2)
  
  list(
    a1 = a1, b1 = b1, a2 = a2, b2 = b2,
    h = h, y1 = y1, y2 = y2
  )
}

# Function to compute mu and sigma terms
compute_mu_sigma <- function(params) {
  with(params, {
    # Square root terms
    a2_b2_sq <- sqrt(a2^2 - b2^2)
    a1_b1_sq <- sqrt(a1^2 - b1^2)
    a2_plus_b2 <- sqrt(a2 + b2)
    a2_minus_b2 <- sqrt(a2 - b2)
    a1_plus_b1 <- sqrt(a1 + b1)
    a1_minus_b1 <- sqrt(a1 - b1)
    
    # Mu terms
    mu1 <- 0.5*log(a2_b2_sq*4*h^2/(h*(a2_plus_b2 + a2_minus_b2) - y2*(a2_plus_b2 - a2_minus_b2))^2)
    mu2 <- 0.5*log(a1_b1_sq*4*h^2/(h*(a1_plus_b1 + a1_minus_b1) - y2*(a1_plus_b1 - a1_minus_b1))^2)
    mu_n_t <- mu1 - (y2/(y1+y2))*mu2
    
    # Variance terms
    cov_term <- 2*log((a1_plus_b1 + a1_minus_b1)*(a2_plus_b2 + a2_minus_b2)/
                        ((a1_plus_b1 + a1_minus_b1)*(a2_plus_b2 + a2_minus_b2) - 
                           (a1_plus_b1 - a1_minus_b1)*(a2_plus_b2 - a2_minus_b2)))
    
    var1 <- 2*log((a2 + a2_b2_sq)/(2*a2_b2_sq))
    var2 <- 2*log((a1 + a1_b1_sq)/(2*a1_b1_sq))
    
    list(
      mu_n_t = mu_n_t,
      sigma_sq = var1 + (y2/(y1+y2))^2*var2 - 2*(y2/(y1+y2))*cov_term
    )
  })
}

compute_test_statistic <- function(Y, t, p, n) {
  # Compute eigenvalues
  ev <- compute_eigenvalues(Y, t, p)
  
  # Dimension ratios
  y1 <- p/(t-1)
  y2 <- p/(n-t-1)
  c1 <- (t-1)/(n-2)
  c2 <- (n-t-1)/(n-2)
  
  # Compute statistic components
  Statistic <- sum(((t-1)/2 * log(ev) - ((n-2)/2 * log(c1*ev + c2))))
  u_term <- compute_u_term(y1, y2)
  
  # Asymptotic parameters
  params <- compute_asymptotic_params(y1, y2)
  mu_sigma <- compute_mu_sigma(params)
  
  # Standardized statistic
  ((-2/(n-2))*Statistic - p*u_term - mu_sigma$mu_n_t)/sqrt(abs(mu_sigma$sigma_sq))
}

simulations_likelihood_ratio <- function(n_simulations, n, p, delta, cp, minseglen) {
  # Initialize results
  results <- matrix(nrow = n_simulations, ncol = 3)
  colnames(results) <- c("max_statistic", "cp", "detection")
  
  # Pre-compute candidate points
  minseglen <- max(4 * p, 30)
  t_0 <- seq(p + minseglen + 1, n - p - minseglen)
  
  # Main simulation loop
  for (sim in 1:n_simulations) {
    # Generate data
    Y <- matrix(rnorm(n * p), nrow = n)
    Y[cp:n, ] <- delta * Y[cp:n, ]
    
    # Compute test statistics
    test_stats <- sapply(t_0, function(t) compute_test_statistic(Y, t, p, n))
    
    # Store results
    max_idx <- which.max(test_stats)
    results[sim, ] <- c(
      max(test_stats),
      t_0[max_idx],
      as.integer(max(test_stats) > qnorm(1 - 0.05/n))
    )
  }
  
  return(results)
}

# set.seed(1299)
# results <- simulations_likelihood_ratio(
#   n_simulations = 10,
#   n = 500,
#   p = 10,
#   delta = 1,
#   cp = 250,
#   minseglen = max(4*10, 30)
# )
# 
# head(results)