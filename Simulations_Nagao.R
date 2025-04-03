simulations_nagao <- function(n_simulations, n, p, delta, cp, minseglen = NULL) {
  # Input validation
  if(length(cp) != 1 || !is.numeric(cp)) {
    stop("cp must be a single numeric value")
  }
  
  # Set default minseglen if not provided
  if(is.null(minseglen)) {
    minseglen <- max(4 * p, 30)
  }
  
  # Initialize results matrix
  results <- matrix(nrow = n_simulations, ncol = 3)
  colnames(results) <- c("max_statistic", "detected_cp", "is_significant")
  
  # Pre-compute candidate points
  t_0 <- seq(minseglen + 1, n - minseglen)
  
  # Critical value for significance (using Bonferroni correction)
  critical_value <- qnorm(1 - 0.05/length(t_0))
  
  # Main simulation loop
  for (sim in 1:n_simulations) {
    # Generate data with change point at cp
    Y <- matrix(rnorm(n * p), nrow = n)
    if(cp < n && cp > 0) {  # Ensure cp is within valid range
      Y[(cp+1):n, ] <- Y[(cp+1):n, ] * delta
    }
    
    # Compute test statistics for all candidate points
    test_stats <- sapply(t_0, function(t) {
      tryCatch({
        # Calculate components
        data_before <- Y[1:t, , drop = FALSE]
        data_after <- Y[(t+1):n, , drop = FALSE]
        
        Sigma_1 <- cov(data_before) + 1e-8 * diag(p)
        Sigma_2 <- cov(data_after) + 1e-8 * diag(p)
        
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
        
        # Calculate statistic
        statistic <- sum(((values - 1)^2)/((eta*values + (1 - eta))^2))
        (statistic - p*u - mu)/sqrt(sigma_sq)
      }, error = function(e) NA)
    })
    
    # Store results
    if(all(is.na(test_stats))) {
      results[sim, ] <- c(NA, NA, 0)
    } else {
      max_idx <- which.max(test_stats)
      max_stat <- max(test_stats, na.rm = TRUE)
      detected_cp <- t_0[max_idx]
      
      results[sim, ] <- c(
        max_stat,
        detected_cp,
        as.integer(max_stat > critical_value)
      )
    }
  }
  
  return(results)
}

# # Example usage
# set.seed(1299)
# n <- 500
# p <- 10
# delta <- 1
# cp <- floor(n / 2)
# minseglen <- max(4 * p, 30)
# 
# # Run simulations
# sim_results <- simulations_nagao(
#   n_simulations = 10,
#   n = n,
#   p = p,
#   delta = delta,
#   cp = cp,
#   minseglen = minseglen
# )
# 
# # View results
# head(sim_results)