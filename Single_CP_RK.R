source("Simulations_Ryan_Killick.R")

analyze_rk <- function(data, minseglen = NULL, alpha = 0.05) {
  # Validate input
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input data must be a dataframe or matrix")
  }
  
  # Remove non-numeric columns if needed
  data <- data[, sapply(data, is.numeric)]
  
  # Set default minseglen if not provided
  p <- ncol(data)
  n <- nrow(data)
  if (is.null(minseglen)) {
    minseglen <- max(4 * p, 30)
  }
  
  # Calculate test statistics
  result <- matrix.dist.test.stat(data, minseglen)
  
  # Find change point
  tau <- which.max(abs(result))
  max_stat <- result[tau]
  
  # Calculate p-value and significance
  critical_value <- bonferoni(n, alpha)
  is_significant <- max_stat > critical_value
  
  # Return comprehensive results
  list(
    test_statistics = result,
    change_point = tau,
    max_statistic = max_stat,
    critical_value = critical_value,
    is_significant = is_significant,
    time_index = if (is_significant) tau else NA,
    minseglen = minseglen,
    dimensions = c(n = n, p = p)
  )
}

# # Usage example with your data:
# data <- read.csv("Monthly_Prices.csv")[,4:22]
# 
# # Single dataset analysis
# rk_results <- analyze_rk(data, minseglen = 20)
# print(rk_results)