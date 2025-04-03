# Computing the integral term ---------------------------------------
# Compute Fy(dx)
fisher.esd = function (y1, y2) 
{c2 = y2
c1 = y1
h = sqrt(c1 + c2 - c1 * c2)
a = ((1 - h)^2)/(1 - c2)^2
b = ((1 + h)^2)/(1 - c2)^2
function(x) {
  result = rep(0, length(x))
  result[x > a & x < b] = ((1 - c2) * sqrt((b - x) * (x - a)))/(2 * pi * x * (c1 + x * c2))
  return(result)
}
}

# Computes the support for the Fisher ESD
construct.fisher.support = function (y1, y2) 
{ c2 = y2
c1 = y1
h = sqrt(c1 + c2 - c1 * c2)
a = ((1 - h)^2)/(1 - c2)^2
b = ((1 + h)^2)/(1 - c2)^2
return(c(a, b))
}

# Calculates f*(x)
function.prod = function(f1,f2){
  return(function(x){ return(f1(x)*f2(x))})
}  

### Calculates the second term after minus in the test statistic T tilda
calculate.expected.trace = function (y1, y2) 
{ asymptotic.pdf = fisher.esd(y1, y2)
integrand = function.prod(function(x) {
  (1 - x)^2 + (1 - 1/x)^2
}, asymptotic.pdf)
asymptotic.supports = construct.fisher.support(y1, y2)
safe_integral = safely(integrate)
integral = exec(safe_integral, integrand, `!!!`(asymptotic.supports))[[1]]
return(integral)
}
# Calculating the T(S1,S2)--------------------------------------
# Covariance distance estimator between Sigma1 and Sigma2 for tau
covariance.distance.estimator <- function (data, epsilon, f) {
  products = purrr::map(as.data.frame(t(data)), ~.x %*% t(.x))
  forward.cumsum = purrr::accumulate(products, ~.x + .y)
  backward.cumsum = purrr::map(forward.cumsum, ~-.x + forward.cumsum[[nrow(data)]]) 
  function(tau) { 
    sigma1 = (1/tau) * forward.cumsum[[tau]] + epsilon * diag(ncol(data)) 
    sigma2 = (1/(nrow(data) - tau)) * backward.cumsum[[tau]] + epsilon * diag(ncol(data)) 
    output = tryCatch(cov.dist(sigma1, sigma2, f), error = function(error_message) { 
      return(NA)
    })
    return(output)
  }
}

# Computes the distance between two covariance matrices using eigenvalues
cov.dist= function(sigma1, sigma2, f){
  A = geigen::geigen(sigma2, sigma1, symmetric=TRUE)$values 
  return(rlang::exec(f,(A-1)^2) + rlang::exec(f,(1/A-1)^2))
}


# Calculating mu(gamma) and sigma_square ----------------------------------
# Calculates mu(gamma)
asymptotic.bias = function (y1, y2) 
{ h = sqrt(y1 + y2 - y1 * y2)
K_1 = 2 * h * (1 + h^2)/(1 - y2)^4 - 2 * h/(1 - y2)^2
J_1 = 2 * h * (1 + h^2)/(1 - y1)^4 - 2 * h/(1 - y1)^2
return(2 * (h^2 - y2^2)/(1 - y2)^4 + 2 * K_1 * y2/h + 2 * 
         (h^2 - y1^2)/(1 - y1)^4 + 2 * J_1 * y1/h)
}

# calculate the sigma square function, Wrong formula in the return
asymptotic.variance = function (y1, y2) 
{ h = sqrt(y1 + y2 - y1 * y2)
K_21 = 2 * h * (1 + h^2)/(1 - y2)^4 - 2 * h/(1 - y2)^2
K_22 = 2 * h * (1 + h^2)/(1 - y1)^4 - 2 * h/(1 - y1)^2
K_31 = h^2/(1 - y2)^4
K_32 = h^2/(1 - y1)^4
J_1 = -2 * (1 - y2)^2
J_2 = (1 - y2)^4
var_x = K_21^2 + 2 * K_31^2
var_y = K_22^2 + 2 * K_32^2
cov_xy = J_1 * K_21/h + J_1 * K_21/(h * (h^2 - 1)) + (-J_1 * 
                                                        K_31 * (h^2 + 1)/h^2) + (-J_1 * K_31/(h^2 * (h^2 - 1))) 
J_2 * K_21 * h/(h^2 - 1)^3 + J_2 * K_31/h^2 + J_2 * + 
  K_31 * ((1 - 3 * h^2)/(h^2 * (h^2 - 1)^3))
return(2 * (var_x + var_y + 1 * cov_xy)) 
return(limiting.var)
}

# Test Statistic ----------------------------------------------------------
# Calculates the statistic for change point detection
matrix.dist.test.stat = function (data, minseglen) 
{ p = ncol(data)
n = nrow(data)
t = seq(p + minseglen + 1, n - p - minseglen) 
estimate.covariance.distance = covariance.distance.estimator(data, 0, mean) 
test.stat = purrr::map_dbl(t, estimate.covariance.distance)
dimension.over.length = purrr::map(t, ~c(p/.x, p/(n - .x)))
trace = map(dimension.over.length, ~exec(calculate.expected.trace, !!!.x))
values = map_lgl(trace, ~length(.x[[1]]) == 0)
bias = map_dbl(dimension.over.length, ~exec(asymptotic.bias, !!!.x)) 
variance = map_dbl(dimension.over.length, ~exec(asymptotic.variance, !!!.x)) 
trace = map_if(trace, values, ~NA) %>% map_dbl(~.x[[1]]) 
bias = map_if(bias, values, ~NA)  
variance = map_if(variance, values, ~NA)  
test.stat = pmap_dbl(list(test.stat, trace, bias, variance), 
                     ~(p * (..1 - ..2) - ..3)/sqrt(..4)) 
return(c(rep(NA, p + minseglen), (test.stat), rep(NA, p + minseglen)))
}
# Results -----------------------------------------------------------------
# for different possible change points
limiting.variance = function(n,p){
  t = seq(p+1,n-p)
  dimension.over.length = purrr::map(t, ~c(p/.x, p/(n-.x)))
  variance = map_dbl(dimension.over.length, ~exec(asymptotic.variance,!!!.x) ) 
  return(variance)
}

# calculates from where to start the sequenc eof possible change points
calculate.minseglen = function(n, p, func, alpha=2, constraint=1/n){
  if(func == "linear"){
    return(alpha*p)
  }  else if( func == "log-linear"){
    return(log(n)*p)
  } else if( func == "constrained-gradient"){
    grad = limiting.variance(n,p) %>% diff %>% abs 
    return(p + min(which(grad < constraint)))
  }
}

# Obtain a critical value
bonferoni = function(n,alpha){ qnorm(1-.05/n)}
library(tidyverse)

simulations_ryan_killick <- function(n_simulations, delta, n, p, cp, minseglen = NULL) {
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
  
  # Calculate candidate points
  t_0 <- seq(p + minseglen + 1, n - p - minseglen)
  
  # Critical value (using Bonferroni correction)
  critical_value <- bonferoni(length(t_0))
  
  # Main simulation loop
  for (i in 1:n_simulations) {
    # Data generation with change point
    Y <- matrix(rnorm(n*p, 0, 1), nrow = n)
    if(cp < n && cp > 0) {  # Ensure cp is within valid range
      Y[cp:n, ] <- delta * Y[cp:n, ]
    }
    
    # Calculate test statistics
    test_stats <- matrix.dist.test.stat(Y, minseglen)
    test_stats <- test_stats[!is.na(test_stats)]  # Remove NAs
    
    if(length(test_stats) == 0) {
      results[i, ] <- c(NA, NA, 0, NA)
      next
    }
    
    # Find maximum statistic
    max_idx <- which.max(abs(test_stats))
    max_stat <- test_stats[max_idx]
    detected_cp <- t_0[max_idx]
    
    # Determine significance
    is_sig <- as.integer(abs(max_stat) > critical_value)
    
    # Store results
    results[i, ] <- c(
      max_stat,
      detected_cp,
      is_sig
    )
  }
  
  return(results)
}

# # Example usage
# set.seed(1299)
# sim_results <- simulations_ryan_killick(
#   n_simulations = 10,
#   n = 500,
#   p = 10,
#   delta = 1,
#   cp = 250,  # floor(500/2)
#   minseglen = max(4*p, 30)
# )
# 
# # View results
# head(sim_results)
