# Source the R scripts that contain the functions
source("Simulations_Likelihood_Ratio.R")
source("Simulations_Nagao.R")
source("Simulations_Ryan_Killick.R")

# Set your parameters
n <- 500 ; p <- 10; delta <- 1; cp <- floor(n / 2); minseglen <- max(4 * p, 30); n_simulations <- 10  

# Common helper function for change point detection
check_cp_detection <- function(statistic, n, alpha = 0.05) {
  critical_value <- qnorm(1 - alpha/n)
  ifelse(statistic > critical_value, "C.P. detected", "No C.P. detected")
}

# Run the simulations for Likelihood Ratio
set.seed(1299)  # For reproducibility
sim_results_Likelihood_Ratio <- simulations_likelihood_ratio(n_simulations = n_simulations,  n = n, p = p, delta = delta, cp = cp, minseglen = minseglen)
head(sim_results_Likelihood_Ratio)

# Run the Simulations for Nagao Trace
set.seed(1299)
sim_results_nagao <- simulations_nagao(n_simulations = n_simulations, n = n, p = p, delta = delta, cp = cp, minseglen = minseglen)
head(sim_results_nagao)

# Run the Simulations for Ryan and Killick
set.seed(1299)
sim_results_Ryan_Killick <- simulations_ryan_killick(n_simulations = n_simulations, n = n, p = p, delta = delta, cp = cp, minseglen = minseglen)
head(sim_results_Ryan_Killick)