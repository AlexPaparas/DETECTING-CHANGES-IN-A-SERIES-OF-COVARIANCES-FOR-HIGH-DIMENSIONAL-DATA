source("Single_CP_LR.R")
source("Single_CP_NT.R")
source("Single_CP_RK.R")

data = read.csv("Monthly_Prices.csv")[,4:22]

lr_results <- analyze_likelihood_ratio(data, minseglen = 20)
print(lr_results)
nagao_results <- analyze_nagao(data, minseglen = 20)
print(nagao_results)
rk_results <- analyze_rk(data, minseglen = 20)
print(rk_results)