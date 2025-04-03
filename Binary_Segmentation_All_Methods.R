source("Binary_Segmentation_LR.R")
source("Binary_Segmentation_NT.R")
source("Binary_Segmentation_RK.R")

data = read.csv("Monthly_Prices.csv")[,4:22]

LR_results <- LR_binary_segmentation(data, threshold = 3, minseglen = 20)
LR_cpts <- LR_results$cpts
LR_cpts

NT_results <- NT_binary_segmentation(data, minseglen = 20)
NT_cpts <- NT_results$cpts
NT_cpts

RK_results <- RK_binary_segmentation(data, minseglen = 20)
RK_cpts <- RK_results$cpts
RK_cpts