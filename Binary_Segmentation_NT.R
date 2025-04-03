# Binary Segmentation for Nagao's Trace Test (NT Method) -------------------------
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
                   2*(y1 + y2)*(y1^2 + y2^2 - y1*y2) + y1*y2*(2*y1 - y2)*(2*y2 - y1))
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

# Helper function to analyze a single segment for NT method
analyze_nagao_segment <- function(Y, minseglen, alpha) {
  p <- ncol(Y)
  n <- nrow(Y)
  
  # Candidate change points within this segment
  t_0 <- seq(minseglen + 1, n - minseglen)
  
  # Pre-allocate results
  results <- matrix(nrow = length(t_0), ncol = 3)
  colnames(results) <- c("candidate_cp", "standardized_stat", "raw_stat")
  
  # Compute statistics for each candidate point
  for(i in seq_along(t_0)) {
    t <- t_0[i]
    components <- compute_nagao_components(
      data_before = Y[1:t, , drop = FALSE],
      data_after = Y[(t+1):n, , drop = FALSE],
      p = p, n = n, t = t
    )
    
    results[i, ] <- c(
      t,
      compute_nagao_statistic(components, p),
      sum(((components$values - 1)^2)/((components$eta*components$values + (1-components$eta))^2))
    )
  }
  return(results)
}

# Binary segmentation function for NT method
bin.seg.NT <- function(data, parent, threshold, minseglen, alpha = 0.05) {
  n <- nrow(data)
  p <- ncol(data)
  
  # Only proceed if segment is long enough
  if(n > 2*(minseglen + p)) {
    # Compute NT statistics for all candidate points in this segment
    results <- analyze_nagao_segment(data, minseglen, alpha)
    
    # Find the maximum statistic and its location
    max_stat <- max(results[, "standardized_stat"], na.rm = TRUE)
    candidate <- results[which.max(results[, "standardized_stat"]), "candidate_cp"]
    
    if(max_stat > threshold) {
      # Split the data at the detected change point
      data1 <- data[1:candidate, , drop = FALSE]
      data2 <- data[(candidate+1):n, , drop = FALSE]
      
      # Update parent boundaries
      parent1 <- c(parent[1], parent[1] + candidate)
      parent2 <- c(parent[1] + candidate, parent[2])
      
      # Recursively apply binary segmentation
      cpt1 <- bin.seg.NT(data1, parent1, threshold, minseglen, alpha)
      cpt2 <- bin.seg.NT(data2, parent2, threshold, minseglen, alpha)
      
      return(list(
        c(parent, max_stat, parent[1] + candidate),  # Current change point info
        cpt1,  # Left child results
        cpt2   # Right child results
      ))
    } else {
      return(list(c(parent, max_stat, NA), 0))  # No significant change point
    }
  } else {
    return(list(c(parent, NA, NA), 0))  # Segment too short
  }
}

# Binary Segmentation Output Processing Function ---------------------------
bin.seg.to.cpt <- function(result, threshold) {
  cpts <- matrix(ncol = 4, nrow = 0)
  colnames(cpts) <- c("start", "end", "statistic", "cp")
  cpts.out <- matrix(ncol = 4, nrow = 0)
  thresholds <- c()
  
  # Recursive processing function
  process_segment <- function(segment) {
    if(!is.na(segment[[1]][3])) {  # If there's a candidate CP
      if(segment[[1]][3] > threshold) {
        # This is a significant CP - add to results and process children
        cpts <<- rbind(cpts, segment[[1]])
        thresholds <<- c(thresholds, segment[[1]][3])
        
        # Process left and right children if they exist
        if(length(segment) > 2) {
          process_segment(segment[[2]])
          process_segment(segment[[3]])
        }
      } else {
        # Not significant - add to non-CP segments
        cpts.out <<- rbind(cpts.out, segment[[1]])
      }
    } else {
      # No CP found in this segment
      cpts.out <<- rbind(cpts.out, segment[[1]])
    }
  }
  
  process_segment(result)
  
  # Adjust segment start points (R is 1-indexed)
  if(nrow(cpts.out) > 0) {
    cpts.out[,1] <- cpts.out[,1] + 1
  }
  
  return(list(
    cpts = if(nrow(cpts) > 0) cpts[,4] else numeric(0),
    segments = cpts.out,
    thresholds = thresholds
  ))
}

# Modified filtering function that keeps the stronger signal
filter_change_points <- function(cpts, stats, min_distance = 30) {
  # Combine into data frame and sort by position
  cp_df <- data.frame(position = cpts, statistic = stats)
  cp_df <- cp_df[order(cp_df$position), ]
  
  if(nrow(cp_df) == 0) return(numeric(0))
  
  filtered <- cp_df[1, , drop = FALSE]
  
  for(i in 2:nrow(cp_df)) {
    last_pos <- filtered$position[nrow(filtered)]
    current_pos <- cp_df$position[i]
    
    if(current_pos - last_pos < min_distance) {
      # If too close, keep the one with higher statistic
      if(cp_df$statistic[i] > filtered$statistic[nrow(filtered)]) {
        filtered[nrow(filtered), ] <- cp_df[i, ]
      }
    } else {
      filtered <- rbind(filtered, cp_df[i, ])
    }
  }
  
  return(filtered$position)
}

# Bonferroni correction function
bonferoni <- function(n, alpha = 0.05) {
  qnorm(1 - alpha/n)
}

# Wrapper function for NT binary segmentation
NT_binary_segmentation <- function(data, threshold = NULL, minseglen = NULL, alpha = 0.05, bonferroni = TRUE) {
  if(is.null(minseglen)) minseglen <- max(4*ncol(data), 30)
  
  # Calculate threshold with or without Bonferroni correction
  if(is.null(threshold)) {
    if(bonferroni) {
      # Estimate number of tests (approximate for binary segmentation)
      n_tests <- ceiling(log2(nrow(data)/minseglen)) * (nrow(data)/minseglen)
      threshold <- bonferoni(n_tests, alpha)
    } else {
      threshold <- qnorm(1 - alpha)
    }
  }
  
  initial_parent <- c(1, nrow(data))
  result <- bin.seg.NT(data, initial_parent, threshold, minseglen, alpha)
  
  # Process results
  processed <- bin.seg.to.cpt(result, threshold)
  
  # Apply filtering that considers statistic strength
  if(length(processed$cpts) > 0) {
    filtered <- filter_change_points(processed$cpts, processed$thresholds, 30)
    processed$cpts <- filtered
  }
  
  # Add threshold information to output
  processed$threshold_used <- threshold
  processed$threshold_method <- ifelse(bonferroni, "Bonferroni", "Standard")
  
  return(processed)
}

# Example:
#data <- read.csv("Monthly_Prices.csv")[,4:22]
#result <- NT_binary_segmentation(data, minseglen = 20)
#filtered_cpts <- result$cpts
#filtered_cpts