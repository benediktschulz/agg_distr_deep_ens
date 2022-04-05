## Case study: Script 2
# Aggregation of deep ensembles
# !! Script based on the wind gust data from the case study !!

#### Housekeeping ####
rm(list=ls())
gc()

#### Settings ####
# Load package
library(lubridate)
library(scoringRules)
library(crch)
library(dplyr)
library(doParallel)
library(foreach)

# Path to auxiliary data
data_r_path <- "-"

# Path of deep ensembles
data_in_path <- "-"

# Path of aggregated forecasts
data_out_path <- "-"

# Load functions
source(file = paste0(get_wd(), "/fn_eval.R"))

#### Initialize ####
# Load some auxiliary information on the wind gust data (DATA NOT GIVEN)
load(file = paste0(data_r_path, "general_pars.RData"))

# Cores to use
numCores <- 20

# Network variants
nn_vec <- c("drn", "bqn", "hen")

# Aggregation methods
agg_meths_ls <- list()
agg_meths_ls[["drn"]] <- c("lp", "vi", "vi-q", "vi-w", "vi-a", "vi-aw")
agg_meths_ls[["bqn"]] <- 
  agg_meths_ls[["hen"]] <- c("lp", "vi", "vi-w", "vi-a", "vi-aw")

# lp -> Linear pool
# vi -> Vincentization
# vi-q -> Quasi-Vincentization
# vi-w -> (Quasi-)Vincentization with weight estimation
# vi-a -> (Quasi-)Vincentization with intercept estimation
# vi-aw -> (Quasi-)Vincentization with weight and intercept estimation

# Number of simulated ensemble members
n_sim <- 100

# Number of aggregation repetitions
n_rep <- 20

# Ensemble sizes to be combined
n_ens_vec <- seq(from = 2,
                 to = 40,
                 by = 2)

# Size of LP mixture samples
n_lp_samples <- 6*(n_ens + 1) - 1
# n_lp_samples <- 1000

# Size of BQN quantile samples
n_q_samples <- 6*(n_ens + 1) - 1

# Vector of initialization hours
temp_hr_vec <- 0

# Vector of forecast steps
temp_step_vec <- c(0, 6, 12, 18)

# Quantile levels for evaluation
q_levels <- seq(from = 1/(n_q_samples + 1),
                to = 1 - 1/(n_q_samples + 1),
                by = 1/(n_q_samples + 1))

#### Functions for weight estimation (and uPIT) ####
# CRPS of VI forecast in case of DRN
fn_vi_drn <- function(a, w, f_sum, y){
  ###-----------------------------------------------------------------------------
  ###Input
  #a.......Intercept (scalar)
  #w.......Equal weight (positive scalar)
  #f_sum...Sum of single DRN forecasts (n_train x 2 matrix)
  #y.......Observations (n_train vector)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...CRPS of VI forecast
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Penalty for non-positive weights
  if(w <= 0){ return(1e+6) }
  
  #### Calculation ####
  # Calculate weighted average
  f <- w*f_sum
  
  # Add intercept term (only to location)
  f[,1] <- a + f[,1]
  
  # Calculate CRPS
  res <- mean(crps_tlogis(y = y, 
                          location = f[,1], 
                          scale = f[,2],
                          lower = 0))
  
  # Output
  return(res)
}

# CRPS of VI forecast in case of BQN
fn_vi_bqn <- function(a, w, f_sum, y){
  ###-----------------------------------------------------------------------------
  ###Input
  #a.......Intercept (scalar)
  #w.......Weight (positive scalar)
  #f_sum...Sum of single BQN coefficients (n_train x 2 matrix)
  #y.......Observations (n_train vector)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...CRPS of VI forecast
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Penalty for non-positive weights
  if(w <= 0){ return(1e+6) }
  
  #### Calculation ####
  # Calculate weighted average
  alpha <- w*f_sum
  
  # Calculate quantiles
  q <- a + bern_quants(alpha = alpha,
                       q_levels = q_levels)
  
  # Calculate CRPS
  res <- mean(crps_sample(y = y, 
                          dat = q))
  
  # Output
  return(res)
}

# CRPS of VI forecast in case of HEN
fn_vi_hen <- function(a, w, f, bin_edges_sum, y){
  ###-----------------------------------------------------------------------------
  ###Input
  #a...............Intercept (scalar)
  #w...............Weight (positive scalar)
  #f_ls............List of bin probabilities (list of n_train x n_bins matrices)
  #bin_edges_sum...Accumulated bin edges ((n_bins + 1) vector)
  #y...............Observations (n_train vector)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...CRPS of VI forecast
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Penalty for non-positive weights
  if(w <= 0){ return(1e+6) }
  
  #### Calculation ####
  # Generate bin edges for each forecast
  bin_edges_f <- lapply(1:length(bin_edges_sum), function(i){
    a + w*bin_edges_sum[[i]] })
  
  # Calculate CRPS
  res <- mean(crps_hd(y = y,
                      f = f,
                      bin_edges = bin_edges_f))
  
  # Output
  return(res)
}

# Function to transform ranks to uPIT-values
fn_upit <- function(ranks, max_rank){
  ###-----------------------------------------------------------------------------
  ###Input
  #ranks......Ranks (n vector)
  #max_rank...Maximal rank (positive integer)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...uPIt-values (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Calculation ####
  # Transform to uPIT
  res <- ranks/max_rank - runif(n = length(ranks),
                                min = 0,
                                max = 1/max_rank)
  
  # Output
  return(res)
}

#### Initialize parallel computing ####
# Grid for parallel computing
grid_par <- expand.grid(1:length(nn_vec), 
                         sort(n_ens_vec, decreasing = TRUE),
                         temp_hr_vec, 
                         temp_step_vec,
                         1:n_rep)

#### Parallel-Function ####
# Function for parallel computing
fn_mc <- function(i_nn, temp_hr, temp_step, n_ens, i_rep){
  #### Input ####
  #i_nn........Network type indicator
  #temp_hr.....Initialization hour
  #temp_step...Lead time
  #n_ens.......Size of network ensemble
  #i_rep.......Repetition
  
  #### Initialization ####
  # Read out network type
  temp_nn <- nn_vec[i_nn]
  
  # Get aggregation methods
  agg_meths <- agg_meths_ls[[temp_nn]]
  
  # Create list for ensemble member forecasts on test and validation set
  f_ls <- f_valid_ls <- list()
  
  #### Get data ####
  # Draw indeces of ensemble members to aggregate (same for all aggregation methods)
  i_vec <- sort(sample(x = 1:n_sim,
                       size = n_ens,
                       replace = FALSE))
  
  # For-Loop over ensemble member
  for(i_sim in i_vec){
    # Load ensemble member
    load(file = paste0(data_in_path, temp_nn, 
                       "_hr", temp_hr, "_step", temp_step, "_sim", i_sim, ".RData"))
    
    # Get indices of validation and test set
    i_valid <- 1:length(y_valid)
    i_test <- length(y_valid) + 1:length(y_test)
    
    # Save forecasts
    if(temp_nn == "drn"){ 
      # Validation set
      f_valid_ls[[paste0("f", i_sim)]] <- pred_nn[["f"]][i_valid,]
      
      # Test set
      f_ls[[paste0("f", i_sim)]] <- pred_nn[["f"]][i_test,]
    }
    else if(temp_nn == "bqn"){
      # Validation set
      f_valid_ls[[paste0("f", i_sim)]] <- pred_nn[["f"]][i_valid,]
      f_valid_ls[[paste0("alpha", i_sim)]] <- pred_nn[["alpha"]][i_valid,]
      
      # Test set
      f_ls[[paste0("f", i_sim)]] <- pred_nn[["f"]][i_test,]
      f_ls[[paste0("alpha", i_sim)]] <- pred_nn[["alpha"]][i_test,]
    }
    else if(temp_nn == "hen"){ 
      # Validation set
      f_valid_ls[[paste0("f", i_sim)]] <- pred_nn[["f"]][i_valid,]
      f_valid_ls[[paste0("bin_edges", i_sim)]] <- pred_nn[["bin_edges_f"]]
      
      # Test set
      f_ls[[paste0("f", i_sim)]] <- pred_nn[["f"]][i_test,]
      f_ls[[paste0("bin_edges", i_sim)]] <- pred_nn[["bin_edges_f"]]
    }
  }
  
  # Average forecasts depending on method
  if(temp_nn == "drn"){
    # Average parameters
    f_sum <- apply(simplify2array(f_ls[paste0("f", i_vec)]), 1:2, sum)
    
    # Calculate sum of validation set forecasts
    f_valid_sum <- apply(simplify2array(f_valid_ls[paste0("f", i_vec)]), 1:2, sum)
  }
  else if(temp_nn == "bqn"){
    # Average forecasts
    f_sum <- apply(simplify2array(f_ls[paste0("f", i_vec)]), 
                   1:2, sum)
    
    # Average coefficients
    alpha_sum <- apply(simplify2array(f_ls[paste0("alpha", i_vec)]), 
                       1:2, sum)
    
    # Calculate sum of validation set forecasts
    f_valid_sum <- apply(simplify2array(f_valid_ls[paste0("alpha", i_vec)]), 1:2, sum)
  }
  else if(temp_nn == "hen"){
    # Matrix of all cumulative probabilities
    p_cum_valid <- do.call(cbind, lapply(i_vec, function(i_sim){
      # Cumulative sums of probabilities
      return(t(apply(f_valid_ls[[paste0("f", i_sim)]], 1, function(x){ pmin(1, cumsum(x)) }))) }))
    
    # Sort by row (round to reduce number of bins and increase increments)
    p_cum_sort_valid <- lapply(1:nrow(p_cum_valid), function(i){
      unique(c(0, round(sort(p_cum_valid[i,]), 4))) })
    
    # Generate corresponding bin probabilities
    f_valid <- lapply(1:length(p_cum_sort_valid), 
                      function(i) diff(p_cum_sort_valid[[i]]))
    
    # Generate bin edges for each forecast
    bin_edges_valid_sum <- lapply(1:length(p_cum_sort_valid), function(i){
      unique( rowSums( sapply(i_vec, function(i_sim){
        quant_hd(tau = p_cum_sort_valid[[i]],
                 probs = f_valid_ls[[paste0("f", i_sim)]][i,],
                 bin_edges = f_valid_ls[[paste0("bin_edges", i_sim)]]) }) )) })
    
    # Matrix of all cumulative probabilities
    p_cum <- do.call(cbind, lapply(i_vec, function(i_sim){
      # Cumulative sums of probabilities
      return(t(apply(f_ls[[paste0("f", i_sim)]], 1, 
                     function(x){ pmin(1, cumsum(x)) })))
    }))
    
    # Sort by row (round to reduce number of bins and increase increments)
    p_cum_sort <- lapply(1:nrow(p_cum), function(i){ unique(c(0, round(sort(p_cum[i,]), 4))) })
    
    # Probabilities of averaged bins
    f_bins <- lapply(1:length(p_cum_sort), function(i) diff(p_cum_sort[[i]]))
    
    # Accumulate bin edges
    bin_edges_sum <- lapply(1:length(p_cum_sort), function(i){
      unique(rowSums(sapply(i_vec, function(i_sim){
        quant_hd(tau = p_cum_sort[[i]],
                 probs = f_ls[[paste0("f", i_sim)]][i,],
                 bin_edges = f_ls[[paste0("bin_edges", i_sim)]])
      })))
    })
  }
  
  #### Aggregation ####
  # For-Loop over aggregation methods
  for(temp_agg in agg_meths){ 
    # Create list
    pred_agg <- list()
    
    # Different cases
    if((temp_nn == "drn") & (temp_agg == "lp")){
      # Function for mixture ensemble
      fn_apply <- function(i){
        # Sample individual distribution
        i_run <- sample(x = i_vec,
                        size = n_lp_samples,
                        replace = TRUE)
        
        # Get distributional parameters
        temp_f <- sapply(i_run, function(j){ f_ls[[paste0("f", j)]][i,] })
        
        # Draw from individual distributions
        res <- crch::qtlogis(p = runif(n = n_lp_samples),
                             location = temp_f[1,],
                             scale = temp_f[2,],
                             left = 0)
        
        # Output
        return(res)
      }
      
      # Simulate ensemble for mixture
      pred_agg[["f"]] <- t(sapply(1:length(y_test), fn_apply))
      
      # Calculate evaluation measure of simulated ensemble (mixture)
      pred_agg[["scores"]] <- fn_scores_ens(ens = pred_agg[["f"]],
                                            y = y_test)
      
      # Transform ranks to PIT
      pred_agg[["scores"]][["pit"]] <- fn_upit(ranks = pred_agg[["scores"]][["rank"]],
                                               max_rank = (n_lp_samples + 1))
      
      # No ranks
      pred_agg[["scores"]][["rank"]] <- NULL
    }
    else if((temp_nn == "drn") & (temp_agg == "vi")){
      # Calculate quantiles and average them
      pred_agg[["f"]] <- apply(simplify2array(lapply(i_vec, function(i){
        sapply(q_levels, function(p){
          crch::qtlogis(p = p, 
                        location = f_ls[[paste0("f", i)]][,1], 
                        scale = f_ls[[paste0("f", i)]][,2],
                        left = 0) })
      }) ), 1:2, mean)
      
      # Scores
      pred_agg[["scores"]] <- fn_scores_ens(ens = pred_agg[["f"]],
                                            y = y_test)
      
      # Transform ranks to PIT
      pred_agg[["scores"]][["pit"]] <- fn_upit(ranks = pred_agg[["scores"]][["rank"]],
                                               max_rank = (n_q_samples + 1))
      
      # No ranks
      pred_agg[["scores"]][["rank"]] <- NULL
    }
    else if((temp_nn == "drn") & (temp_agg == "vi-q")){
      # Average parameters
      pred_agg[["f"]] <- f_sum/n_ens
      
      # Scores
      pred_agg[["scores"]] <-  fn_scores_distr(f = pred_agg[["f"]],
                                               y = y_test,
                                               distr = "tlogis")
    }
    else if((temp_nn == "drn") & (temp_agg == "vi-w")){
      # Wrapper function
      fn_optim <- function(par){ fn_vi_drn(a = 0,
                                           w = par,
                                           f_sum = f_valid_sum,
                                           y = y_valid) }
      
      # Optimize
      est <- optim(par = 1/n_ens,
                   fn = fn_optim,
                   method = "Nelder-Mead")
      
      # Read out weights
      pred_agg[["w"]] <- est[["par"]]
      
      # Calculate optimally weighted VI
      pred_agg[["f"]] <- pred_agg[["w"]]*f_sum
      
      # Scores
      pred_agg[["scores"]] <-  fn_scores_distr(f = pred_agg[["f"]],
                                               y = y_test,
                                               distr = "tlogis")
    }
    else if((temp_nn == "drn") & (temp_agg == "vi-a")){
      # Wrapper function
      fn_optim <- function(par){ fn_vi_drn(a = par,
                                           w = 1/n_ens,
                                           f_sum = f_valid_sum,
                                           y = y_valid) }
      
      # Optimize
      est <- optim(par = 0,
                   fn = fn_optim,
                   method = "Nelder-Mead")
      
      # Read out intercept
      pred_agg[["a"]] <- est[["par"]]
      
      # Calculate equally weighted VI
      pred_agg[["f"]] <- f_sum/n_ens
      
      # Add intercept term (only to location)
      pred_agg[["f"]][,1] <- pred_agg[["a"]] + pred_agg[["f"]][,1]
      
      # Scores
      pred_agg[["scores"]] <-  fn_scores_distr(f = pred_agg[["f"]],
                                               y = y_test,
                                               distr = "tlogis")
    }
    else if((temp_nn == "drn") & (temp_agg == "vi-aw")){
      # Wrapper function
      fn_optim <- function(par){ fn_vi_drn(a = par[1],
                                           w = par[2],
                                           f_sum = f_valid_sum,
                                           y = y_valid) }
      
      # Optimize
      est <- optim(par = c(0, 1/n_ens),
                   fn = fn_optim,
                   method = "Nelder-Mead")
      
      # Read out intercept and weight
      pred_agg[["a"]] <- est[["par"]][1]
      pred_agg[["w"]] <- est[["par"]][2]
      
      # Calculate optimally weighted VI
      pred_agg[["f"]] <- pred_agg[["w"]]*f_sum
      
      # Add intercept term (only to location)
      pred_agg[["f"]][,1] <- pred_agg[["a"]] + pred_agg[["f"]][,1]
      
      # Scores
      pred_agg[["scores"]] <-  fn_scores_distr(f = pred_agg[["f"]],
                                               y = y_test,
                                               distr = "tlogis")
    }
    else if((temp_nn == "bqn") & (temp_agg == "lp")){
      # Function for mixture ensemble
      fn_apply <- function(i){
        # Sample individual distribution
        i_run <- sample(x = i_vec,
                        size = n_lp_samples,
                        replace = TRUE)
        
        # Draw from individual distributions
        res <- sapply(i_run, function(j){ 
          bern_quants(q_levels = runif(n = 1),
                      alpha = matrix(data = f_ls[[paste0("alpha", j)]][i,],
                                     nrow = 1)) })
        
        # Output
        return(res)
      }
      
      # Simulate ensemble for mixture
      pred_agg[["f"]] <- t(sapply(1:length(y_test), fn_apply))
      
      # Calculate evaluation measure of simulated ensemble (mixture)
      pred_agg[["scores"]] <- fn_scores_ens(ens = pred_agg[["f"]],
                                            y = y_test)
      
      # Transform ranks to PIT
      pred_agg[["scores"]][["pit"]] <- fn_upit(ranks = pred_agg[["scores"]][["rank"]],
                                               max_rank = (n_lp_samples + 1))
      
      # No ranks
      pred_agg[["scores"]][["rank"]] <- NULL
    }
    else if((temp_nn == "bqn") & (temp_agg == "vi")){
      # Average parameters
      pred_agg[["alpha"]] <- alpha_sum/n_ens
      pred_agg[["f"]] <- f_sum/n_ens
      
      # Scores
      pred_agg[["scores"]] <-  fn_scores_ens(ens = pred_agg[["f"]],
                                             y = y_test,
                                             skip_evals = c("e_me"))
      
      # Calculate bias of mean forecast (formula given)
      pred_agg[["scores"]][["e_me"]] <- rowMeans(pred_agg[["alpha"]]) - y_test
      
      # Transform ranks to PIT
      pred_agg[["scores"]][["pit"]] <- fn_upit(ranks = pred_agg[["scores"]][["rank"]],
                                               max_rank = (n_q_samples + 1))
      
      # No ranks
      pred_agg[["scores"]][["rank"]] <- NULL
    }
    else if((temp_nn == "bqn") & (temp_agg == "vi-w")){
      # Wrapper function
      fn_optim <- function(par){ 
        fn_vi_bqn(a = 0,
                  w = par,
                  f_sum = f_valid_sum,
                  y = y_valid) }
      
      # Optimize
      est <- optim(par = 1/n_ens,
                   fn = fn_optim,
                   method = "Nelder-Mead")
      
      # Read out weight
      pred_agg[["w"]] <- est[["par"]]
      
      # Optimally weighted parameters
      pred_agg[["alpha"]] <- pred_agg[["w"]]*alpha_sum
      pred_agg[["f"]] <- pred_agg[["w"]]*f_sum
      
      # Scores
      pred_agg[["scores"]] <-  fn_scores_ens(ens = pred_agg[["f"]],
                                             y = y_test,
                                             skip_evals = c("e_me"))
      
      # Calculate bias of mean forecast (formula given)
      pred_agg[["scores"]][["e_me"]] <- rowMeans(pred_agg[["alpha"]]) - y_test
      
      # Transform ranks to PIT
      pred_agg[["scores"]][["pit"]] <- fn_upit(ranks = pred_agg[["scores"]][["rank"]],
                                               max_rank = (n_q_samples + 1))
      
      # No ranks
      pred_agg[["scores"]][["rank"]] <- NULL
    }
    else if((temp_nn == "bqn") & (temp_agg == "vi-a")){
      # Wrapper function
      fn_optim <- function(par){ 
        fn_vi_bqn(a = par,
                  w = 1/n_ens,
                  f_sum = f_valid_sum,
                  y = y_valid) }
      
      # Optimize
      est <- optim(par = 0,
                   fn = fn_optim,
                   method = "Nelder-Mead")
      
      # Read out intercept
      pred_agg[["a"]] <- est[["par"]]
      
      # Optimally weighted parameters
      pred_agg[["alpha"]] <- alpha_sum/n_ens
      pred_agg[["f"]] <- pred_agg[["a"]] + f_sum/n_ens
      
      # Scores
      pred_agg[["scores"]] <-  fn_scores_ens(ens = pred_agg[["f"]],
                                             y = y_test,
                                             skip_evals = c("e_me"))
      
      # Calculate bias of mean forecast (formula given)
      pred_agg[["scores"]][["e_me"]] <- (pred_agg[["a"]] + rowMeans(pred_agg[["alpha"]])) - y_test
      
      # Transform ranks to PIT
      pred_agg[["scores"]][["pit"]] <- fn_upit(ranks = pred_agg[["scores"]][["rank"]],
                                               max_rank = (n_q_samples + 1))
      
      # No ranks
      pred_agg[["scores"]][["rank"]] <- NULL
    }
    else if((temp_nn == "bqn") & (temp_agg == "vi-aw")){
      # Wrapper function
      fn_optim <- function(par){ 
        fn_vi_bqn(a = par[1],
                  w = par[2],
                  f_sum = f_valid_sum,
                  y = y_valid) }
      
      # Optimize
      est <- optim(par = c(0, 1/n_ens),
                   fn = fn_optim,
                   method = "Nelder-Mead")
      
      # Read out intercept and weight
      pred_agg[["a"]] <- est[["par"]][1]
      pred_agg[["w"]] <- est[["par"]][2]
      
      # Optimally weighted parameters
      pred_agg[["alpha"]] <- pred_agg[["w"]]*alpha_sum
      pred_agg[["f"]] <- pred_agg[["a"]] + pred_agg[["w"]]*f_sum
      
      # Scores
      pred_agg[["scores"]] <-  fn_scores_ens(ens = pred_agg[["f"]],
                                             y = y_test,
                                             skip_evals = c("e_me"))
      
      # Calculate bias of mean forecast (formula given)
      pred_agg[["scores"]][["e_me"]] <- (pred_agg[["a"]] + rowMeans(pred_agg[["alpha"]])) - y_test
      
      # Transform ranks to PIT
      pred_agg[["scores"]][["pit"]] <- fn_upit(ranks = pred_agg[["scores"]][["rank"]],
                                               max_rank = (n_q_samples + 1))
      
      # No ranks
      pred_agg[["scores"]][["rank"]] <- NULL
    }
    else if((temp_nn == "hen") & (temp_agg == "lp")){
      # Average parameters
      pred_agg[["f"]] <- apply(simplify2array(f_ls[paste0("f", i_vec)]), 
                               1:2, mean)
      
      # Calculate scores
      pred_agg[["scores"]] <- fn_scores_hd(f = pred_agg[["f"]],
                                           y = y_test,
                                           bin_edges = f_ls[["bin_edges1"]])
    }
    else if((temp_nn == "hen") & (temp_agg == "vi")){
      # Generate corresponding bin probabilities
      pred_agg[["f"]] <- f_bins
      
      # Generate bin edges for each forecast
      pred_agg[["bin_edges_f"]] <- lapply(1:length(bin_edges_sum), function(i){
        bin_edges_sum[[i]]/n_ens })
      
      # Calculate scores
      pred_agg[["scores"]] <- fn_scores_hd(f = pred_agg[["f"]],
                                           y = y_test,
                                           bin_edges = pred_agg[["bin_edges_f"]])
    }
    else if((temp_nn == "hen") & (temp_agg == "vi-w")){
      # Wrapper function
      fn_optim <- function(par){ 
        fn_vi_hen(a = 0,
                  w = par,
                  f = f_valid,
                  bin_edges_sum = bin_edges_valid_sum,
                  y = y_valid) }
      
      # Optimize
      est <- optim(par = 1/n_ens,
                   fn = fn_optim,
                   method = "Nelder-Mead")
      
      # Read out weight
      pred_agg[["w"]] <- est[["par"]]
      
      # Generate corresponding bin probabilities
      pred_agg[["f"]] <- f_bins
      
      # Generate bin edges for each forecast
      pred_agg[["bin_edges_f"]] <- lapply(1:length(bin_edges_sum), function(i){
        pred_agg[["w"]]*bin_edges_sum[[i]] })
      
      # Calculate scores
      pred_agg[["scores"]] <- fn_scores_hd(f = pred_agg[["f"]],
                                           y = y_test,
                                           bin_edges = pred_agg[["bin_edges_f"]])
    }
    else if((temp_nn == "hen") & (temp_agg == "vi-a")){
      # Wrapper function
      fn_optim <- function(par){ 
        fn_vi_hen(a = par,
                  w = 1/n_ens,
                  f = f_valid,
                  bin_edges_sum = bin_edges_valid_sum,
                  y = y_valid) }
      
      # Optimize
      est <- optim(par = 0,
                   fn = fn_optim,
                   method = "Nelder-Mead")
      
      # Read out intercept
      pred_agg[["a"]] <- est[["par"]]
      
      # Generate corresponding bin probabilities
      pred_agg[["f"]] <- f_bins
      
      # Generate bin edges for each forecast
      pred_agg[["bin_edges_f"]] <- lapply(1:length(bin_edges_sum), function(i){
        pred_agg[["a"]] + bin_edges_sum[[i]]/n_ens })
      
      # Calculate scores
      pred_agg[["scores"]] <- fn_scores_hd(f = pred_agg[["f"]],
                                           y = y_test,
                                           bin_edges = pred_agg[["bin_edges_f"]])
    }
    else if((temp_nn == "hen") & (temp_agg == "vi-aw")){
      # Wrapper function
      fn_optim <- function(par){ 
        fn_vi_hen(a = par[1],
                  w = par[2],
                  f = f_valid,
                  bin_edges_sum = bin_edges_valid_sum,
                  y = y_valid) }
      
      # Optimize
      est <- optim(par = c(0, 1/n_ens),
                   fn = fn_optim,
                   method = "Nelder-Mead")
      
      # Read out intercept and weight
      pred_agg[["a"]] <- est[["par"]][1]
      pred_agg[["w"]] <- est[["par"]][2]
      
      # Generate corresponding bin probabilities
      pred_agg[["f"]] <- f_bins
      
      # Generate bin edges for each forecast
      pred_agg[["bin_edges_f"]] <- lapply(1:length(bin_edges_sum), function(i){
        pred_agg[["a"]] + pred_agg[["w"]]*bin_edges_sum[[i]] })
      
      # Calculate scores
      pred_agg[["scores"]] <- fn_scores_hd(f = pred_agg[["f"]],
                                           y = y_test,
                                           bin_edges = pred_agg[["bin_edges_f"]])
    }
    
    # Name of file
    file_name <- paste0(temp_nn, "_hr", temp_hr, "_step", temp_step, 
                        "_", temp_agg, "_ens", n_ens, "_rep", i_rep, ".RData")
    
    # Save aggregated forecasts and scores
    save(file = paste0(data_out_path, file_name),
         list = c("pred_agg", "i_vec"))
    
    # Delete and clean
    rm(pred_agg)
    gc()
  }
  
  # Delete
  if(temp_nn == "drn"){ rm(f_sum, f_valid_sum) }
  else if(temp_nn == "bqn"){ rm(f_sum, alpha_sum, f_valid_sum) }
  else if(temp_nn == "hen"){ rm(p_cum_valid, p_cum_sort_valid, f_valid, bin_edges_valid_sum, 
                                p_cum, p_cum_sort, f_bins, bin_edges_sum) }
  gc()
  
  # Dummy
  return(0)
}

#### Parallel-Loop ####
# Maximum number of cores
numCores <- min(numCores, nrow(grid_par))

# Use given number of cores
registerDoParallel(numCores)

# Call via for-each
foreach (i = 1:nrow(grid_par), 
         .errorhandling = c("remove"),
         .packages = c("lubridate", "dplyr", "scoringRules", "crch")) %dopar% {
           # Run function
           dummy <- fn_mc(i_nn = grid_par[i, 1],
                          n_ens = grid_par[i, 2],
                          temp_hr = grid_par[i, 3],
                          temp_step = grid_par[i, 4],
                          i_rep = grid_par[i, 5])
}