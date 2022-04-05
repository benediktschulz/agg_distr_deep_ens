## Simulation study: Script 3
# Scores of deep ensembles and aggregated forecasts


#### Housekeeping ####
rm(list=ls())
gc()

#### Settings ####
# Load package
library(lubridate)
library(scoringRules)
library(dplyr)
library(RColorBrewer)

# Path of simulated data
data_sim_path <- "-"

# Path of deep ensemble data
data_ens_path <- "-"

# Path of aggregated deep ensembles
data_agg_path <- "-"

# Load functions
source(file = paste0(getwd(), "/fn_eval.R"))

#### Initialize ####
# Models considered
scenario_vec <- 1:6

# Number of simulations
n_sim <- 50

# Ensemble sizes to be combined
n_ens_vec <- seq(from = 2,
                 to = 40,
                 by = 2)

# (Maximum) number of network ensembles
n_rep <- max(n_ens_vec)

# Network types
nn_vec <- c("drn", "bqn", "hen")

# Aggregation methods
agg_meths_ls <- list()
agg_meths_ls[["drn"]] <- c("lp", "vi", "vi-w", "vi-a", "vi-aw")
agg_meths_ls[["bqn"]] <-
  agg_meths_ls[["hen"]] <- c("lp", "vi", "vi-w", "vi-a", "vi-aw")

# lp -> Linear pool
# vi -> Vincentization
# vi-w -> Vincentization with weight estimation
# vi-a -> Vincentization with intercept estimation
# vi-aw -> Vincentization with weight and intercept estimation

# To evaluate
sr_eval <- c("crps", "logs", "lgt", "cov", "mae", "me", "rmse")

# Skill scores
sr_skill <- c("crps", "logs", "mae", "rmse")

# Vector of column names
col_vec_pp <- c("model", "n_sim", "n_test", "n_train", "n_valid",
                "nn", "type", "n_ens", "n_rep", "a", "w", 
                sr_eval, paste0(sr_skill, "s"))

#### Create data frame ####
# Data frame for scores
df_scores <- as.data.frame(matrix(ncol = length(col_vec_pp),
                                  nrow = 0))
colnames(df_scores) <- col_vec_pp

# For-Loop over network types
for(temp_nn in nn_vec){ 
  # Get aggregation methods
  agg_meths <- agg_meths_ls[[temp_nn]]

  # For-Loop over scenarios and simulations
  for(i_scenario in scenario_vec){ for(i_sim in 1:n_sim){
    # Optimal forecast for first network variant
    if(temp_nn == nn_vec[1]){
      # Load data
      load(paste0(data_sim_path, "model", i_scenario, "_sim", i_sim, ".RData"))
      
      # Create data frame
      res <- as.data.frame(matrix(ncol = length(col_vec_pp),
                                  nrow = 1))
      colnames(res) <- col_vec_pp
      
      # Make entry for optimal forecast
      res[["model"]] <- i_scenario
      res[["n_sim"]] <- i_sim
      res[["nn"]] <- "ref"
      res[["type"]] <- "ref"
      
      # For-Loop over evaluation measures
      for(temp_sr in sr_eval){
        # Depending on measure
        if(temp_sr == "mae"){ 
          res[[temp_sr]] <- mean(abs(scores_opt[["e_md"]])) }
        else if(temp_sr == "me"){ 
          res[[temp_sr]] <- mean(scores_opt[["e_md"]]) }
        else if(temp_sr == "rmse"){
          res[[temp_sr]] <- sqrt(mean((scores_opt[["e_me"]])^2)) }
        else if(temp_sr == "cov"){ 
          res[[temp_sr]] <- mean(fn_cover(scores_opt[["pit"]])) }
        else{ res[[temp_sr]] <- mean(scores_opt[[temp_sr]]) }
      }
      
      # Append to data frame
      df_scores <- rbind(df_scores, res)
    }
    
    # For-Loop over repetitions
    for(i_rep in 1:n_rep){
      # Create data frame
      res <- as.data.frame(matrix(ncol = length(col_vec_pp),
                                  nrow = 1))
      colnames(res) <- col_vec_pp
      
      # Write in data frame
      res[["model"]] <- i_scenario
      res[["n_sim"]] <- i_sim
      res[["nn"]] <- temp_nn
      res[["type"]] <- "ind"
      res[["n_ens"]] <- 1
      res[["n_rep"]] <- i_rep
      
      # Load ensemble member
      load(file = paste0(data_ens_path, "model", i_scenario, 
                         "/model", i_scenario, "_", temp_nn, 
                         "_sim", i_sim, "_ens", i_rep, ".RData"))
      
      # Get set size of first repetition
      for(temp in c("n_train", "n_valid", "n_test")){
        # Save set sizes
        if(i_rep == 1){ 
          # Read out set sizes
          assign(x = temp,
                 value = pred_nn[[temp]]) 
          
          # Calculate actual test set size
          if(temp == "n_test"){ n_test <- n_test - n_valid }
        } 
        
        # Read out set sizes
        res[[temp]] <- get(temp)
      }
      
      # Cut validation data from scores
      pred_nn[["scores"]] <- pred_nn[["scores"]][-(1:pred_nn[["n_valid"]]),]
      
      # For-Loop over evaluation measures
      for(temp_sr in sr_eval){
        # Depending on measure
        if(temp_sr == "mae"){ 
          res[[temp_sr]] <- mean(abs(pred_nn[["scores"]][["e_md"]])) }
        else if(temp_sr == "me"){ 
          res[[temp_sr]] <- mean(pred_nn[["scores"]][["e_md"]]) }
        else if(temp_sr == "rmse"){
          res[[temp_sr]] <- sqrt(mean((pred_nn[["scores"]][["e_me"]])^2)) }
        else if(temp_sr == "cov"){ 
          res[[temp_sr]] <- mean(fn_cover(pred_nn[["scores"]][["pit"]])) }
        else{ res[[temp_sr]] <- mean(pred_nn[["scores"]][[temp_sr]]) }
      }
      
      # Append to data frame
      df_scores <- rbind(df_scores, res)
    }
    
    # For-Loop over aggregation methods
    for(temp_agg in agg_meths){
      # For-Loop over number of aggregated members
      for(i_ens in n_ens_vec){
        # Create data frame
        res <- as.data.frame(matrix(ncol = length(col_vec_pp),
                                    nrow = 1))
        colnames(res) <- col_vec_pp
        
        # Write in data frame
        res[["model"]] <- i_scenario
        res[["n_sim"]] <- i_sim
        res[["nn"]] <- temp_nn
        res[["type"]] <- temp_agg
        res[["n_ens"]] <- i_ens
        res[["n_rep"]] <- 0
        
        # Get set sizes
        for(temp in c("n_train", "n_valid", "n_test")){
          # Read out set sizes
          res[[temp]] <- get(temp)
        }
        
        # Load aggregated forecasts
        load(file = paste0(data_agg_path, "model", i_scenario, 
                           "/model", i_scenario, "_", temp_nn, 
                           "_sim", i_sim, "_", temp_agg, "_ens", i_ens, ".RData"))
        
        # # Estimated weights and intercept
        if(temp_agg == "vi-w"){ res[["w"]] <- pred_agg[["w"]] }
        else if(temp_agg == "vi-a"){ res[["a"]] <- pred_agg[["a"]] }
        else if(temp_agg == "vi-aw"){
          res[["a"]] <- pred_agg[["a"]]
          res[["w"]] <- pred_agg[["w"]]
        }
        
        # For-Loop over evaluation measures
        for(temp_sr in sr_eval){
          # Depending on measure
          if(temp_sr == "mae"){ 
            res[[temp_sr]] <- mean(abs(pred_agg[["scores"]][["e_md"]])) }
          else if(temp_sr == "me"){ 
            res[[temp_sr]] <- mean(pred_agg[["scores"]][["e_md"]]) }
          else if(temp_sr == "rmse"){
            res[[temp_sr]] <- sqrt(mean((pred_agg[["scores"]][["e_me"]])^2)) }
          else if(temp_sr == "cov"){ 
            res[[temp_sr]] <- mean(fn_cover(pred_agg[["scores"]][["pit"]])) }
          else{ res[[temp_sr]] <- mean(pred_agg[["scores"]][[temp_sr]]) }
        }
        
        # For-Loop over skill scores
        for(temp_sr in sr_skill){
          # Reference is given by mean score of network ensemble members
          s_ref <- mean(subset(df_scores, (model == i_scenario) &
                                 (n_sim == i_sim) &
                                 (nn == temp_nn) &
                                 (type == "ind") &
                                 (n_rep <= i_ens))[[temp_sr]])
          
          # Score of optimal forecast
          s_opt <- subset(df_scores, (model == i_scenario) & 
                            (n_sim == i_sim) & 
                            (nn == "ref"))[[temp_sr]]
          
          # Calculate skill
          res[[paste0(temp_sr, "s")]] <- (s_ref - res[[temp_sr]])/(s_ref - s_opt)
        }
        
        # Append to data frame
        df_scores <- rbind(df_scores, res)
      }
    }
  }}
}

#### Save ####
save(file = paste0(getwd(), "/data/eval_ss.RData"),
     list = c("df_scores"))
