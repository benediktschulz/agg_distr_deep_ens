## Case study: Script 2
# Scores of deep ensembles and aggregated forecasts
# !! Script based on the wind gust data from the case study !!


#### Housekeeping ####
rm(list=ls())
gc()

#### Settings ####
# Load package
library(lubridate)
library(scoringRules)
library(dplyr)
library(RColorBrewer)

# Path to auxiliary data
data_r_path <- "-"

# Path of deep ensembles
data_ens_path <- "-"

# Path of aggregated forecasts
data_agg_path <- "-"

# Load functions
source(file = paste0(getwd(), "/fn_eval.R"))

#### Initialize ####
# Load some auxiliary information on the wind gust data (DATA NOT GIVEN)
load(file = paste0(data_r_path, "general_pars.RData"))

# Number of simulated ensemble members
n_sim <- 100

# Number of aggregation repetitions
n_rep <- 20

# Ensemble sizes considered for aggregation
n_ens_vec <- seq(from = 2,
                 to = 40,
                 by = 2)

# Network types
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

# Vector of initialization hours
temp_hr_vec <- 0

# Vector of forecast steps
temp_step_vec <- c(0, 6, 12, 18)

# To evaluate
sr_eval <- c("crps", "logs", "lgt", "cov", "mae", "me", "rmse")

# Skill scores
sr_skill <- c("crps", "logs", "mae", "rmse")

# Vector of column names
col_vec_pp <- c("init_hr", "fc_step", "n_test", "n_train", "n_valid",
                "nn", "type", "n_ens", "i_sim", "i_rep", "a", "w", 
                sr_eval, paste0(sr_eval, "_ref"), paste0(sr_eval, "_refmin"),
                paste0(sr_skill, "s"), paste0(sr_skill, "s_min"))

#### Create data frame ####
# Data frame for scores
df_scores <- as.data.frame(matrix(ncol = length(col_vec_pp),
                                  nrow = 0))
colnames(df_scores) <- col_vec_pp

# For-Loop over network types
for(temp_nn in nn_vec){ 
  # Get aggregation methods
  agg_meths <- agg_meths_ls[[temp_nn]]

  # For-Loop over initialization hours and lead times
  for(temp_hr in temp_hr_vec){ for(temp_step in temp_step_vec){
    # For-Loop over simulated ensemble members
    for(i_sim in 1:n_sim){
      # Create data frame
      res <- as.data.frame(matrix(ncol = length(col_vec_pp),
                                  nrow = 1))
      colnames(res) <- col_vec_pp
      
      # Write in data frame
      res[["init_hr"]] <- temp_hr
      res[["fc_step"]] <- temp_step
      res[["nn"]] <- temp_nn
      res[["type"]] <- "ind"
      res[["n_ens"]] <- 1
      res[["i_sim"]] <- i_sim
      res[["i_rep"]] <- 0
      
      # Load ensemble member
      load(file = paste0(data_ens_path, temp_nn, 
                         "_hr", temp_hr, "_step", temp_step, "_sim", i_sim, ".RData"))
      
      # Get set size of first repetition
      for(temp in c("n_train", "n_valid", "n_test")){
        # Save set sizes
        if(i_sim == 1){ 
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
      pred_nn[["scores_pp"]] <- pred_nn[["scores_pp"]][-(1:pred_nn[["n_valid"]]),]
      
      # For-Loop over evaluation measures
      for(temp_sr in sr_eval){
        # Depending on measure
        if(temp_sr == "mae"){ 
          res[[temp_sr]] <- mean(abs(pred_nn[["scores_pp"]][["e_md"]])) }
        else if(temp_sr == "me"){ 
          res[[temp_sr]] <- mean(pred_nn[["scores_pp"]][["e_md"]]) }
        else if(temp_sr == "rmse"){
          res[[temp_sr]] <- sqrt(mean((pred_nn[["scores_pp"]][["e_me"]])^2)) }
        else if(temp_sr == "cov"){ 
          res[[temp_sr]] <- mean(fn_cover(pred_nn[["scores_pp"]][["pit"]])) }
        else{ res[[temp_sr]] <- mean(pred_nn[["scores_pp"]][[temp_sr]]) }
      }
      
      # Append to data frame
      df_scores <- rbind(df_scores, res)
    }
    
    # For-Loop over aggregation methods
    for(temp_agg in agg_meths){
      # For-Loop over number of aggregated members and repetition
      for(n_ens in n_ens_vec){ for(i_rep in 1:n_rep){
        # Create data frame
        res <- as.data.frame(matrix(ncol = length(col_vec_pp),
                                    nrow = 1))
        colnames(res) <- col_vec_pp
        
        # Write in data frame
        res[["init_hr"]] <- temp_hr
        res[["fc_step"]] <- temp_step
        res[["nn"]] <- temp_nn
        res[["type"]] <- temp_agg
        res[["n_ens"]] <- n_ens
        res[["i_rep"]] <- i_rep
        res[["i_sim"]] <- 0
        
        # Get set sizes
        for(temp in c("n_train", "n_valid", "n_test")){
          # Read out set sizes
          res[[temp]] <- get(temp)
        }
        
        # Load aggregated forecasts
        load(file = paste0(data_agg_path, temp_nn, "_hr", temp_hr, "_step", temp_step, 
                           "_", temp_agg, "_ens", n_ens, "_rep", i_rep, ".RData"))
        
        # # Estimated weights and intercept
        if(temp_agg == "vi-w"){ res[["w"]] <- pred_agg[["w"]] }
        else if(temp_agg == "vi-a"){ res[["a"]] <- pred_agg[["a"]] }
        else if(temp_agg == "vi-aw"){
          res[["a"]] <- pred_agg[["a"]]
          res[["w"]] <- pred_agg[["w"]]
        }
        
        # Get subset of ensemble members for aggregation repetition
        df_rep <- subset(df_scores, (type == "ind") &
                           (init_hr == temp_hr) &
                           (fc_step == temp_step) &
                           (nn == temp_nn) &
                           is.element(i_sim, i_vec))
        
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
          
          # Reference scores mean
          if(temp_sr == "rmse"){
            res[[paste0(temp_sr, "_ref")]] <- sqrt(mean((df_rep[[temp_sr]])^2)) }
          else{ res[[paste0(temp_sr, "_ref")]] <- mean(df_rep[[temp_sr]]) }
          
          # Reference scores minimum
          if(temp_sr == "me"){
            res[[paste0(temp_sr, "_refmin")]] <- df_rep[[temp_sr]][which.min(abs(df_rep[[temp_sr]]))] }
          else if(temp_sr == "cov"){
            res[[paste0(temp_sr, "_refmin")]] <- df_rep[[temp_sr]][which.min(abs(df_rep[[temp_sr]] - 19/21*100))] }
          else{ res[[paste0(temp_sr, "_refmin")]] <- min(df_rep[[temp_sr]]) }
        }
        
        # For-Loop over skill scores
        for(temp_sr in sr_skill){
          # Calculate skill w.r.t. mean score
          res[[paste0(temp_sr, "s")]] <- 1 - (res[[temp_sr]]/res[[paste0(temp_sr, "_ref")]])
          
          # Calculate skill w.r.t. minimum score
          res[[paste0(temp_sr, "s_min")]] <- 1 - (res[[temp_sr]]/res[[paste0(temp_sr, "_refmin")]])
        }
        
        # Append to data frame
        df_scores <- rbind(df_scores, res)
      }}
    }
  }}
}

#### Save data ####
# Save
save(file = paste0(getwd(), "/data/eval_cs.RData"),
     list = c("df_scores"))