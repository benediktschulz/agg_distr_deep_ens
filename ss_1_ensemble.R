## Simulation study: Script 1
# Generation of deep ensembles

#### Housekeeping ####
rm(list=ls())
gc()

#### Settings ####
# Load package
library(lubridate)
library(scoringRules)
library(dplyr)

# Path of simulated data
data_in_path <- "-"

# Path of deep ensemble forecasts
data_out_path <- "-"

# Load functions
source(file = paste0(getwd(), "/fn_cne_ss.R"))

#### Initialize ####
# Networks
nn_vec <- c("drn", "bqn", "hen")

# Models considered
scenario_vec <- 1:6

# Number of simulated runs
n_sim <- 50

# Size of network ensembles
n_ens <- 40

# # Choose number of cores
numCores <- parallel::detectCores()/2 - 1
if(pc == "simu"){ numCores <- 10 }

#### Functions ####
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

#### Loop over scenarios ####
# For-Loop over scenarios and simulations
for(i_scenario in scenario_vec){ for(i_sim in 1:n_sim){ 
  #### Get data ####
  # Load corresponding data
  load(file = paste0(data_in_path, "model", i_scenario, "_sim", i_sim, ".RData"))
  
  # Indeces of validation set
  if(i_scenario == 6){ i_valid <- 2501:3000 }
  else{ i_valid <- 5001:6000 }
  
  # Observations of validation set
  y_valid <- y_train[i_valid]
  
  #### Loop over network variants ####
  # For-Loop over network variants
  for(temp_nn in nn_vec){
    # Read out function
    fn_pp <- get(paste0(temp_nn, "_pp"))
    
    # Set seed (same for each network variant)
    set.seed(123 + 10*i_scenario + 100*i_sim)
    
    # For-Loop over ensemble member
    for(i_ens in 1:n_ens){
      # Postprocessing
      pred_nn <- fn_pp(train = X_train,
                       test = rbind(as.matrix(X_train[i_valid,]), X_test),
                       y_train = y_train,
                       y_test = c(y_valid, y_test),
                       i_valid = i_valid,
                       n_cores = numCores)
      
      # Transform ranks
      if(is.element("rank", names(pred_nn[["scores"]]))){
        pred_nn[["scores"]][["pit"]] <- fn_upit(ranks = pred_nn[["scores"]][["rank"]],
                                                max_rank = max(pred_nn[["scores"]][["rank"]]))
        
        # Omit ranks
        pred_nn[["scores"]][["rank"]] <- NULL
      }
      
      # Save ensemble member
      save(file = paste0(data_out_path, "model", i_scenario, 
                         "/model", i_scenario, "_", temp_nn, 
                         "_sim", i_sim, "_ens", i_ens, ".RData"),
           list = c("pred_nn", "y_valid", "y_test"))
    }
  }
}}

