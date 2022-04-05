## Case study: Script 1
# Generation of deep ensembles
# !! Script based on the wind gust data from the case study !!

#### Housekeeping ####
rm(list=ls())
gc()

#### Settings ####
# Load package
library(lubridate)
library(scoringRules)
library(dplyr)

# Path to auxiliary data
data_r_path <- "-"

# Path to wind gust forecast data
data_in_path <- "-"

# Path to storage of ensemble member
data_out_path <- "-"

# Load functions
source(file = paste0(get_wd(), "/fn_basic.R"))
source(file = paste0(get_wd(), "/fn_eval.R"))
source(file = paste0(get_wd(), "/fn_nn_cs.R"))

#### Initialize ####
# Load some auxiliary information on the wind gust data (DATA NOT GIVEN)
load(file = paste0(data_r_path, "general_pars.RData"))

# Networks
nn_vec <- c("drn", "bqn", "hen")

# Size of network ensembles
n_sim <- 100

# Choose number of cores
numCores <- parallel::detectCores()/2 - 1

# Vector of initialization hours
temp_hr_vec <- 0

# Vector of forecast steps
temp_step_vec <- c(0, 6, 12, 18)

# Location ID's
temp_loc_vec <- loc_vec

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

#### Predictors ####
# Predictors of wind gust ensemble (Load all)
pred_ens_stats <- c("ens_mean", "ens_sd")
pred_ens <- paste0("ens_", 1:n_ens)

# Temporal predictors
pred_tm <- c("yday")

# Spatial predictors (for global postprocessing)
pred_loc <- c("location", "lat", "lon", "altitude", "orog_diff",
              "loc_cover", "loc_bias")

# Additional variables
vars_add <- c(fc_var_vec[!grepl("_MS|_LS", fc_var_vec) & !grepl("VMAX_10M", fc_var_vec)],
              paste0("wind", c("_10M", "500", "700", "850", "950", "1000")))

# Additional predictors
pred_add <- c(paste0(vars_add, "_mean"))

#### Loop over hours and steps ####
# For-Loop over initialization hours and forecast steps
for(temp_hr in temp_hr_vec){ for(temp_step in temp_step_vec){
  #### Get data ####
  # Load corresponding data (DATA NOT GIVEN)
  load(file = paste0(data_in_path, "df_pp_data_test2016_init", 
                     temp_hr, "_step", temp_step,".RData"))
  
  # Prepare data (FUNCTION NOT GIVEN)
  df_train <- prep_data(df = df_train,
                        pred_ens = c(pred_ens_stats, pred_ens),
                        loc_path = data_r_path,
                        pred_tm = pred_tm,
                        pred_loc = pred_loc,
                        pred_add = pred_add,
                        sort_ens = any(grepl("ens_1", pred_ens, fixed = TRUE)),
                        return_ens = FALSE)
  df_test <- prep_data(df = df_test,
                       pred_ens = c(pred_ens_stats, pred_ens),
                       loc_path = data_r_path,
                       pred_tm = pred_tm,
                       pred_loc = pred_loc,
                       pred_add = pred_add,
                       sort_ens = any(grepl("ens_1", pred_ens, fixed = TRUE)))

  # Assign spatial predictors to test set (FUNCTION NOT GIVEN)
  df_test <- assign_spatial_preds(preds = pred_loc,
                                  df_train = df_train,
                                  df_test = df_test)
  
  # Indeces of validation set
  i_valid <- which(year(df_train$init_tm) == 2015)
  
  # Observations of validation set
  y_valid <- df_train[i_valid, "obs"]
  
  # Observations
  y_test <- df_test[["obs"]]
  
  #### Loop over network variants ####
  # For-Loop over network variants
  for(temp_nn in nn_vec){
    # Read out function
    fn_pp <- get(paste0(temp_nn, "_pp"))
    
    # General predictors for global post-processing
    if(temp_nn == "bqn"){ pred_vars <- c(pred_ens, pred_tm, pred_loc, pred_add) }
    else{ pred_vars <- c(pred_ens_stats, pred_tm, pred_loc, pred_add) }
    
    # Set seed (same for each network variant)
    set.seed(123)
    
    # For-Loop over ensemble member
    for(i_sim in 1:n_sim){
      # Postprocessing
      pred_nn <- fn_pp(train = df_train,
                       X = rbind(df_train[i_valid,], df_test),
                       i_valid = i_valid,
                       loc_id_vec = temp_loc_vec,
                       pred_vars = pred_vars,
                       n_cores = numCores,
                       scores_ens = FALSE,
                       scores_pp = TRUE)
      
      # Generate uPIT-values from ranks
      if(is.element("rank", names(pred_nn[["scores_pp"]]))){
        pred_nn[["scores_pp"]][["pit"]] <- fn_upit(ranks = pred_nn[["scores_pp"]][["rank"]],
                                                   max_rank = max(pred_nn[["scores_pp"]][["rank"]]))
        
        # Omit ranks
        pred_nn[["scores_pp"]][["rank"]] <- NULL
      }
      
      # Save ensemble member
      save(file = paste0(data_out_path, temp_nn, 
                         "_hr", temp_hr, "_step", temp_step, "_sim", i_sim, ".RData"),
           list = c("pred_nn", "y_valid", "y_test"))
    }
  }
}}

