## Simulation study: Script 0
# Simulate underlying data

#### Housekeeping ####
rm(list=ls())
gc()

#### Settings ####
# Path for simulated data
data_out_path <- "-"

# Load functions
source(file = paste0(getwd(), "/fn_eval.R"))

#### Initialize ####
# Models considered
scenario_vec <- 1:6

# Number of simulations
n_sim <- 50

#### Generate data ####
# For-Loop over simulations and scenarios
for(i_scenario in scenario_vec){ for(i_sim in 1:n_sim){
  #### Initiation ####
  # Set seed
  set.seed(123*i_sim)
  
  # Size of training sets
  if(i_scenario == 6){ n_train <- 3000 }
  else{ n_train <- 6000 }
  
  # Size of test sets
  n_test <- 1e+4
  
  # Indices for training and test set
  i_train <- 1:n_train
  i_test <- n_train + 1:n_test
  
  # Number of predictors for each model
  if(i_scenario == 3){ n_preds <- 1 }
  else{ n_preds <- 5 }
  
  #### Generate data ####
  # Differentiate models
  if(i_scenario == 1){
    # Coefficients
    beta_1 <- rnorm(n = n_preds,
                    mean = 0,
                    sd = 1)
    beta_2 <- rnorm(n = n_preds,
                    mean = 0,
                    sd = 0.45)
    
    # Predictors
    X <- matrix(data = rnorm(n = (n_train + n_test)*n_preds,
                             mean = 0,
                             sd = 1),
                      nrow = (n_train + n_test),
                      ncol = n_preds)
    
    # Observational error
    eps <- rnorm(n = n_train + n_test,
                 mean = 0,
                 sd = 1)
    
    # Calculate observations
    y <- as.vector(X %*% beta_1 + exp(X %*% beta_2)*eps)
  }
  else if(is.element(i_scenario, c(2, 5, 6))){
    # Predictors
    if(i_scenario == 5){
      # Generate covariance matrix
      sigma <- matrix(data = 1,
                      nrow = n_preds,
                      ncol = n_preds)
      
      # For-Loop over rows and columns
      for(i in 1:nrow(sigma)){ for(j in 1:ncol(sigma)){
        sigma[i,j] <- 0.5^(abs(i - j)) 
      }}
      
      # Draw correlated variables
      X <- MultiRNG::draw.d.variate.uniform(
        no.row = n_train + n_test,
        d = n_preds,
        cov.mat = sigma
      )
    }
    else{
      # Draw iid uniform random variables
      X <- matrix(data = runif(n = (n_train + n_test)*n_preds,
                               min = 0,
                               max = 1),
                  nrow = (n_train + n_test),
                  ncol = n_preds)
    }
    
    # Bernoulli variable
    bn <- rbinom(n = n_train + n_test,
                 size = 1,
                 prob = 0.5)
    
    # Observational error
    eps_1 <- rnorm(n = n_train + n_test,
                   mean = 0,
                   sd = 1.5) # var = 2.25
    eps_2 <- rnorm(n = n_train + n_test,
                   mean = 0,
                   sd = 1)
    
    # Calculate observations
    y <- bn*( 10*sin(2*pi*X[,1]*X[,2]) + 10*X[,4] + eps_1 ) +
      (1 - bn)*( 20*(X[,3] - 0.5)^2 + 5*X[,5] + eps_2 )
  }
  else if(i_scenario == 3){
    # Predictors
    X <- matrix(data = runif(n = (n_train + n_test)*n_preds,
                             min = 0,
                             max = 10),
                nrow = (n_train + n_test),
                ncol = n_preds)
    
    # Bernoulli variable
    bn <- rbinom(n = n_train + n_test,
                 size = 1,
                 prob = 0.5)
    
    # Observational error
    eps_1 <- rnorm(n = n_train + n_test,
                   mean = 0,
                   sd = 0.3) # var = 0.09
    eps_2 <- rnorm(n = n_train + n_test,
                   mean = 0,
                   sd = 0.8) # var = 0.64
    
    # Calculate observations
    y <- bn*( sin(X[,1]) + eps_1 ) +
      (1 - bn)*( 2*sin(1.5*X[,1] + 1) + eps_2 )
  }
  else if(i_scenario == 4){
    # Predictors
    X <- matrix(data = runif(n = (n_train + n_test)*n_preds,
                             min = 0,
                             max = 1),
                nrow = (n_train + n_test),
                ncol = n_preds)
    
    # Observational error (skew normal)
    eps <- sn::rsn(n = n_train + n_test,
                   xi = 0, # location
                   omega = 1, # scale
                   alpha = -5) # skew
    
    # Calculate observations
    y <- 10*sin(2*pi*X[,1]*X[,2]) + 20*(X[,3] - 0.5)^2 + 10*X[,4] + 5*X[,5] + eps
  }

  #### Data partition ####
  # Split in training and testing
  X_train <- as.matrix(X[i_train,])
  X_test <- as.matrix(X[i_test,])
  y_train <- y[i_train]
  y_test <- y[i_test]
  
  #### Optimal forecast ####
  # Generate matrix for parameter forecasts / sample
  if(i_scenario == 4){
    # Number of samples to draw
    n_sample <- 1e+3
  }
  else{ 
    # Normal distribution
    f_opt <- matrix(ncol = 2,
                    nrow = length(y_test))
    colnames(f_opt) <- c("loc", "scale")
  }
  
  # Differentiate scenarions
  if(i_scenario == 1){
    # Location parameter
    f_opt[,1] <- as.vector(X_test %*% beta_1)
    
    # Scale parameter (standard deviation)
    f_opt[,2] <- as.vector(exp(X_test %*% beta_2))
  }
  else if(is.element(i_scenario, c(2, 5, 6))){
    ## Assumption Bernoulli variable is known (elsewise multimodal)
    # Location
    f_opt[,1] <- bn[i_test]*( 10*sin(2*pi*X_test[,1]*X_test[,2]) + 10*X_test[,4] ) +
      (1 - bn[i_test])*( 20*(X_test[,3] - 0.5)^2 + 5*X_test[,5] )
    
    # Scale parameter (standard deviation)
    f_opt[,2] <- bn[i_test]*1.5 + (1 - bn[i_test])*1
  }
  else if(i_scenario == 3){
    ## Assumption Bernoulli variable is known (elsewise multimodal)
    # Location
    f_opt[,1] <- bn[i_test]*( sin(X_test[,1]) ) +
      (1 - bn[i_test])*( 2*sin(1.5*X_test[,1] + 1) )
    
    # Scale parameter (standard deviation)
    f_opt[,2] <- bn[i_test]*0.3 + (1 - bn[i_test])*0.8
  }
  else if(i_scenario == 4){
    # Draw samples from a skewed normal
    f_opt <- t(apply(X_test, 1, function(x){
      sn::rsn(n = n_sample,
              xi = 10*sin(2*pi*x[1]*x[2]) + 
                20*(x[3] - 0.5)^2 + 10*x[4] + 5*x[5], # location
              omega = 1, # scale
              alpha = -5) # skew (?)
    }))
  }
  
  # Optimal scores
  if(i_scenario == 4){
    # Number of samples to draw
    scores_opt <- fn_scores_ens(ens = f_opt,
                                y = y_test)
    
    # Transform ranks to uPIT
    scores_opt[["pit"]] <- 
      scores_opt[["rank"]]/(n_sample + 1) - runif(n = n_test,
                                                  min = 0,
                                                  max = 1/(n_sample + 1))
    
    # Omit ranks
    scores_opt[["rank"]] <- NULL
  }
  else{ 
    # Normal distribution
    scores_opt <- fn_scores_distr(f = f_opt,
                                  y = y_test,
                                  distr = "norm")
  }
  
  #### Save data ####
  # Save ensemble member
  save(file = paste0(data_out_path, "model", i_scenario, "_sim", i_sim, ".RData"),
       list = c("X_train", "y_train", 
                "X_test", "y_test",
                "f_opt", "scores_opt"))
}}