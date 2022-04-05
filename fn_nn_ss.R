## Function file
# Network variants for the simulation study

#### Import ####
# Import basic functions
source(paste0(getwd(), "/fn_basic.R"))
source(paste0(getwd(), "/fn_eval.R"))

#### DRN ####
## Distribution: Normal distribution
## Estimation: CRPS

# Function for estimation and prediction
drn_pp <- function(train, test, y_train, y_test, 
                   i_valid = NULL, nn_ls = list(), n_cores = NULL){
  ###-----------------------------------------------------------------------------
  ###Input
  #train/test......Training/Test data including predictors (n_train/test x n_preds data.frame)
  #y_train/test....Training/Test observations (n_train/test vector)
  #i_valid.........Indices of validation data (n_valid vector)
  #................Default: NULL -> Use fraction of training data
  #nn_ls...........List that may contain the following variables:
  #...lr_adam......Learning rate in Adam-Optimizer (scalar)
  #................Default: 5e-4
  #...n_epochs.....Number of epochs (integer)
  #................Default: 150
  #...n_patience...Patience in callbacks (integer)
  #................Default: 10
  #...n_batch......Size of batches is equal to n_train/n_batch (input parameter batch_size)
  #................Default: 64
  #...lay1.........Number of nodes in first and second (use half) hidden layer (integer)
  #................Default: 64 -> 32 in second
  #...actv.........Activation function of non-output layers (string)
  #................Default: "softplus"
  #...nn_verbose...Query, if network output should be shown (logical)
  #................Default: 0
  #n_cores.........Number of cores used in keras (integer)
  #................Default: NULL -> Use one less than available
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f...............Distributional forecasts (i.e. parameters) (n_test x 2 matrix)
  #......nn_ls...........Hyperparameters (list)
  #......n_train.........Number of training samples (integer)
  #......n_valid.........Number of validation samples (integer)
  #......n_test..........Number of test samples (integer)
  #......runtime_est.....Estimation time (numeric)
  #......runtime_pred....Prediction time (numeric)
  #......scores..........Data frame containing (n x 6 data frame):
  #.........pit..........PIT values of DRN forecasts (n vector)
  #.........crps.........CRPS of DRN forecasts (n vector)
  #.........logs.........Log-Score of DRN forecasts (n vector)
  #.........lgt..........Length of DRN prediction interval (n vector)
  #.........e_md.........Bias of median forecast (n vector)
  #.........e_me.........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Disable GPU
  Sys.setenv("CUDA_VISIBLE_DEVICES" = -1)
  
  # Load packages
  library(keras)
  library(tensorflow)
  
  # Number of cores
  if(is.null(n_cores)){ n_cores <- parallel::detectCores()/2 - 1 }
  
  # Set number of cores
  n_cores_config <- tf$compat$v1$ConfigProto(intra_op_parallelism_threads = as.integer(n_cores), 
                                             inter_op_parallelism_threads = as.integer(n_cores))
  tf$compat$v1$keras$backend$set_session(tf$compat$v1$Session(config = n_cores_config))
  
  #### Hyperparameter ####
  # Hyperparameters and their default values
  hpar_ls <- list(lr_adam = 5e-4, # -1 for Adam-default
                  n_epochs = 150,
                  n_patience = 10,
                  n_batch = 64,
                  lay1 = 64,
                  actv = "softplus",
                  nn_verbose = 0)
  
  # Update hyperparameters
  nn_ls <- update_hpar(hpar_ls = hpar_ls,
                       in_ls = nn_ls)
  
  # Custom optimizer
  if(nn_ls$lr_adam == -1){ custom_opt <- "adam" }
  else{ custom_opt <- optimizer_adam(lr = nn_ls$lr_adam) }
  
  # Get tensorflow probability
  tfp <- tf_probability()$distributions

  # Get standard logistic distribution
  tfd <- tfp$Normal(loc = 0, scale = 1)
  
  # Custom loss function
  custom_loss <- function(y_true, y_pred){
    # Get location and scale
    mu  <- k_dot(y_pred, k_constant(as.numeric(rbind(1, 0)), shape = c(2, 1)))
    sigma  <- k_dot(y_pred, k_constant(as.numeric(rbind(0, 1)), shape = c(2, 1)))
    
    # Standardization
    z_y <- (y_true - mu)/sigma
    
    # Calculate CRPS
    res <- sigma*( z_y*(2*tfd$cdf(z_y) - 1) + 2*tfd$prob(z_y) - 1/sqrt(pi) )
    
    # Calculate mean
    res <- k_mean(res)
    
    # Return mean CRPS
    return(res)
  }
  
  #### Data preparation ####
  # Divide data in training and validation set
  i_train <- (1:nrow(train))[-i_valid]
  
  # Read out set sizes
  n_train <- length(i_train)
  n_valid <- length(i_valid)
  n_test <- nrow(test)
  
  # Scale training data
  X_train <- scale(train[i_train,])
  
  # Save center and scale parameters
  tr_center <- attr(X_train, "scaled:center")
  tr_scale <- attr(X_train, "scaled:scale")
  
  # Scale validation data with training data attributes
  X_valid <- scale(train[i_valid,],
                   center = tr_center,
                   scale = tr_scale)
  
  # Input for fit
  X_train <- list(input = as.matrix(X_train))
  X_valid <- list(as.matrix(X_valid))
  
  # Scale data for prediction
  X_pred <- list(as.matrix(scale(test, 
                                 center = tr_center, 
                                 scale = tr_scale)))
  
  #### Build network ####
  # Input
  input <- layer_input(shape = ncol(X_train[["input"]]), name = "input")
  
  # Hidden layers
  hidden <- input %>%
    layer_dense(units = nn_ls$lay1, activation = nn_ls$actv) %>%
    layer_dense(units = nn_ls$lay1/2, activation = nn_ls$actv)
  
  # Different activation functions for output
  loc_out <- layer_dense(object = hidden, units = 1)
  scale_out <- layer_dense(object = hidden, units = 1, activation = "softplus")
  
  # Concatenate output
  output <- layer_concatenate(c(loc_out, scale_out))
  
  # Define model
  model <- keras_model(inputs = input, outputs = output)
  
  #### Estimation ####
  # Compile model
  model %>% compile(
    optimizer = custom_opt,
    loss = custom_loss
  )
  
  # Take time
  start_tm <- Sys.time()
  
  # Fit model
  history <- model %>% fit(
    x = X_train,
    y = y_train[i_train],
    epochs = nn_ls$n_epochs,
    batch_size = nn_ls$n_batch,
    validation_data = list(X_valid, y_train[i_valid]),
    verbose = nn_ls$nn_verbose,
    callbacks = callback_early_stopping(patience = nn_ls$n_patience,
                                        restore_best_weights = TRUE,
                                        monitor = "val_loss")
  )
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime_est <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  # Delete history
  rm(history)
  
  #### Prediction ####
  # Take time
  start_tm <- Sys.time()
  
  # Predict parameters of distributional forecasts (on scaled data)
  f <- predict(model, X_pred)
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime_pred <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  # Delete model
  rm(input, hidden, loc_out, scale_out, output, model)
  
  # Clear memory and session
  gc()
  k_clear_session()
  
  # Delete data
  rm(X_train, X_valid, X_pred)
  
  #### Evaluation ####
  # Calculate evaluation measures of DRN forecasts
  scores <- fn_scores_distr(f = f,
                            y = y_test,
                            distr = "norm")
  
  #### Output ####
  # Output
  return(list(f = f,
              nn_ls = nn_ls,
              scores = scores, 
              n_train = n_train,
              n_valid = n_valid,
              n_test = n_test,
              runtime_est = runtime_est,
              runtime_pred = runtime_pred))
}

#### BQN ####
## Distribution: Quantile function of Bernstein polynomials
## Estimation: Quantile loss

# Function including estimation and prediction
bqn_pp <- function(train, test, y_train, y_test, i_valid = NULL, 
                   q_levels = NULL, nn_ls = list(), n_cores = NULL){
  ###-----------------------------------------------------------------------------
  ###Input
  #train/test......Training/Test data including predictors (n_train/test x n_preds data.frame)
  #y_train/test....Training/Test observations (n_train/test vector)
  #i_valid.........Indices of validation data (n_valid vector)
  #................Default: NULL -> Use fraction of training data
  #q_levels........Quantile levels used for output and evaluation (n_q probability vector)
  #................Default: NULL -> 99 member, incl. median 
  #nn_ls...........List that may contain the following variables:
  #...p_degree.....Degree of Bernstein polynomials (integer)
  #................Default: 12
  #...n_q..........Number of equidistant quantile levels used in loss function (integer)
  #................Default: 99 (steps of 1%)
  #...lr_adam......Learning rate in Adam-Optimizer (scalar)
  #................Default: 5e-4
  #...n_epochs.....Number of epochs (integer)
  #................Default: 150
  #...n_patience...Patience for early stopping (integer)
  #................Default: 10
  #...n_batch......Size of batches is equal to n_train/n_batch (input parameter batch_size)
  #................Default: 64
  #...lay1.........Number of nodes in first hidden layer (integer)
  #................Default: 48
  #...actv.........Activation function of non-output layers (string)
  #................Default: "softplus"
  #...actv_out.....Activation function of output layer exlcuding alpha_0 (string)
  #................Default: "softplus"
  #...nn_verbose...Query, if network output should be shown (logical)
  #................Default: 0
  #n_cores.........Number of cores used in keras (integer)
  #................Default: NULL -> Use one less than available
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f...............BQN forecasts (i.e. quantiles) based on q_levels (n x n_q matrix)
  #......alpha...........BQN coefficients (n x p_degree matrix)
  #......nn_ls...........Hyperparameters (list)
  #......n_train.........Number of training samples (integer)
  #......n_valid.........Number of validation samples (integer)
  #......n_test..........Number of test samples (integer)
  #......runtime_est.....Estimation time (numeric)
  #......runtime_pred....Prediction time (numeric)
  #......scores..........Data frame containing (n x 6 data frame):
  #.........rank.........Ranks of observations in BQN forecasts (n vector)
  #.........crps.........CRPS of BQN forecasts (n vector)
  #.........logs.........Log-Score of BQN forecasts (n vector)
  #.........lgt..........Length of BQN prediction interval (n vector)
  #.........e_md.........Bias of median forecast (n vector)
  #.........e_me.........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Disable GPU
  Sys.setenv("CUDA_VISIBLE_DEVICES" = -1)
  
  # Load packages
  library(keras)
  library(tensorflow)
  
  # Number of cores
  if(is.null(n_cores)){ n_cores <- parallel::detectCores()/2 - 1 }
  
  # Set number of cores
  n_cores_config <- tf$compat$v1$ConfigProto(intra_op_parallelism_threads = as.integer(n_cores), 
                                             inter_op_parallelism_threads = as.integer(n_cores))
  tf$compat$v1$keras$backend$set_session(tf$compat$v1$Session(config = n_cores_config))
  
  # If not given use equidistant quantiles (multiple of ensemble coverage, incl. median)
  if(is.null(q_levels)){ q_levels <- seq(from = 1/100,
                                         to = 1,
                                         by = 1/100) }
  
  #### Hyperparameter ####
  # Hyperparameters and their default values
  hpar_ls <- list(p_degree = 12,
                  n_q = 99,
                  lr_adam = 5e-4, # -1 for Adam-default
                  n_epochs = 150,
                  n_patience = 10,
                  n_batch = 64,
                  lay1 = 48,
                  actv = "softplus",
                  actv_out = "softplus",
                  nn_verbose = 0)
  
  # Update hyperparameters
  nn_ls <- update_hpar(hpar_ls = hpar_ls,
                       in_ls = nn_ls)
  
  # Calculate equidistant quantile levels for loss function
  q_levels_loss <- seq(from = 1/(nn_ls$n_q + 1),
                       to = nn_ls$n_q/(nn_ls$n_q + 1),
                       by = 1/(nn_ls$n_q + 1))
  
  # Basis of Bernstein polynomials evaluated at quantile levels
  B <- sapply(0:nn_ls$p_degree, 
              dbinom, size = nn_ls$p_degree, prob = q_levels_loss)
  
  # Quantile loss functions (for neural network)
  qt_loss <- function(y_true, y_pred){
    # Quantiles calculated via basis and increments
    q  <- k_dot(k_cumsum(y_pred, axis = 0),
                k_constant(as.numeric(B), shape = c(nn_ls$p_degree + 1, nn_ls$n_q)))
    
    # Calculate individual quantile scores
    err  <- y_true - q
    e1   <- err * k_constant(q_levels_loss, shape = c(1, nn_ls$n_q))
    e2   <- err * k_constant(q_levels_loss - 1, shape = c(1, nn_ls$n_q))
    
    # Find correct values (max) and return mean
    return(k_mean( k_maximum(e1, e2), axis = 2 ))
  }
  
  # Custom optimizer
  if(nn_ls$lr_adam == -1){ custom_opt <- "adam" }
  else{ custom_opt <- optimizer_adam(lr = nn_ls$lr_adam) }
  
  #### Data preparation ####
  # Divide data in training and validation set
  i_train <- (1:nrow(train))[-i_valid]
  
  # Read out set sizes
  n_train <- length(i_train)
  n_valid <- length(i_valid)
  n_test <- nrow(test)
  
  # Scale training data of direct predictors
  X_train <- scale(train[i_train,])
  
  # Save center and scale parameters
  tr_center <- attr(X_train, "scaled:center")
  tr_scale <- attr(X_train, "scaled:scale")
  
  # Scale validation data with training data attributes
  X_valid <- scale(train[i_valid,],
                   center = tr_center,
                   scale = tr_scale)
  
  # Input for fit
  X_train <- list(input = as.matrix(X_train))
  X_valid <- list(as.matrix(X_valid))
  
  # Scale data for prediction
  X_pred <- list(as.matrix(scale(test,
                                 center = tr_center,
                                 scale = tr_scale)))
  
  #### Build network ####
  # Input
  input <- layer_input(shape = ncol(train), name = "input")
  
  # Hidden layers
  hidden <- input %>%
    layer_dense(units = nn_ls$lay1, activation = nn_ls$actv) %>%
    layer_dense(units = nn_ls$lay1/2, activation = nn_ls$actv)
  
  # Different activation functions for output (alpha_0 and positive increments)
  alpha0_out <- layer_dense(object = hidden, units = 1)
  alphai_out <- layer_dense(object = hidden, units = nn_ls$p_degree, activation = "softplus")
  
  # Concatenate output
  output <- layer_concatenate(c(alpha0_out, alphai_out))
  
  # Define model
  model <- keras_model(inputs = input, outputs = output)
  
  #### Estimation ####
  # Compile model
  model %>% compile(
    optimizer = custom_opt,
    loss = qt_loss
  )
  
  # Take time
  start_tm <- Sys.time()
  
  # Fit model
  history <- model %>% fit(
    x = X_train,
    y = y_train[i_train],
    epochs = nn_ls$n_epochs,
    batch_size = nn_ls$n_batch,
    validation_data = list(X_valid, y_train[i_valid]),
    verbose = nn_ls$nn_verbose,
    callbacks = callback_early_stopping(patience = nn_ls$n_patience,
                                        restore_best_weights = TRUE,
                                        monitor = "val_loss")
  )
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime_est <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  # Delete history
  rm(history)
  
  #### Prediction ####
  # Take time
  start_tm <- Sys.time()
  
  # Predict coefficients of Bernstein polynomials
  coeff_bern <- predict(model, X_pred)
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime_pred <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  # Delete model
  rm(input, hidden, alpha0_out, alphai_out, output, model)
  
  # Clear memory and session
  gc()
  k_clear_session()

  # Delete data
  rm(X_train, X_valid, X_pred)
  
  # Accumulate increments
  coeff_bern <- t(apply(coeff_bern, 1, cumsum))
  
  #### Evaluation ####
  # Sum up calculated quantiles (Sum of basis at quantiles times coefficients)
  q <- bern_quants(alpha = coeff_bern,
                   q_levels = q_levels)
  
  # Calculate evaluation measure of BQN forecasts
  scores <- fn_scores_ens(ens = q,
                          y = y_test,
                          skip_evals = c("e_me"),
                          scores_ens = TRUE)
  
  # Transform ranks to n_(ens + 1) bins (for multiples of (n_ens + 1) exact)
  if(ncol(q) != n_ens){ scores[["rank"]] <- ceiling(scores[["rank"]]*(n_ens + 1)/(ncol(q) + 1)) }
  
  # Calculate bias of mean forecast (formula given)
  scores[["e_me"]] <- rowMeans(coeff_bern) - y_test
  
  #### Output ####
  return(list(f = q, 
              alpha = coeff_bern,
              nn_ls = nn_ls,
              scores = scores, 
              n_train = n_train,
              n_valid = n_valid,
              n_test = n_test,
              runtime_est = runtime_est,
              runtime_pred = runtime_pred))
}

#### HEN ####
## Distribution: Piecewise uniform distribution
## Estimation: MLE

# Function for pp including estimation and prediction
hen_pp <- function(train, test, y_train, y_test, i_valid = NULL, 
                   nn_ls = list(), n_cores = NULL){
  ###-----------------------------------------------------------------------------
  ###Input
  #train/test......Training/Test data including predictors (n_train/test x n_preds data.frame)
  #y_train/test....Training/Test observations (n_train/test vector)
  #i_valid.........Indices of validation data (n_valid vector)
  #................Default: NULL -> Use fraction of training data
  #nn_ls...........List that may contain the following variables:
  #...bin_edges....Boundaries of bins starting from 0 ((n_bins + 1) vector)
  #................Default: NULL -> Use default partition
  #...lr_adam......Learning rate in Adam-Optimizer (scalar)
  #................Default: 5e-4
  #...n_epochs.....Number of epochs (integer)
  #................Default: 150
  #...n_patience...Patience in callbacks (integer)
  #................Default: 10
  #...n_batch......Size of batches is equal to n_train/n_batch (input parameter batch_size)
  #................Default: 64
  #...lay1.........Number of nodes in first hidden layer (integer)
  #................Default: 64
  #...actv.........Activation function of non-output layers (string)
  #................Default: "relu"
  #...nn_verbose...Query, if network output should be shown (logical)
  #................Default: 0
  #n_cores.........Number of cores used in keras (integer)
  #................Default: NULL -> Use one less than available
  ###-----------------------------------------------------------------------------
  ###Output
  #res...List containing:
  #......f...............HEN forecasts (n x n_bins matrix)
  #......nn_ls...........Hyperparameters (list)
  #......n_train.........Number of training samples (integer)
  #......n_valid.........Number of validation samples (integer)
  #......n_test..........Number of test samples (integer)
  #......runtime_est.....Estimation time (numeric)
  #......runtime_pred....Prediction time (numeric)
  #......scores_ens/pp...Data frame containing (n x 6 data frame):
  #.........rank/pit.....Ranks of ensemble / PIT values of HEN forecasts (n vector)
  #.........crps.........CRPS of ensemble/HEN forecasts (n vector)
  #.........logs.........Log-Score of ensemble/HEN forecasts (n vector)
  #.........lgt..........Ensemble range / Length of HEN prediction interval (n vector)
  #.........e_md.........Bias of median forecast (n vector)
  #.........e_me.........Bias of mean forecast (n vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Disable GPU
  Sys.setenv("CUDA_VISIBLE_DEVICES" = -1)
  
  # Load packages
  library(keras)
  library(tensorflow)
  
  # Number of cores
  if(is.null(n_cores)){ n_cores <- parallel::detectCores()/2 - 1 }
  
  # Set number of cores
  n_cores_config <- tf$compat$v1$ConfigProto(intra_op_parallelism_threads = as.integer(n_cores), 
                                             inter_op_parallelism_threads = as.integer(n_cores))
  tf$compat$v1$keras$backend$set_session(tf$compat$v1$Session(config = n_cores_config))
  
  # Function to create bin edges
  get_edges <- function(obs, N){
    # Quantile-based binning
    bin_edges <- unique(round(quantile(x = obs,
                                       probs = (0:(N - 1))/(N - 1)), 
                              digits = 2))
    
    # Output
    return(bin_edges)
  }
  
  #### Hyperparameter ####
  # Hyperparameters and their default values
  hpar_ls <- list(bin_edges = get_edges(obs = y_train,
                                        N = 50),
                  lr_adam = 5e-4, # -1 for Adam-default
                  n_epochs = 150,
                  n_patience = 10,
                  n_batch = 64,
                  lay1 = 64,
                  actv = "softplus",
                  nn_verbose = 0)
  
  # Update hyperparameters
  nn_ls <- update_hpar(hpar_ls = hpar_ls,
                       in_ls = nn_ls)

  # Custom optimizer
  if(nn_ls$lr_adam == -1){ custom_opt <- "adam" }
  else{ custom_opt <- optimizer_adam(lr = nn_ls$lr_adam) }
  
  # Choice of loss function
  custom_loss <- 'categorical_crossentropy'
  
  # Number of bins
  n_bins <- length(nn_ls$bin_edges) - 1
  
  #### Data preparation ####
  # Divide data in training and validation set
  i_train <- (1:nrow(train))[-i_valid]
  
  # Read out set sizes
  n_train <- length(i_train)
  n_valid <- length(i_valid)
  n_test <- nrow(test)
  
  # Calculate bin of each observation (requires lowest edge to be leq than smallest obs.)
  # (Note: Python starts indexing at 0!)
  y_train <- sapply(y_train, function(z) sum(z >= nn_ls$bin_edges[-length(nn_ls$bin_edges)])) - 1
  
  # Generate categorical matrices (via Keras function)
  y_train <- to_categorical(y_train, n_bins)
  
  # Scale training data of direct predictors
  X_train <- scale(train[i_train,])
  
  # Save center and scale parameters
  tr_center <- attr(X_train, "scaled:center")
  tr_scale <- attr(X_train, "scaled:scale")
  
  # Scale validation data with training data attributes
  X_valid <- scale(train[i_valid,],
                   center = tr_center,
                   scale = tr_scale)
  
  # Input for fit
  X_train <- list(input = as.matrix(X_train))
  X_valid <- list(as.matrix(X_valid))
  
  # Scale data for prediction
  X_pred <- list(as.matrix(scale(test,
                                 center = tr_center,
                                 scale = tr_scale)))
  
  #### Build network ####
  # Input
  input <- layer_input(shape = ncol(train), name = "input")
  
  # Hidden layers
  hidden <- input %>%
    layer_dense(units = nn_ls$lay1, activation = nn_ls$actv) %>%
    layer_dense(units = nn_ls$lay1/2, activation = nn_ls$actv)
  
  # Output
  output <- hidden %>%
    layer_dense(units = n_bins, activation = 'softmax')
  
  # Define model
  model <- keras_model(inputs = input, outputs = output)
  
  #### Estimation ####
  # Compile model
  model %>% compile(
    loss = custom_loss,
    optimizer = custom_opt,
    metrics = c('accuracy')
  )
  
  # Take time
  start_tm <- Sys.time()
  
  # Fit model
  history <- model %>% fit(
    x = X_train,
    y = y_train[i_train,],
    epochs = nn_ls$n_epochs,
    batch_size = nn_ls$n_batch,
    validation_data = list(X_valid, y_train[i_valid,]),
    verbose = nn_ls$nn_verbose,
    callbacks = callback_early_stopping(patience = nn_ls$n_patience,
                                        restore_best_weights = TRUE,
                                        monitor = "val_loss")
  )
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime_est <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  # Delete history
  rm(history)
  
  #### Prediction ####
  # Take time
  start_tm <- Sys.time()
  
  # Predict bin probabilities
  p_bins <- predict(model, X_pred)
  
  # Take time
  end_tm <- Sys.time()
  
  # Time needed
  runtime_pred <- as.numeric(difftime(end_tm, start_tm, units = "mins"))
  
  # Delete model
  rm(input, hidden, output, model)
  
  # Clear memory and session
  gc()
  k_clear_session()
  
  # Delete data
  rm(X_train, X_valid, X_pred)
  
  #### Evaluation ####
  # Calculate scores
  scores <- fn_scores_hd(f = p_bins,
                         y = y_test,
                         bin_edges = nn_ls[["bin_edges"]])
  
  #### Output ####
  return(list(f = p_bins,
              nn_ls = nn_ls,
              scores = scores,
              n_train = n_train,
              n_valid = n_valid,
              n_test = n_test,
              runtime_est = runtime_est,
              runtime_pred = runtime_pred))
}
