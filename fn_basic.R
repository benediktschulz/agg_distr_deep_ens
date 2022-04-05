## Function file
# Helper functions

#### Update hyperparameters ####
update_hpar <- function(hpar_ls, in_ls){
  ###-----------------------------------------------------------------------------
  ###Input
  #hpar_ls...Default hyperparameter (list)
  #in_ls.....Selected hyperparameter given by user (list)
  ###-----------------------------------------------------------------------------
  ###Output
  #hpar_ls...All hyperparameters including users selection (list)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Names of hyperparameters
  hpar_names <- names(hpar_ls)
  
  # Names of hyperparameter to update
  in_names <- names(in_ls)
  
  #### Update ####
  # Loop over given names
  for(temp_hpar in in_names){
    # Update list if correct name is given
    if(is.element(temp_hpar, hpar_names)){ hpar_ls[[temp_hpar]] <- in_ls[[temp_hpar]] }
    else{ print(paste0("Wrong hyperparameter given: ", temp_hpar))}
  }
  
  #### Output ####
  # Return list
  return(hpar_ls)
}

#### Remove constant columns ####
# Function that removes constant columns of data-frame
rm_const <- function(data, cols = NULL, t_c = 0){
  ###-----------------------------------------------------------------------------
  ###Input
  #data...Data to check (data frame)
  #cols...Columns to check (String or integer vector)
  #.......Default: NULL -> Check all
  #t_c....Threshold for (almost) constant column (non-negative scalar)
  #.......Default: 0 -> Constant
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Columns of data that are not constant (String or integer vector)
  ###-----------------------------------------------------------------------------
  
  #### Initiation ####
  # Set cols if not given
  if(is.null(cols)){ cols <- colnames(data) }
  
  # Use only data that is needed
  data <- data[,cols]
  
  #### Remove columns ####
  # Number of samples to check with
  n_check <- min(10, nrow(data))
  
  # Check on sample which rows are candidates (-> computational more feasible)
  # bool_res <- (apply(data[sample(1:nrow(data), n_check),], 2, function(x) sd(as.numeric(x)) ) <= t_c)
  bool_res <- (apply(data[sample(1:nrow(data), n_check),], 2, sd) <= t_c)
  
  # Check if any of candidates is constant (Special case: Apply only on matrices)
  if(sum(bool_res) == 1){ bool_res[bool_res] <- (sd(data[,bool_res]) <= t_c) }
  else if(any(bool_res)){ bool_res[bool_res] <- (apply(data[,bool_res], 2, sd) <= t_c) }
  
  #### Output ####
  # Return columns that are not (almost) constant
  return(cols[!bool_res])
}