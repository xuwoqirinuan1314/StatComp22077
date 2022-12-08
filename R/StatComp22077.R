#' @title Caluclate the numeric hessian of the loglikelihood
#' @description Caluclate the numeric hessian of the loglikelihood
#' @param data_i the parameter vector
#' @param symbolic_hessian penalty coefficient
#' @param param_vals the predictive variable matrix
#' @param dimension the response variable
#' @return The loglikelihood function of beta
#' @export
hessian_i = function(
  data_i, #i^th observation: (Y_i, X_i, Z_i)
  symbolic_hessian, #matrix with elements that are the symbolic expression for the derivatives
  param_vals, #estimates of first stage and second stage parameter estimates (theta and beta)
  dimension #the dimension of the hessian matrix
  ){
  
  env = as.list(c(data_i, param_vals))

  # Iterate through every element in symbolic_hessian and evaluate it
  # at the i^th observation using the given parameter estimates
  hessian = matrix(
    sapply(symbolic_hessian, eval, env = env, enclos = NULL), 
    dimension)
  
  return(hessian)
  
}

