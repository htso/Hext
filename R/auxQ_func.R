#' Compute EM Auxiliary Q function for beta distribution
#'
#' @param theta the shape parameters of beta distribution
#' @param x numeric vector to evaluate Aux Q function
#' @param wt weight vector
#' @description See ref
#' @return numeric value of Aux Q function
#' @export
Qaux.beta = function(theta, x, wt) {
  s1 = theta[1] # shape1
  s2 = theta[2] # shape2
  lb = lgamma(s1+s2) - lgamma(s1) - lgamma(s2) + (s1-1)*log(x) + (s2-1)*log(1-x)
  return(-sum(wt*lb))
}

#' Compute EM Auxiliary Q function for gamma distribution
#'
#' @param theta the shape parameters of gamma distribution
#' @param x numeric vector to evaluate Aux Q function
#' @param wt weight vector
#' @description See ref
#' @return numeric value of Aux Q function
#' @export
Qaux.gamma = function(theta, x, wt) {
  a = theta[1] # shape
  s = theta[2] # scale
  lg = (a-1)*log(x) - x/s - lgamma(a) - a*log(s)
  return(-sum(wt*lg))
}


