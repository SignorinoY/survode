#' Ancillary arguments for controlling survode fitting
#' @param eps the convergence criteria for the optimization algorithm.
#'   Default is 1e-4.
#' @param maxit the maximum number of iterations for the optimization algorithm.
#'   Default is 100.
survode_control <- function(eps = 1e-4, maxit = 100) {
  if (eps <= 0) stop("Invalid convergence criteria")
  if (maxit < 0) stop("Invalid value for iterations")
  list(eps = eps, maxit = as.integer(maxit))
}
