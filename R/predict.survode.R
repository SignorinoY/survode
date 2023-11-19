#' Predict the baseline hazard function
#'
#' @param object Object of class \code{survode} containing the objectted model.
#' @param type Type of prediction. Currently only \code{"hazard"} is supported.
#' @param time Vector of times at which to evaluate the baseline hazard.
#' @param ... Additional arguments passed to \code{\link{predict}}.
#' @return A list containing the following components:
#' \item{time}{Times at which the baseline hazard was evaluated.}
#' \item{basehaz}{Baseline hazard estimates.}
#' \item{cumhaz}{Cumulative baseline hazard estimates.}
#' @export
predict.survode <- function(object, type = "hazard", time, ...) {
  if (missing(object)) stop("a object argument is required")
  if (!inherits(object, "survode")) stop("object must be of class survode")

  type <- match.arg(type, c("hazard"))
  if (type != "hazard") stop("only type = 'hazard' is supported")

  if (missing(time)) stop("a time argument is required")
  if (!is.numeric(time)) stop("time must be numeric")
  if (any(time < 0)) stop("time must be non-negative")

  stime <- sort(time)
  otime <- order(time)
  b <- splines::bs(
    time,
    knots = object$coefficients$gamma$knots,
    Boundary.knots = object$coefficients$gamma$Boundary.knots,
    degree = object$coefficients$gamma$degree,
    intercept = TRUE
  )
  bh <- exp(b %*% object$coefficients$gamma$alpha)
  parms_ode <- object$coefficients$gamma
  cbh <- deSolve::ode(
    y = 0,
    times = c(0, stime),
    func = coxph_ode,
    parms = parms_ode,
    method = "ode45"
  )[-1, ][order(otime), ][, -1]

  list(time = time, basehaz = bh, cumhaz = cbh)
}
