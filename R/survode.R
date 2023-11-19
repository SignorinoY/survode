#' Fit Survial Models via Ordinary Differential Equations
#'
#' @param formula a formula object, with the response on the left of a \code{~}
#'   operator, and the terms on the right. The response must be a survival
#'   object as returned by the \code{\link{Surv}} function.
#' @param data a data frame in which to interpret the variables named in the
#'   \code{formula}.
#' @param df the number of degrees of freedom for the spline.
#' @param degree the degree of the spline. Default is 3.
#' @param knots the knots of the spline, either \code{"uniform"} or
#'   \code{"quantile"}. For \code{"uniform"}, the knots are equally spaced
#'    between 0 and the maximum observed time. For \code{"quantile"}, the knots
#'    are equally spaced between the quantiles of the observed times. Default is
#'     \code{"uniform"}.
#' @param init a vector of initial values for the parameters. Default initial
#'   value is zero for all parameters.
#' @param control Object of class \link{survode_control} containing
#'   control parameters for the fitting algorithm. Default is
#'  \code{survode_control(...)}.
#' @param ... Other arguments passed to \code{\link{survode_control}}.
#' @importFrom survival Surv
#' @export
#' @examples
#' library(simsurv)
#' library(survode)
#'
#' # Create a simple data set
#' set.seed(42)
#' cov <- data.frame(
#'   id = 1:400,
#'   trt = stats::rbinom(400, 1L, 0.5),
#'   hormon = stats::rnorm(400, 0, 1)
#' )
#' dat <- simsurv(
#'   lambdas = 0.1, gammas = 1.5, betas = c(trt = -0.5, hormon = 1.0),
#'   x = cov, maxt = 5
#' )
#' sim <- merge(cov, dat)
#'
#' # Fit a Cox model with a spline for the baseline hazard
#' fit <- survode(Surv(eventtime, status) ~ trt + hormon, data = sim, df = 5)
#' fit$coefficients$beta
#'
#' # Predict the baseline hazard at times (0, 5)
#' bh <- predict(fit, type = "hazard", time = seq(0, 5, length.out = 100))
#' par(mfrow = c(1, 2), cex = 0.85)
#' plot(
#'   bh$time, bh$basehaz,
#'   type = "l", xlab = "Time", ylab = "Baseline hazard"
#' )
#' lines(bh$time, 0.15 * bh$time^0.5, col = "red")
#' plot(
#'   bh$time, bh$cumhaz,
#'   type = "l", xlab = "Time", ylab = "Cumulative baseline hazard"
#' )
#' lines(bh$time, 0.1 * bh$time^1.5, col = "red")
survode <- function(
    formula, data, df, degree = 3, knots = "uniform", init, control, ...) {
  call <- match.call()
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(init)) init <- NULL
  if (missing(control)) control <- survode_control(...)
  if (degree <= 1) stop("degree must be greater than 1")
  degree <- as.integer(degree)
  if (df < (degree + 1)) stop("df must be greater than degree + 1")
  df <- as.integer(df)
  knots <- match.arg(knots, c("uniform", "quantile"))
  nknots <- df - degree - 1

  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)

  if (knots == "quantile") {
    probs <- seq(0, 1, length.out = nknots + 2)
    knots <- stats::quantile(y[, 1], probs = probs)
    knots <- knots[-c(1, nknots + 2)]
  } else {
    knots <- seq(0, max(y[, 1]), length.out = nknots + 2)
    knots <- knots[-c(1, nknots + 2)]
  }
  boundary_knots <- c(0, max(y[, 1]))
  parms_bs <- list(
    knots = knots,
    Boundary.knots = boundary_knots,
    degree = degree
  )

  x <- stats::model.matrix(formula, data)
  x <- x[, !colnames(x) %in% "(Intercept)"]
  if (!is.matrix(x)) x <- as.matrix(x)
  time <- y[, 1]
  status <- y[, 2]

  n <- nrow(x)
  nvar <- ncol(x)
  nevent <- sum(status)
  if (!missing(init) && length(init) > 0) {
    if (length(init) != nvar + df) {
      stop("The length of init must be equal to the number of variables.")
    }
  } else {
    init <- rep(0, nvar + df)
  }

  sorted <- order(time)
  xsorted <- x[sorted, ]
  if (!is.matrix(xsorted)) xsorted <- as.matrix(xsorted)
  survfit <- optimx::optimr(
    init,
    fn = coxph_loss,
    gr = coxph_grad,
    hess = coxph_hess,
    time = time[sorted],
    status = status[sorted],
    x = xsorted,
    parms_bs = parms_bs,
    method = "snewton",
    control = list(
      maxit = control$maxit,
      eps = control$eps
    )
  )
  coefficients <- list(
    gamma = list(
      alpha = survfit$par[1:df],
      knots = parms_bs$knots,
      Boundary.knots = parms_bs$Boundary.knots,
      degree = parms_bs$degree
    ),
    beta = survfit$par[(df + 1):(df + nvar)]
  )
  info <- coxph_hess(survfit$par, time, status, x, parms_bs, lambda = 0)
  var <- MASS::ginv(info / n)

  rval <- list(
    coefficients = coefficients,
    var = var,
    loglik = survfit$value,
    iter = survfit$counts,
    n = n,
    nevent = nevent,
    first = survfit$grad,
    formula = formula,
    call = call
  )
  class(rval) <- "survode"

  return(rval)
}
