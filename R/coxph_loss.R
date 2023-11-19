coxph_loss <- function(parms, time, status, x, parms_bs) {
  p <- ncol(x)
  q <- length(parms_bs$knots) + parms_bs$degree + 1

  alpha <- parms[1:q]
  beta <- parms[(q + 1):(q + p)]

  parms_ode <- parms_bs
  parms_ode$alpha <- alpha
  cbh <- deSolve::ode(
    y = 0,
    times = c(0, time),
    func = coxph_ode,
    parms = parms_ode,
    method = "ode45"
  )[-1, ][, -1]

  b <- splines::bs(
    time,
    knots = parms_bs$knots,
    Boundary.knots = parms_bs$Boundary.knots,
    degree = parms_bs$degree,
    intercept = TRUE
  )
  eta <- x %*% beta

  sum(cbh * exp(eta) - status * (b %*% alpha + eta))
}
