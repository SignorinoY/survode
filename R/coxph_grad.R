coxph_grad <- function(parms, time, status, x, parms_bs) {
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
  dcbh <- deSolve::ode(
    y = rep(0, q),
    times = c(0, time),
    func = coxph_ode_grad,
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
  eta <- as.vector(x %*% beta)
  dalpha <- colSums(dcbh * exp(eta) - b * status)
  dbeta <- colSums(x * (cbh * exp(eta) - status))
  c(dalpha, dbeta)
}
