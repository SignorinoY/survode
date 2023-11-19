coxph_hess <- function(parms, time, status, x, parms_bs, lambda = 0.1) {
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
  d2cbh <- deSolve::ode(
    y = rep(0, q * q),
    times = c(0, time),
    func = coxph_ode_hess,
    parms = parms_ode,
    method = "ode45"
  )[-1, ][, -1]
  eta <- as.vector(x %*% beta)
  d2alpha <- matrix(colSums(sweep(d2cbh, 1, exp(eta), "*")), nrow = q)
  d2beta <- t(x) %*% diag(exp(eta) * cbh) %*% x
  dalph_dabeta <- t(dcbh) %*% diag(exp(eta)) %*% x

  hess <- cbind(rbind(d2alpha, t(dalph_dabeta)), rbind(dalph_dabeta, d2beta))
  hess <- as.matrix(hess)
  if (det(hess) < 1e-10) {
    hess <- hess + lambda * diag(q + p)
  }
  return(hess)
}
