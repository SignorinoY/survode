coxph_ode_hess <- function(t, y, parms) {
  b <- splines::bs(
    t,
    knots = parms$knots,
    Boundary.knots = parms$Boundary.knots,
    degree = parms$degree,
    intercept = TRUE
  )
  list(t(b) %*% b * exp(sum(parms$alpha * b)))
}
