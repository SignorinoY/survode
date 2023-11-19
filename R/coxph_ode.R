coxph_ode <- function(t, y, parms) {
  b <- splines::bs(
    t,
    knots = parms$knots,
    Boundary.knots = parms$Boundary.knots,
    degree = parms$degree,
    intercept = TRUE
  )
  list(exp(sum(parms$alpha * b)))
}
