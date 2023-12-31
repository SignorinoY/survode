% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.survode.R
\name{predict.survode}
\alias{predict.survode}
\title{Predict the baseline hazard function}
\usage{
\method{predict}{survode}(object, type = "hazard", time, ...)
}
\arguments{
\item{object}{Object of class \code{survode} containing the objectted model.}

\item{type}{Type of prediction. Currently only \code{"hazard"} is supported.}

\item{time}{Vector of times at which to evaluate the baseline hazard.}

\item{...}{Additional arguments passed to \code{\link{predict}}.}
}
\value{
A list containing the following components:
\item{time}{Times at which the baseline hazard was evaluated.}
\item{basehaz}{Baseline hazard estimates.}
\item{cumhaz}{Cumulative baseline hazard estimates.}
}
\description{
Predict the baseline hazard function
}
