% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.par.start.R
\name{get.par.start}
\alias{get.par.start}
\title{Get starting values for parameters}
\usage{
get.par.start(formula, data)
}
\arguments{
\item{formula}{the estimation formula}

\item{data}{data matrix with y, u, l, and covariates, including cohort}
}
\value{
Named vector of initial parameter estimates
}
\description{
Uses linear modeling to compute initial values for MLE optimizer
}
