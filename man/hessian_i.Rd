% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatComp22077.R
\name{hessian_i}
\alias{hessian_i}
\title{Caluclate the numeric hessian of the loglikelihood}
\usage{
hessian_i(data_i, symbolic_hessian, param_vals, dimension)
}
\arguments{
\item{data_i}{the parameter vector}

\item{symbolic_hessian}{penalty coefficient}

\item{param_vals}{the predictive variable matrix}

\item{dimension}{the response variable}
}
\value{
The loglikelihood function of beta
}
\description{
Caluclate the numeric hessian of the loglikelihood
}
