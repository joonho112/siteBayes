% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_rubin_normal.R
\name{interpret_prior}
\alias{interpret_prior}
\title{Interpret a user-facing prior object into an integer-coded family}
\usage{
interpret_prior(prior_obj)
}
\arguments{
\item{prior_obj}{a list from e.g. normal_prior(), cauchy_prior(), etc.}
}
\value{
a list with \code{family_id} and \code{param_vector}
}
\description{
This function is used internally to transform e.g. \code{list(dist="normal", values=c(0,10))}
into e.g. \code{family_id=1, param_vector=c(0,10,0)} for the Stan model to parse.
}
\keyword{internal}
