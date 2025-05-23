% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prior_rubin_normal.R
\name{student_t_prior}
\alias{student_t_prior}
\title{Create a Student-t prior object}
\usage{
student_t_prior(df, loc, scale)
}
\arguments{
\item{df}{degrees of freedom (>0)}

\item{loc}{location}

\item{scale}{scale (>0)}
}
\value{
A list with dist="student_t" and values=c(df, loc, scale)
}
\description{
Create a Student-t prior object
}
