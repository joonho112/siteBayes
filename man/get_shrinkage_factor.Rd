% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_shrinkage_factor.R
\name{get_shrinkage_factor}
\alias{get_shrinkage_factor}
\title{Compute a "shrinkage factor" based on average site size, cross-site variation, and sampling assumptions}
\usage{
get_shrinkage_factor(nj_mean = 10, sigma_tau = 0.15, varY = 1, p = 0.5, R2 = 0)
}
\arguments{
\item{nj_mean}{Positive numeric (default 10). Interpreted as the \emph{average site size}
across all sites, used to approximate a typical within-site sampling variance.}

\item{sigma_tau}{Positive numeric (default 0.15). The cross-site standard deviation
of the true effects.}

\item{varY}{Positive numeric (default 1). The assumed variance of the outcome \eqn{Y}.
Often set to 1 if working in effect-size units.}

\item{p}{Numeric in (0,1) (default 0.5). Proportion of units treated in a typical site.}

\item{R2}{Numeric in [0,1) (default 0.0). Proportion of variance in \eqn{Y} explained by
covariates (level-1).}
}
\value{
A numeric value (between 0 and 1) representing the shrinkage factor.
}
\description{
This function calculates a typical shrinkage factor of the form:
\deqn{\frac{\sigma_{\tau}^2}{\sigma_{\tau}^2 + GM(\widehat{se}^2_j)}.}
under the assumption that:
\itemize{
\item The \emph{average} site size is \eqn{\bar{n}_j},
\item The cross-site standard deviation of the true effects is \eqn{\sigma_{\tau}},
\item The outcome variance is \eqn{\mathrm{var}(Y)} (often 1 if working in effect-size units),
\item A fixed proportion \eqn{p} of units are treated,
\item A proportion \eqn{R^2} of outcome variance is explained by covariates.
}
}
\details{
Specifically, we define:
\deqn{\kappa = \mathrm{var}(Y) \times
      \left(\frac{1}{p} + \frac{1}{1 - p}\right) \times (1 - R^2),}
and compute:
\deqn{GM(\widehat{se}^2_j) = \kappa \times \frac{1}{\bar{n}_j}.}
The function then returns:
\deqn{\frac{\sigma_{\tau}^2}{\sigma_{\tau}^2 + GM(\widehat{se}^2_j)}.}
}
\examples{
# Example 1: Average site size is ~10, cross-site SD=0.15, half units treated,
# outcome variance=1, R2=0
get_shrinkage_factor(nj_mean=10, sigma_tau=0.15, varY=1, p=0.5, R2=0)

# Example 2: If the average site size is larger (e.g., 50),
# then we have a different typical within-site variance => a different shrinkage factor
get_shrinkage_factor(nj_mean=50, sigma_tau=0.15, varY=1, p=0.5, R2=0)

}
