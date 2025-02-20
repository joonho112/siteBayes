% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_multisite_data.R
\name{sim_multisite_data}
\alias{sim_multisite_data}
\title{Simulate a complete multi-site dataset (no final n_j), with optional rank correlation}
\usage{
sim_multisite_data(
  true_dist = "Gaussian",
  J = 50,
  sigma_tau = 0.25,
  tau = 0,
  variance = 1,
  nu = NULL,
  slant = NULL,
  rho = NULL,
  delta = NULL,
  eps = NULL,
  ups = NULL,
  formula = NULL,
  beta = NULL,
  data = NULL,
  nj_mean = 10,
  cv = 0.3,
  nj_min = 5,
  p = 0.5,
  R2 = 0,
  precision_dependence = FALSE,
  rank_corr = 0.5,
  max_iter = 20000,
  tol = 0.01,
  set_seed = FALSE
)
}
\arguments{
\item{true_dist}{A character string specifying the type of distribution to use.
Options are "Gaussian", "T", "Skew", "ALD", or "Mixture".}

\item{J}{A positive integer indicating the number of sites (default: 50). If \code{formula}
and \code{data} are provided, \code{J} should match \code{nrow(data)} or can be
left as default for consistency check.}

\item{sigma_tau}{A positive number to scale the generated treatment effects (default: 0.25).}

\item{tau}{A numeric value specifying an overall offset for the site means (default: 0).}

\item{variance}{A positive number specifying the baseline variance for the distribution (default: 1).}

\item{nu}{Degrees of freedom for the T distribution. Must be provided if \code{true_dist} is "T".}

\item{slant}{A numeric value for the skew parameter in the Skewed Normal distribution.
Must be provided if \code{true_dist} is "Skew".}

\item{rho}{A numeric value for the parameter in the Asymmetric Laplace distribution.
Must be provided if \code{true_dist} is "ALD".}

\item{delta}{A numeric value for the Mixture distribution.}

\item{eps}{A numeric value (typically between 0 and 1) for the Mixture distribution.}

\item{ups}{A positive numeric value for the Mixture distribution.}

\item{formula}{An optional model formula describing how covariates in \code{data} should be used
to form the design matrix. For example, \code{~ X1 + X2} or \code{~ 0 + X1 + X2}.
If provided, you must also supply \code{beta} and \code{data}.}

\item{beta}{A numeric vector of regression coefficients corresponding to the columns of
\code{model.matrix(formula, data)}. If \code{formula} is provided, \code{beta} must
have length equal to \code{ncol(model.matrix(...))}.}

\item{data}{An optional data frame containing site-level covariates to be used in \code{formula}.
Must have \code{nrow(data) == J} if provided.}

\item{nj_mean}{Positive numeric. Desired truncated-gamma mean of site sizes.}

\item{cv}{Numeric >= 0. Coefficient of variation for the truncated Gamma.}

\item{nj_min}{Positive numeric. Lower bound for truncation (default 5).}

\item{p}{Numeric in (0,1). Proportion assigned to treatment (defaults to 0.5).}

\item{R2}{Numeric in [0,1). Proportion of variance explained by level-1 covariates
(defaults to 0.0).}

\item{precision_dependence}{Logical indicating whether to enforce a partial correlation
(default \code{FALSE}). If \code{FALSE}, we random-shuffle \code{se2_j} to break
any existing correlation (near-zero result).
If \code{TRUE}, we reorder \code{se2_j} to achieve \code{rank_corr}.}

\item{rank_corr}{Numeric in \eqn{[-1,1]}, the target Spearman rank correlation to achieve
between \eqn{\tau_j} and \eqn{\text{se2}_j} (default 0.50). Only used if
\code{precision_dependence = TRUE}.}

\item{max_iter}{Integer. Maximum number of swap attempts for the iterative
correlation injection (default 20000).}

\item{tol}{Numeric tolerance for stopping early if the current correlation is
within \code{tol} of \code{rank_corr} (default 0.01).}

\item{set_seed}{Logical (default \code{FALSE}). If \code{TRUE}, we do \code{set.seed(123)}
at the start. If \code{FALSE}, no seed is set.}
}
\value{
A tibble with columns:
\describe{
\item{\code{tau_j}}{The true site effects.}
\item{\code{se2_j}}{The final within-site variances.}
\item{\code{tau_j_hat}}{Observed site effects \(\sim N(\tau_j, se2_j)\).}
\item{\code{corr_est}}{Spearman rank correlation between \code{tau_j} and \code{se2_j}.}
}
}
\description{
This function does the following:
\enumerate{
\item Calls \code{\link{gen_priorG2}} to generate \eqn{\tau_j} (true effects).
\item Calls \code{\link{sim_sitesize_withinvar}} to generate \eqn{\text{se2}_j}.
(Internally also generates \eqn{n_j}, but that is not included in the final output.)
\item Calls \code{\link{sim_observed_effects}} to optionally reorder \eqn{\text{se2}_j}
for partial Spearman correlation \eqn{\rho=\text{rank_corr}}, then simulates
\(\hat{\tau}_j \sim N(\tau_j, se2_j)\).
}
}
\details{
The user can specify \code{set_seed=TRUE} to do \code{set.seed(123)} for reproducible results.
By default, \code{set_seed=FALSE} means no seed is set. The final returned tibble has columns
\code{tau_j, se2_j, tau_j_hat, corr_est} in that order, omitting any site size column.
}
\examples{
\dontrun{
# Example:
df_sim <- sim_multisite_data(
  # For gen_priorG2:
  true_dist = "Gaussian", J = 10, variance = 1,
  # For sim_sitesize_withinvar:
  nj_mean = 20, cv = 0.5,
  # partial correlation ~ 0.4
  precision_dependence = TRUE, rank_corr = 0.4,
  set_seed=TRUE
)
df_sim

cor(df_sim$tau_j, df_sim$se2_j, method="spearman") # near 0.4
}

}
