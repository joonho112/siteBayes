% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_priorG2.R
\name{gen_priorG2}
\alias{gen_priorG2}
\title{Generate Site-Specific Treatment Effects from a Prior Distribution}
\usage{
gen_priorG2(
  true_dist,
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
  seed = 123,
  set_seed = TRUE
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

\item{seed}{An optional integer for setting the random seed (default: 123). If \code{NULL}, seed is not set.}

\item{set_seed}{Logical indicating whether to set the seed (default: TRUE).}
}
\value{
A numeric vector of length \code{J} containing the generated treatment effects.
}
\description{
This function generates site-specific treatment effects based on the specified
distribution and parameters. The generated values are scaled by \code{sigma_tau}.
}
\details{
You can specify site-level covariates in two ways:
\enumerate{
\item Provide a \code{formula} and a \code{data} frame, which will be processed via
\code{model.matrix}. The resulting design matrix is multiplied by \code{beta}
to shift each site's mean.
\item Or simply set \code{tau} if you want a constant baseline mean with no additional covariates.
}

The final site-specific means are \code{site_means = tau + model_matrix \%*\% beta} (if a formula is provided),
or simply \code{rep(tau, J)} if no formula and beta are given. Then, for each site, a value is drawn
according to the chosen distribution with that mean and variance \code{variance}.
Finally, all values are multiplied by \code{sigma_tau}.
}
\examples{
\dontrun{
  # 1) No covariates (Gaussian distribution, just a constant mean = tau)
  effects1 <- gen_priorG2("Gaussian", J = 5, tau = 2, variance = 1)
  effects1

  # 2) Using formula & data (Gaussian distribution)
  #    Suppose we have 2 covariates (X1, X2) for 5 sites
  site_data <- data.frame(
    X1 = c(1, 2, 3, 4, 5),
    X2 = c(0, 1, 0, 1, 0)
  )
  # Model formula includes an intercept by default (~ X1 + X2)
  # => design matrix has columns: (Intercept), X1, X2
  # => length(beta) must be 3
  betas <- c(10, 1, -2)
  effects2 <- gen_priorG2("Gaussian", formula = ~ X1 + X2, beta = betas, data = site_data)
  effects2

  # 3) T distribution with covariates
  effects3 <- gen_priorG2("T", formula = ~ 0 + X1 + X2, beta = c(1, -2), data = site_data,
                          nu = 5, variance = 2)
  effects3

  # 4) Asymmetric Laplace (requires LaplacesDemon)
  #    We'll also use tau offset
  if (requireNamespace("LaplacesDemon", quietly = TRUE)) {
    effects4 <- gen_priorG2("ALD", formula = ~ X1 + X2, beta = betas, data = site_data,
                            rho = 2, tau = 5)
    effects4
  }
}

}
