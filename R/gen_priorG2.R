#' Generate Site-Specific Treatment Effects from a Prior Distribution
#'
#' This function generates site-specific treatment effects based on the specified
#' distribution and parameters. The generated values are scaled by \code{sigma_tau}.
#'
#' You can specify site-level covariates in two ways:
#' 1) Provide a \code{formula} and a \code{data} frame, which will be processed via
#'    \code{model.matrix}. The resulting design matrix is multiplied by \code{beta}
#'    to shift each site's mean.
#' 2) Or simply set \code{tau} if you want a constant baseline mean with no additional covariates.
#'
#' The final site-specific means are \code{site_means = tau + model_matrix \%*\% beta} (if a formula is provided),
#' or simply \code{rep(tau, J)} if no formula and beta are given. Then, for each site, a value is drawn
#' according to the chosen distribution with that mean and variance \code{variance}.
#' Finally, all values are multiplied by \code{sigma_tau}.
#'
#' @param true_dist A character string specifying the type of distribution to use.
#'        Options are "Gaussian", "T", "Skew", "ALD", or "Mixture".
#' @param J A positive integer indicating the number of sites (default: 50). If \code{formula}
#'        and \code{data} are provided, \code{J} should match \code{nrow(data)} or can be
#'        left as default for consistency check.
#' @param sigma_tau A positive number to scale the generated treatment effects (default: 0.25).
#' @param tau A numeric value specifying an overall offset for the site means (default: 0).
#' @param variance A positive number specifying the baseline variance for the distribution (default: 1).
#' @param nu Degrees of freedom for the T distribution. Must be provided if \code{true_dist} is "T".
#' @param slant A numeric value for the skew parameter in the Skewed Normal distribution.
#'              Must be provided if \code{true_dist} is "Skew".
#' @param rho A numeric value for the parameter in the Asymmetric Laplace distribution.
#'            Must be provided if \code{true_dist} is "ALD".
#' @param delta A numeric value for the Mixture distribution.
#' @param eps A numeric value (typically between 0 and 1) for the Mixture distribution.
#' @param ups A positive numeric value for the Mixture distribution.
#' @param formula An optional model formula describing how covariates in \code{data} should be used
#'                to form the design matrix. For example, \code{~ X1 + X2} or \code{~ 0 + X1 + X2}.
#'                If provided, you must also supply \code{beta} and \code{data}.
#' @param beta A numeric vector of regression coefficients corresponding to the columns of
#'             \code{model.matrix(formula, data)}. If \code{formula} is provided, \code{beta} must
#'             have length equal to \code{ncol(model.matrix(...))}.
#' @param data An optional data frame containing site-level covariates to be used in \code{formula}.
#'             Must have \code{nrow(data) == J} if provided.
#' @param seed An optional integer for setting the random seed (default: 123). If \code{NULL}, seed is not set.
#' @param set_seed Logical indicating whether to set the seed (default: TRUE).
#'
#' @return A numeric vector of length \code{J} containing the generated treatment effects.
#'
#' @examples
#' \dontrun{
#'   # 1) No covariates (Gaussian distribution, just a constant mean = tau)
#'   effects1 <- gen_priorG2("Gaussian", J = 5, tau = 2, variance = 1)
#'   effects1
#'
#'   # 2) Using formula & data (Gaussian distribution)
#'   #    Suppose we have 2 covariates (X1, X2) for 5 sites
#'   site_data <- data.frame(
#'     X1 = c(1, 2, 3, 4, 5),
#'     X2 = c(0, 1, 0, 1, 0)
#'   )
#'   # Model formula includes an intercept by default (~ X1 + X2)
#'   # => design matrix has columns: (Intercept), X1, X2
#'   # => length(beta) must be 3
#'   betas <- c(10, 1, -2)
#'   effects2 <- gen_priorG2("Gaussian", formula = ~ X1 + X2, beta = betas, data = site_data)
#'   effects2
#'
#'   # 3) T distribution with covariates
#'   effects3 <- gen_priorG2("T", formula = ~ 0 + X1 + X2, beta = c(1, -2), data = site_data,
#'                           nu = 5, variance = 2)
#'   effects3
#'
#'   # 4) Asymmetric Laplace (requires LaplacesDemon)
#'   #    We'll also use tau offset
#'   if (requireNamespace("LaplacesDemon", quietly = TRUE)) {
#'     effects4 <- gen_priorG2("ALD", formula = ~ X1 + X2, beta = betas, data = site_data,
#'                             rho = 2, tau = 5)
#'     effects4
#'   }
#' }
#'
#' @export
gen_priorG2 <- function(true_dist,
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
                        set_seed = TRUE) {

  # 1) Set seed if requested
  if (!is.null(seed) && set_seed) {
    set.seed(seed)
  }

  # 2) Basic checks
  if (!is.numeric(J) || J <= 0 || floor(J) != J) {
    stop("J must be a positive integer.")
  }
  if (!is.numeric(sigma_tau) || sigma_tau <= 0) {
    stop("sigma_tau must be a positive number.")
  }
  if (!is.numeric(variance) || variance <= 0) {
    stop("variance must be a positive number.")
  }

  # 3) Handle formula/data approach if provided
  #    If formula is given, then data & beta are required.
  site_means <- rep(tau, J)  # default: no covariates => all sites have mean = tau
  if (!is.null(formula)) {
    if (is.null(data) || !is.data.frame(data)) {
      stop("If 'formula' is provided, 'data' must be a valid data.frame.")
    }
    mm <- model.matrix(formula, data = data)
    if (nrow(mm) != J) {
      stop("Number of rows in model.matrix(formula, data) must match J.")
    }
    if (is.null(beta)) {
      stop("If 'formula' is provided, 'beta' must be provided as well.")
    }
    if (!is.numeric(beta)) {
      stop("'beta' must be a numeric vector.")
    }
    if (length(beta) != ncol(mm)) {
      stop("Length of 'beta' must match the number of columns in model.matrix(formula, data).")
    }
    # site-wise offset from the formula
    offset <- as.vector(mm %*% beta)
    # combine with tau offset
    site_means <- site_means + offset
  }

  # 4) Distribution-specific internal functions
  gen_gaussian <- function(means, var) {
    rnorm(length(means), mean = means, sd = sqrt(var))
  }

  gen_t <- function(means, var, nu) {
    if (is.null(nu) || !is.numeric(nu) || nu <= 2) {
      stop("For T distribution, 'nu' must be a numeric value > 2.")
    }
    std_t <- rt(length(means), df = nu) * sqrt((nu - 2) / nu)
    means + sqrt(var) * std_t
  }

  gen_skew <- function(means, var, slant) {
    if (is.null(slant)) {
      stop("For Skewed Normal distribution, 'slant' must be provided.")
    }
    if (!requireNamespace("sn", quietly = TRUE)) {
      stop("Package 'sn' is required for the Skew distribution. Please install it.")
    }
    out <- numeric(length(means))
    for (i in seq_along(means)) {
      delta_s <- slant / sqrt(1 + slant^2)
      scale_s <- sqrt(var / (1 - (2 * delta_s^2 / pi)))
      xi_i <- means[i] - scale_s * sqrt(2 / pi) * delta_s
      out[i] <- sn::rsn(n = 1, xi = xi_i, omega = scale_s, alpha = slant)
    }
    out
  }

  gen_ald <- function(means, var, rho) {
    if (is.null(rho)) {
      stop("For ALD distribution, 'rho' must be provided.")
    }
    if (!requireNamespace("LaplacesDemon", quietly = TRUE)) {
      stop("Package 'LaplacesDemon' is required for the ALD distribution. Please install it.")
    }
    out <- numeric(length(means))
    for (i in seq_along(means)) {
      scale_a <- sqrt((2 * rho^2 * var) / (1 + rho^4))
      location_a <- means[i] - (scale_a * (1 / rho - rho)) / sqrt(2)
      out[i] <- LaplacesDemon::ralaplace(n = 1, location = location_a, scale = scale_a, kappa = rho)
    }
    out
  }

  gen_mixture <- function(means, var, delta, eps, ups) {
    if (is.null(delta) || is.null(eps) || is.null(ups)) {
      stop("For Mixture distribution, 'delta', 'eps', and 'ups' must be provided.")
    }
    out <- numeric(length(means))
    for (i in seq_along(means)) {
      a <- sqrt((1 - eps) + eps * ups^2 + eps * (1 - eps) * delta^2)
      ind <- runif(1) < (1 - eps)
      m <- if (ind) {
        -eps * delta / a
      } else {
        (1 - eps) * delta / a
      }
      s <- if (ind) {
        1 / a
      } else {
        ups / a
      }
      X_i <- rnorm(1, mean = m, sd = s)
      out[i] <- means[i] + sqrt(var) * X_i
    }
    out
  }

  # 5) Generate site-specific treatment effects
  effects <- switch(
    true_dist,
    "Gaussian" = gen_gaussian(site_means, variance),
    "T"        = gen_t(site_means, variance, nu),
    "Skew"     = gen_skew(site_means, variance, slant),
    "ALD"      = gen_ald(site_means, variance, rho),
    "Mixture"  = gen_mixture(site_means, variance, delta, eps, ups),
    stop("Invalid true_dist. Choose from 'Gaussian', 'T', 'Skew', 'ALD', or 'Mixture'.")
  )

  # 6) Rescale by sigma_tau
  return(effects * sigma_tau)
}
