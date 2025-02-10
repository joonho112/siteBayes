#' Generate Site-Specific Treatment Effects from a Prior Distribution
#'
#' This function generates site-specific treatment effects based on the specified
#' distribution and parameters. The generated values are scaled by \code{sigma_tau}.
#'
#' For each distribution, the function is parameterized so that the generated values
#' have the specified location (\code{tau}) and baseline variance (\code{variance}).
#' Finally, the entire vector is multiplied by \code{sigma_tau}. That is, the final
#' treatment effects have mean \code{tau * sigma_tau} and variance \code{variance * sigma_tau^2}.
#'
#' @param true_dist A character string specifying the type of distribution to use.
#'        Options are "Gaussian", "T", "Skew", "ALD", or "Mixture".
#' @param J A positive integer indicating the number of sites (default: 50).
#' @param sigma_tau A positive number to scale the generated treatment effects (default: 0.25).
#' @param tau A numeric value specifying the baseline mean for the distribution (default: 0).
#' @param variance A positive number specifying the baseline variance for the distribution (default: 1).
#' @param nu Degrees of freedom for the T distribution. Must be provided if \code{true_dist} is "T".
#' @param slant A numeric value for the skew parameter in the Skewed Normal distribution.
#'              Must be provided if \code{true_dist} is "Skew".
#' @param rho A numeric value for the parameter in the Asymmetric Laplace distribution.
#'            Must be provided if \code{true_dist} is "ALD".
#' @param delta A numeric value for the Mixture distribution.
#' @param eps A numeric value (typically between 0 and 1) for the Mixture distribution.
#' @param ups A positive numeric value for the Mixture distribution.
#' @param seed An optional integer for setting the random seed (default: 123). If \code{NULL}, seed is not set.
#' @param set_seed Logical indicating whether to set the seed (default: TRUE).
#'
#' @return A numeric vector of length \code{J} containing the generated treatment effects.
#'
#' @examples
#' \dontrun{
#'   # Gaussian distribution with default parameters
#'   gen_priorG("Gaussian", J = 50, sigma_tau = 0.25, tau = 0, variance = 1)
#'
#'   # T distribution with 5 degrees of freedom
#'   gen_priorG("T", J = 50, sigma_tau = 0.25, tau = 0, variance = 1, nu = 5)
#'
#'   # Skewed Normal distribution
#'   gen_priorG("Skew", J = 50, sigma_tau = 0.25, tau = 0, variance = 1, slant = 5)
#'
#'   # Asymmetric Laplace distribution
#'   gen_priorG("ALD", J = 50, sigma_tau = 0.25, tau = 0, variance = 1, rho = 2)
#'
#'   # Mixture distribution
#'   gen_priorG("Mixture", J = 50, sigma_tau = 0.25, tau = 0, variance = 1,
#'              delta = 1, eps = 0.3, ups = 2)
#' }
#' @export
gen_priorG <- function(true_dist,
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
                       seed = 123,
                       set_seed = TRUE) {
  # Set seed if provided and requested
  if (!is.null(seed) && set_seed) {
    set.seed(seed)
  }

  # Input validation
  if (!is.numeric(J) || J <= 0 || floor(J) != J) {
    stop("J must be a positive integer.")
  }
  if (!is.numeric(sigma_tau) || sigma_tau <= 0) {
    stop("sigma_tau must be a positive number.")
  }
  if (!is.numeric(variance) || variance <= 0) {
    stop("variance must be a positive number.")
  }

  # Helper function: Gaussian distribution (standard normal scaled and shifted)
  gen_gaussian <- function(J, tau, variance) {
    tau + sqrt(variance) * rnorm(J)
  }

  # Helper function: T distribution (standardized to have unit variance)
  gen_t <- function(J, tau, variance, nu) {
    if (is.null(nu) || !is.numeric(nu) || nu <= 2) {
      stop("For T distribution, nu must be provided as a numeric value greater than 2.")
    }
    # rt(J, df = nu) has variance nu/(nu-2); multiplying by sqrt((nu-2)/nu) standardizes it.
    tau + sqrt(variance) * (rt(J, df = nu) * sqrt((nu - 2) / nu))
  }

  # Helper function: Skewed Normal distribution using the 'sn' package
  gen_skew <- function(J, tau, variance, slant) {
    if (is.null(slant)) {
      stop("slant must be provided for the Skewed Normal distribution.")
    }
    if (!requireNamespace("sn", quietly = TRUE)) {
      stop("Package 'sn' is required for the Skew distribution. Please install it.")
    }
    # Compute parameters so that the resulting distribution has mean 'tau' and variance 'variance'
    delta_s <- slant / sqrt(1 + slant^2)
    scale_s <- sqrt(variance / (1 - (2 * delta_s^2 / pi)))
    location_s <- tau - scale_s * sqrt(2 / pi) * delta_s
    sn::rsn(n = J, xi = location_s, omega = scale_s, alpha = slant)
  }

  # Helper function: Asymmetric Laplace distribution using the 'LaplacesDemon' package
  gen_ald <- function(J, tau, variance, rho) {
    if (is.null(rho)) {
      stop("rho must be provided for the Asymmetric Laplace distribution.")
    }
    if (!requireNamespace("LaplacesDemon", quietly = TRUE)) {
      stop("Package 'LaplacesDemon' is required for the ALD distribution. Please install it.")
    }
    scale_a <- sqrt((2 * rho^2 * variance) / (1 + rho^4))
    location_a <- tau - (scale_a * (1 / rho - rho)) / sqrt(2)
    LaplacesDemon::ralaplace(n = J, location = location_a, scale = scale_a, kappa = rho)
  }

  # Helper function: Mixture distribution (mixture of two normals)
  gen_mixture <- function(J, tau, variance, delta, eps, ups) {
    if (is.null(delta) || is.null(eps) || is.null(ups)) {
      stop("delta, eps, and ups must be provided for the Mixture distribution.")
    }
    a <- sqrt((1 - eps) + eps * ups^2 + eps * (1 - eps) * delta^2)
    ind <- runif(J) < (1 - eps)
    # Generate a mixture of normals that (by construction) has zero mean and unit variance.
    X <- rnorm(J,
               mean = ifelse(ind, -eps * delta / a, (1 - eps) * delta / a),
               sd = ifelse(ind, 1 / a, ups / a))
    tau + sqrt(variance) * X
  }

  # Generate site-specific treatment effects based on the chosen distribution
  tau_j <- switch(true_dist,
                  "Gaussian" = gen_gaussian(J, tau, variance),
                  "T"        = gen_t(J, tau, variance, nu),
                  "Skew"     = gen_skew(J, tau, variance, slant),
                  "ALD"      = gen_ald(J, tau, variance, rho),
                  "Mixture"  = gen_mixture(J, tau, variance, delta, eps, ups),
                  stop("Invalid true_dist value. Choose from 'Gaussian', 'T', 'Skew', 'ALD', or 'Mixture'.")
  )

  # Rescale the generated values by sigma_tau
  return(tau_j * sigma_tau)
}
