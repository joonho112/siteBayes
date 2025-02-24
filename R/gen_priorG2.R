#' Generate Site-Specific Treatment Effects from a Prior Distribution
#'
#' This function generates site-specific treatment effects based on the specified
#' distribution and parameters. The generated values are then scaled by \code{sigma_tau}.
#'
#' You can specify site-level covariates in two ways:
#' \enumerate{
#'   \item Provide a \code{formula} and a \code{data} frame, which will be processed via
#'         \code{model.matrix}. The resulting design matrix is multiplied by \code{beta}
#'         to shift each site's mean.
#'   \item Or simply set \code{tau} if you want a constant baseline mean with no additional covariates.
#' }
#'
#' The final site-specific means are
#' \code{site_means = tau + model_matrix \%*\% beta} (if \code{formula} is provided),
#' or simply \code{rep(tau, J)} if no \code{formula} and \code{beta} are given.
#' Then, for each site, a value is drawn from the chosen distribution
#' (e.g. Gaussian, T, skew-normal, etc.) with that mean and variance \code{variance}.
#' Finally, all values are multiplied by \code{sigma_tau}.
#'
#' The function \strong{returns a data frame} with:
#' \itemize{
#'   \item \code{site_index} (integer from \eqn{1} to \eqn{J})
#'   \item \code{site_effect} (the generated treatment effect for each site)
#'   \item If \code{formula} and \code{data} were provided, all original covariate columns
#'         appear to the \emph{right} of \code{site_effect}.
#'   \item Otherwise, only \code{site_index} and \code{site_effect}.
#' }
#'
#' @param true_dist A character string specifying the type of distribution to use.
#'        Options are \code{"Gaussian"}, \code{"T"}, \code{"Skew"}, \code{"ALD"}, or \code{"Mixture"}.
#' @param J A positive integer indicating the number of sites (default: 50). If \code{formula}
#'        and \code{data} are provided, then \code{J} must match \code{nrow(data)}
#'        (or be explicitly set to that).
#' @param sigma_tau A positive number to scale the generated treatment effects (default: 0.25).
#' @param tau A numeric value specifying an overall offset for the site means (default: 0).
#' @param variance A positive number specifying the baseline variance for the distribution (default: 1).
#' @param nu Degrees of freedom for the T distribution. Must be provided if \code{true_dist = "T"}.
#' @param slant A numeric value for the skew parameter in the Skew Normal distribution.
#'              Must be provided if \code{true_dist = "Skew"}.
#' @param rho A numeric value for the parameter in the Asymmetric Laplace distribution (ALD).
#'            Must be provided if \code{true_dist = "ALD"}.
#' @param delta A numeric value for the Mixture distribution (must be provided if \code{true_dist = "Mixture"}).
#' @param eps A numeric value (typically between 0 and 1) for the Mixture distribution.
#' @param ups A positive numeric value for the Mixture distribution.
#' @param formula An optional model formula describing how covariates in \code{data} should be used
#'                to form the design matrix. For example, \code{~ X1 + X2} or \code{~ 0 + X1 + X2}.
#'                If provided, you must also supply \code{beta} and \code{data}.
#' @param beta A numeric vector of regression coefficients corresponding to the columns of
#'             \code{model.matrix(formula, data)}. If \code{formula} is provided,
#'             \code{beta} must have length equal to \code{ncol(model.matrix(...))}.
#' @param data An optional data frame containing site-level covariates to be used in \code{formula}.
#'             Must have \code{nrow(data) == J} if provided.
#' @param seed An optional integer for setting the random seed (default: 123). If \code{NULL}, no seed is set.
#' @param set_seed Logical indicating whether to call \code{set.seed(seed)} (default: TRUE).
#'
#' @return A \strong{data frame} with \code{J} rows. It always contains
#' \code{site_index} and \code{site_effect}, plus covariate columns if you provided them.
#' Covariates appear to the \emph{right} of these two columns.
#'
#' @examples
#' \dontrun{
#'   ##################################################
#'   # Example 1) No covariates, Gaussian distribution
#'   ##################################################
#'   # A simple example: 5 sites, baseline mean=2, variance=1
#'   effects1_df <- gen_priorG2("Gaussian", J = 5, tau = 2, variance = 1)
#'   effects1_df
#'
#'   # Another distribution with no covariates:
#'   # "T" distribution, df=5, variance=2
#'   effectsT_df <- gen_priorG2("T", J=6, nu=5, variance=2)
#'   effectsT_df
#'
#'   # Or "Mixture" (a simple 2-component mixture):
#'   effMix_df <- gen_priorG2("Mixture", J=6, delta=0.5, eps=0.3, ups=1.0)
#'   effMix_df
#'
#'   ##################################################
#'   # Example 2) Covariates: binomial random X1,X2,X3
#'   ##################################################
#'   # Suppose we have 100 sites
#'   J <- 100
#'   p1 <- 0.5
#'   p2 <- 0.3
#'   p3 <- 0.7
#'
#'   X1 <- rbinom(J, size=1, prob=p1)
#'   X2 <- rbinom(J, size=1, prob=p2)
#'   X3 <- rbinom(J, size=1, prob=p3)
#'
#'   # Put covariates in a data.frame
#'   site_data <- data.frame(X1 = X1, X2 = X2, X3 = X3)
#'
#'   # formula = ~ X1 + X2 + X3 => design matrix has 4 columns: Intercept + X1 + X2 + X3
#'   beta_vec <- c(2, 1.5, -1, 0.5)
#'     # Intercept=2, X1=1.5, X2=-1, X3=0.5
#'
#'   # Now generate site effects from e.g. "Gaussian"
#'   effCov_df <- gen_priorG2(
#'     true_dist = "Gaussian",
#'     J         = J,
#'     formula   = ~ X1 + X2 + X3,
#'     beta      = beta_vec,
#'     data      = site_data,
#'     variance  = 1,
#'     sigma_tau = 0.5
#'   )
#'   head(effCov_df)
#'   # We see site_index, site_effect, then X1, X2, X3
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

  # 1) Possibly set random seed
  if (!is.null(seed) && set_seed) {
    set.seed(seed)
  }

  # 2) Basic argument checks
  if (!is.numeric(J) || J <= 0 || floor(J) != J) {
    stop("`J` must be a positive integer.")
  }
  if (!is.numeric(sigma_tau) || sigma_tau <= 0) {
    stop("`sigma_tau` must be a positive number.")
  }
  if (!is.numeric(variance) || variance <= 0) {
    stop("`variance` must be a positive number.")
  }

  # 3) If formula+data => must match row count
  using_covariates <- !is.null(formula)
  if (using_covariates) {
    if (is.null(data) || !is.data.frame(data)) {
      stop("If `formula` is provided, `data` must be a valid data.frame.")
    }
    if (nrow(data) != J) {
      stop("Number of rows in `data` must match `J` when `formula` is used.")
    }
  }

  # 4) Build site_means
  site_means <- rep(tau, J)
  if (using_covariates) {
    mm <- model.matrix(formula, data = data)
    if (is.null(beta)) {
      stop("If `formula` is provided, `beta` must be provided.")
    }
    if (!is.numeric(beta) || length(beta) != ncol(mm)) {
      stop("Length of `beta` must match the columns of model.matrix(formula, data).")
    }
    offset <- as.vector(mm %*% beta)
    site_means <- site_means + offset
  }

  # 5) Generate unscaled random deviate for each site
  unscaled_effects <- switch(
    true_dist,
    "Gaussian" = rnorm(J, mean = site_means, sd = sqrt(variance)),
    "T" = {
      if (is.null(nu) || nu <= 2) {
        stop("For `T`, argument `nu` must be numeric and > 2.")
      }
      devs <- rt(J, df = nu)
      site_means + sqrt(variance) * devs
    },
    "Skew" = {
      if (is.null(slant)) {
        stop("For `Skew`, argument `slant` must be provided.")
      }
      if (!requireNamespace("sn", quietly = TRUE)) {
        stop("Package `sn` is required for the Skew distribution.")
      }
      vapply(seq_len(J), function(i) {
        sn::rsn(n = 1, xi = site_means[i], omega = sqrt(variance), alpha = slant)
      }, numeric(1))
    },
    "ALD" = {
      if (is.null(rho)) {
        stop("For `ALD`, argument `rho` must be provided.")
      }
      if (!requireNamespace("LaplacesDemon", quietly = TRUE)) {
        stop("Package `LaplacesDemon` is required for the ALD distribution.")
      }
      scale_val <- sqrt(variance)
      vapply(seq_len(J), function(i) {
        LaplacesDemon::ralaplace(1, location = site_means[i],
                                 scale = scale_val, kappa = rho)
      }, numeric(1))
    },
    "Mixture" = {
      if (any(sapply(list(delta, eps, ups), is.null))) {
        stop("For `Mixture`, you must provide `delta`, `eps`, and `ups`.")
      }
      sdev <- sqrt(variance)
      vapply(seq_len(J), function(i) {
        if (runif(1) < (1 - eps)) {
          rnorm(1, mean = site_means[i] - delta, sd = sdev)
        } else {
          rnorm(1, mean = site_means[i] + ups, sd = sdev)
        }
      }, numeric(1))
    },
    stop("Invalid `true_dist`. Must be one of 'Gaussian','T','Skew','ALD','Mixture'.")
  )

  # 6) Scale by sigma_tau
  final_effects <- unscaled_effects * sigma_tau

  # 7) Construct final data frame
  if (using_covariates) {
    # If not present, add a site_index column
    if (!"site_index" %in% names(data)) {
      data$site_index <- seq_len(J)
    }
    # Add site_effect
    data$site_effect <- final_effects

    # Reorder columns: put site_index first, then site_effect second, then the rest
    old_cols <- names(data)
    # remove them from old_cols
    remainder <- setdiff(old_cols, c("site_index", "site_effect"))
    new_order <- c("site_index", "site_effect", remainder)
    data <- data[new_order]

    return(data)

  } else {
    # Minimal data frame: site_index & site_effect
    df_out <- data.frame(
      site_index  = seq_len(J),
      site_effect = final_effects
    )
    return(df_out)
  }
}
