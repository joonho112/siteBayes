#' @title Prior distributions for the Rubin Normal partial pooling model
#'
#' @name priors_rubin_normal
#' @description This file provides a list of possible distributions
#' that can be used to specify priors in a Rubin Normal partial pooling model,
#' as well as helper functions to interpret these distributions into a Stan-friendly format.
#'
#' We allow the user to write the priors in a natural syntax
#' (e.g. \code{normal(0,10)}, \code{cauchy(0,5)}) and then
#' transform them internally into an integer-coded family and parameter vector
#' for Stan. The partial pooling model typically needs priors for:
#' - \code{mu} (the overall mean)
#' - \code{tau} (the between-site standard deviation)
#' - \code{beta} (fixed effects coefficients, if any)
#'
#' Currently supported distributions include:
#' \itemize{
#'   \item \code{normal()}
#'   \item \code{cauchy()}
#'   \item \code{lognormal()}
#'   \item \code{uniform()}
#'   \item \code{student_t()}
#'   \item \code{half_student_t()} (for strictly positive scale)
#'   \item \code{half_cauchy()} (for strictly positive scale)
#' }
#' Also included are placeholders for \code{multinormal()} and \code{lkj()}
#' if you want to handle correlation structures or vector hypermeans.
#'
#' @details
#' - The function \code{interpret_prior()} will parse the user-provided prior object and return
#'   an integer family ID plus a 3-element numeric vector. The partial pooling Stan code
#'   uses those to define the actual priors inside \code{model\{\}} or \code{transformed parameters\{\}}.
#' - The function \code{prepare_prior_settings_rubin_normal()} merges \code{hypermean}, \code{hypersd}, and
#'   \code{beta} prior objects for passing along to \code{rubin_normal.stan} data.
#'
#' @author
#'   The original baggr design was by Witold Wiecek and Rachael Meager.
#'   Adapted here for the Rubin Normal model with additional families
#'   (half-student-t, half-cauchy).
#'
#' @references
#'   Gelman, Andrew. (2006). "Prior distributions for variance parameters in hierarchical models
#'   (comment on article by Browne and Draper)."
#'   \emph{Bayesian Analysis} 1(3), 515-534.
#'
#'   Lewandowski, Daniel, Dorota Kurowicka, and Harry Joe.
#'   "Generating Random Correlation Matrices Based on Vines and Extended Onion Method."
#'   \emph{Journal of Multivariate Analysis} 100, no. 9 (October 1, 2009): 1989-2001.
#'
#' @examples
#' normal_prior(0,10)
#' cauchy_prior(0,5)
#' half_student_t_prior(df=3, scale=2)
#' half_cauchy_prior(scale=5)
#'
#' custom_priors <- prepare_prior_settings_rubin_normal(normal_prior(0,10),
#'                                                      half_cauchy_prior(0,5),
#'                                                      normal_prior(0,10))
#' str(custom_priors)
NULL

########################
# 0) Basic checks
########################

check_scalar <- function(x) {
  if(length(x) != 1)
    stop("Argument must be scalar")
  if(!is.numeric(x))
    stop("Argument must be numeric")
}

########################
# 1) User-facing prior constructors
########################

#' Create a normal prior object
#'
#' @param loc Mean (location) of the normal distribution
#' @param scale Standard deviation of the normal distribution (must be > 0)
#'
#' @return A list with elements dist="normal" and values=c(loc, scale).
#' @export
normal_prior <- function(loc, scale) {
  check_scalar(loc)
  check_scalar(scale)
  if(scale <= 0)
    stop("Scale must be positive for normal_prior()")
  list(dist = "normal", values = c(loc, scale))
}


#' Create a cauchy prior object
#'
#' @param loc location parameter
#' @param scale scale parameter (must be > 0)
#'
#' @return A list with dist="cauchy" and values=c(loc, scale)
#' @export
cauchy_prior <- function(loc, scale) {
  check_scalar(loc)
  check_scalar(scale)
  if(scale <= 0)
    stop("Scale must be positive for cauchy_prior()")
  list(dist = "cauchy", values = c(loc, scale))
}


#' Create a lognormal prior object
#'
#' @param mu mean of ln(X) for lognormal
#' @param sigma standard deviation of ln(X) for lognormal, must be >0
#'
#' @return A list with dist="lognormal" and values=c(mu, sigma)
#' @export
lognormal_prior <- function(mu, sigma) {
  check_scalar(mu)
  check_scalar(sigma)
  if(sigma <= 0)
    stop("sigma must be positive for lognormal_prior()")
  list(dist = "lognormal", values = c(mu, sigma))
}


#' Create a uniform prior object
#'
#' @param lower lower bound
#' @param upper upper bound
#'
#' @return A list with dist="uniform" and values=c(lower, upper)
#' @export
uniform_prior <- function(lower, upper) {
  check_scalar(lower)
  check_scalar(upper)
  if(lower > upper)
    stop("In uniform_prior(a,b), a <= b is required.")
  list(dist = "uniform", values=c(lower, upper))
}


#' Create a Student-t prior object
#'
#' @param df degrees of freedom (>0)
#' @param loc location
#' @param scale scale (>0)
#'
#' @return A list with dist="student_t" and values=c(df, loc, scale)
#' @export
student_t_prior <- function(df, loc, scale) {
  check_scalar(df)
  check_scalar(loc)
  check_scalar(scale)
  if(df <= 0)
    stop("degrees of freedom must be >0")
  if(scale <= 0)
    stop("scale must be >0")
  list(dist="student_t", values=c(df, loc, scale))
}


#' Create a half-Student-t prior for scale parameters
#'
#' Often recommended for hierarchical standard deviations like tau (Gelman 2006).
#'
#' @param df degrees of freedom
#' @param scale scale parameter
#'
#' @details Internally you can treat half-Student-t by bounding a Student-t at 0.
#' For your Stan code, you might do real<lower=0> tau;  tau ~ student_t(...);
#' or an explicit T[0,] truncation. This function just signals to the code that
#' it's a half_t.
#'
#' @return A list with dist="half_student_t" and values=c(df, scale, 0)
#' @export
half_student_t_prior <- function(df, scale) {
  check_scalar(df)
  check_scalar(scale)
  if(df <= 0)
    stop("degrees of freedom must be > 0 for half_student_t_prior()")
  if(scale <= 0)
    stop("scale must be > 0 for half_student_t_prior()")
  list(dist = "half_student_t", values=c(df, scale, 0))
}


#' Create a half-Cauchy prior for scale parameters
#'
#' This is another popular choice for scale parameters in hierarchical models,
#' recommended for example in Gelman (2006).
#'
#' @param loc Typically 0 if it's strictly positive and centered at 0
#' @param scale Must be >0
#'
#' @return A list with dist="half_cauchy" and values=c(loc, scale)
#' @export
half_cauchy_prior <- function(loc=0, scale=5) {
  check_scalar(loc)
  check_scalar(scale)
  if(scale <=0)
    stop("scale must be >0 for half_cauchy_prior()")
  list(dist="half_cauchy", values=c(loc, scale, 0))
}


########################
# 2) Extended distributions for vectors/correlations
########################

#' A placeholder for multinormal prior
#'
#' This can handle vector hypermeans. Typically you won't use it for scalar mu/tau,
#' but for multi-dimensional hypermeans if your model is extended.
#'
#' @param location numeric vector of means
#' @param Sigma variance-covariance matrix
#' @return A list with dist="multinormal", mean=location, scale=Sigma
#' @export
multinormal <- function(location, Sigma) {
  if(length(location) == 1)
    stop("For a single dimension, use normal() instead.")
  x <- try(chol(Sigma), silent = TRUE)
  if(inherits(x, "try-error") || !isSymmetric(Sigma))
    stop("Sigma must be positive semi-definite & symmetric")
  if(ncol(Sigma) != length(location))
    stop("Dimensions of location and Sigma must match")
  list(dist="multinormal", mean=location, scale=Sigma, dimension=length(location))
}


#' An LKJ prior for correlation matrices
#'
#' @param shape shape parameter for the LKJ distribution
#' @param order dimension of the correlation matrix (may be optional)
#'
#' @return list with dist="lkj", values=c(shape)
#' @export
lkj <- function(shape, order=NULL) {
  check_scalar(shape)
  list(dist="lkj", values=c(shape), dimension=order)
}


########################
# 3) The interpret_prior & prepare_prior_settings_rubin_normal
########################

#' Interpret a user-facing prior object into an integer-coded family
#'
#' This function is used internally to transform e.g. \code{list(dist="normal", values=c(0,10))}
#' into e.g. \code{family_id=1, param_vector=c(0,10,0)} for the Stan model to parse.
#'
#' @param prior_obj a list from e.g. normal_prior(), cauchy_prior(), etc.
#'
#' @return a list with \code{family_id} and \code{param_vector}
#' @keywords internal
interpret_prior <- function(prior_obj) {
  d <- prior_obj$dist
  v <- prior_obj$values

  # We can assign codes, e.g.:
  # 1=normal,2=cauchy,5=lognormal,0=uniform,
  # 3=student_t, 6=half_student_t, 7=half_cauchy
  # 9=multinormal, 10=lkj
  if(d == "normal"){
    return(list(family_id=1, param_vector=c(v[1], v[2], 0)))
  } else if(d == "cauchy"){
    return(list(family_id=2, param_vector=c(v[1], v[2], 0)))
  } else if(d == "lognormal"){
    return(list(family_id=5, param_vector=c(v[1], v[2], 0)))
  } else if(d == "uniform"){
    return(list(family_id=0, param_vector=c(v[1], v[2], 0)))
  } else if(d == "student_t"){
    return(list(family_id=3, param_vector=c(v[1], v[2], v[3])))
  } else if(d == "half_student_t"){
    return(list(family_id=6, param_vector=c(v[1], v[2], 0)))
  } else if(d == "half_cauchy"){
    return(list(family_id=7, param_vector=c(v[1], v[2], 0)))
  } else if(d == "multinormal"){
    return(list(family_id=9, param_vector=c(0,0,0)))
  } else if(d == "lkj"){
    return(list(family_id=10, param_vector=c(v[1],0,0)))
  } else {
    stop("Unrecognized distribution type: ", d)
  }
}


#' Merge hypermean, hypersd, and beta prior objects for Rubin Normal partial pooling
#'
#' This function merges the three prior objects (for overall mean \code{mu},
#' between-site SD \code{tau}, and optional fixed-effect \code{beta}) into
#' a single list that can be directly passed to the \code{rubin_normal.stan} model.
#'
#' @param hypermean A prior object for the overall mean \(\mu\).
#'   Defaults to \code{normal_prior(0,10)} (a fairly wide normal).
#' @param hypersd A prior object for the between-site standard deviation \(\tau\).
#'   Defaults to \code{half_student_t_prior(df=3, scale=2)},
#'   which is a commonly recommended weakly-informative choice (Gelman 2006).
#' @param beta A prior object for fixed-effect coefficients.
#'   Defaults to \code{normal_prior(0,10)}.
#'
#' @details
#' **Recommended defaults**:
#' \itemize{
#'   \item \strong{hypermean (mu)}: \code{normal_prior(0,10)}
#'   \item \strong{hypersd (tau)}:  \code{half_student_t_prior(df=3, scale=2)}
#'   \item \strong{beta}: \code{normal_prior(0,10)}
#' }
#' These reflect a weakly informative approach often used in Bayesian hierarchical
#' models (Gelman 2006).
#'
#' You can override these defaults with your own prior objects
#' (e.g. \code{cauchy_prior(0,5)}, \code{half_cauchy_prior(0,10)},
#' \code{student_t_prior(4,0,5)}, etc.).
#'
#' This function calls \code{\link{interpret_prior}} internally to transform
#' each prior object into an integer-coded family plus a parameter vector recognized
#' by the \code{rubin_normal.stan} code.
#'
#' @return A named list with fields:
#' \itemize{
#'   \item \code{prior_hypermean_fam}, \code{prior_hypermean_val}
#'   \item \code{prior_hypersd_fam},   \code{prior_hypersd_val}
#'   \item \code{prior_beta_fam},      \code{prior_beta_val}
#' }
#'
#' @seealso \code{\link{normal_prior}}, \code{\link{cauchy_prior}},
#'   \code{\link{half_student_t_prior}}, \code{\link{interpret_prior}}
#'
#' @examples
#' # Using all defaults
#' default_priors <- prepare_prior_settings_rubin_normal()
#' str(default_priors)
#'
#' # Overriding just the hypersd prior
#' alternative <- prepare_prior_settings_rubin_normal(
#'   hypersd = half_cauchy_prior(loc=0, scale=5)
#' )
#' str(alternative)
#'
#' @export
prepare_prior_settings_rubin_normal <- function(
    hypermean = normal_prior(0, 10),
    hypersd   = half_student_t_prior(df=3, scale=2),
    beta      = normal_prior(0, 10)
) {
  hm <- interpret_prior(hypermean)
  hs <- interpret_prior(hypersd)
  b  <- interpret_prior(beta)

  list(
    prior_hypermean_fam = hm$family_id,
    prior_hypermean_val = hm$param_vector,
    prior_hypersd_fam   = hs$family_id,
    prior_hypersd_val   = hs$param_vector,
    prior_beta_fam      = b$family_id,
    prior_beta_val      = b$param_vector
  )
}


#' Generate data-dependent priors for Rubin Normal model (baggr-style), including se2_col checks
#'
#' This function examines the site-level effect estimates (\code{tau_hat_col}) in \code{df},
#' and optionally any covariates (\code{covariates}), to produce three prior objects:
#' \code{hypermean}, \code{hypersd}, and \code{beta}. The resulting list can be supplied to
#' \code{\link{prepare_prior_settings_rubin_normal}()} or \code{\link{build_standata_rubin_normal}()}.
#'
#' The logic here is adapted from the baggr package (for the "rubin" model):
#' \itemize{
#'   \item \code{hypermean} (\(\mu\)) = \(\text{Normal}(0,\;10 \times \max|\tau_j|)\)
#'   \item \code{hypersd} (\(\tau\)) = \(\text{Uniform}(0,\;10 \times \mathrm{sd}(\tau_j))\)
#'         (if \code{J>1}). If there is only one site, a placeholder (\(\text{Normal}(0,1)\)) is used.
#'   \item \code{beta} (covariate effects) = \(\text{Normal}(0,\;10 \times \max \mathrm{sd}(\mathbf{X}))\)
#'         if covariates are present, or \(\text{Uniform}(0,1)\) if none.
#' }
#'
#' Although this function also checks for \code{se2_col}, it does not incorporate the
#' squared standard errors into the automatic prior logic here (consistent with baggr).
#'
#' @param df A data frame containing columns for site-level effects and their squared standard errors
#' @param tau_hat_col Character name of the column for observed site-level effects \(\tau_j\).
#'   Defaults to \code{"tau_j_hat"}.
#' @param se2_col Character name of the column for squared standard errors. Defaults to \code{"se2_hat"}.
#'   It is checked for existence but not otherwise used in prior calculations.
#' @param covariates Character vector of covariate column names, or \code{NULL} if none.
#'
#' @return A named list with three elements:
#' \itemize{
#'   \item \code{hypermean}: Prior object for the hypermean \(\mu\).
#'   \item \code{hypersd}:   Prior object for the hyper-SD \(\tau\).
#'   \item \code{beta}:      Prior object for covariate effects.
#' }
#'
#' @examples
#' \dontrun{
#' # Example data frame with site-level effects and squared SE, plus one covariate
#' df_example <- data.frame(
#'   tau_j_hat = c(0.10, 0.30, -0.20, 0.00),
#'   se2_hat   = c(0.01, 0.02, 0.015, 0.005),
#'   X1        = c(1.2, 0.9, 1.1, 1.3)
#' )
#'
#' # Generate data-dependent priors
#' auto_priors <- auto_prior_rubin_normal(
#'   df          = df_example,
#'   tau_hat_col = "tau_j_hat",
#'   se2_col     = "se2_hat",
#'   covariates  = "X1"
#' )
#' str(auto_priors)
#'
#' # Then pass them to build_standata_rubin_normal (for instance):
#' standata <- build_standata_rubin_normal(
#'   df              = df_example,
#'   tau_hat_col     = "tau_j_hat",
#'   se2_col         = "se2_hat",
#'   covariates      = "X1",
#'   prior_hypermean = auto_priors$hypermean,
#'   prior_hypersd   = auto_priors$hypersd,
#'   prior_beta      = auto_priors$beta
#' )
#'
#' # Now 'standata' is ready for rstan::sampling(...)
#' }
#'
#' @export
auto_prior_rubin_normal <- function(df,
                                    tau_hat_col = "tau_j_hat",
                                    se2_col     = "se2_hat",
                                    covariates  = NULL) {
  # 1) Check that tau_hat_col exists
  if (!tau_hat_col %in% names(df)) {
    stop("Column '", tau_hat_col, "' not found in data.")
  }

  # 2) Check that se2_col exists (even though we don't use it in the prior logic)
  if (!se2_col %in% names(df)) {
    stop("Column '", se2_col, "' not found in data.")
  }

  # 3) Extract site-level effect estimates
  tau_vals <- df[[tau_hat_col]]
  J <- length(tau_vals)

  # 4) Prepare the covariate matrix if specified
  if (!is.null(covariates) && length(covariates) > 0) {
    missing_cov <- setdiff(covariates, names(df))
    if (length(missing_cov) > 0) {
      stop("Missing covariate columns: ", paste(missing_cov, collapse=", "))
    }
    X  <- as.matrix(df[covariates])
    Nc <- ncol(X)
  } else {
    X  <- NULL
    Nc <- 0
  }

  # 5) hypermean: Normal(0, 10 * max(|tau_j|))
  max_tau <- max(abs(tau_vals), na.rm=TRUE)
  if (max_tau < 1e-8) {
    # fallback if near-zero
    max_tau <- 1.0
  }
  auto_hypermean <- normal_prior(0, 10 * max_tau)

  # 6) hypersd: Uniform(0, 10 * sd(tau_j)), fallback if J=1
  if (J < 2) {
    # only 1 site => partial pooling not meaningful
    auto_hypersd <- normal_prior(0, 1)
  } else {
    tau_sd <- sd(tau_vals, na.rm=TRUE)
    if (tau_sd < 1e-8) {
      tau_sd <- 1.0
    }
    auto_hypersd <- uniform_prior(0, 10 * tau_sd)
  }

  # 7) beta: if covariates, Normal(0, 10 * max sd(X)); if none, Uniform(0,1)
  if (Nc > 0) {
    x_sd <- apply(X, 2, sd)
    if (all(x_sd < 1e-8)) {
      x_sd_val <- 10  # fallback
    } else {
      x_sd_val <- 10 * max(x_sd)
    }
    auto_beta <- normal_prior(0, x_sd_val)
  } else {
    # no covariates
    auto_beta <- uniform_prior(0, 1)
  }

  # 8) Return the prior objects in a named list
  list(
    hypermean = auto_hypermean,
    hypersd   = auto_hypersd,
    beta      = auto_beta
  )
}
