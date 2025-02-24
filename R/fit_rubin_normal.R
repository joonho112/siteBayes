#' Fit the Rubin Normal partial pooling model via Stan
#'
#' This function:
#' \enumerate{
#'   \item Interprets the data frame (\code{df}) columns for site-level effects and SEs.
#'   \item Optionally sets priors automatically if \code{auto_prior=TRUE}.
#'   \item Builds a Stan-ready data list using \code{\link{build_standata_rubin_normal}}.
#'   \item Calls \code{\link{compile_rubin_normal_mod}} internally to compile the
#'         \code{rubin_normal.stan} model if not already compiled.
#'   \item Invokes \code{rstan::sampling} on that compiled model.
#'   \item Returns the resulting \code{stanfit} object.
#' }
#'
#' By default, this function uses \code{auto_prior=TRUE} with all three prior arguments
#' set to \code{NULL}, meaning it will automatically choose data-dependent priors
#' (baggr-style). If you wish to override any of them (\code{prior_hypermean},
#' \code{prior_hypersd}, \code{prior_beta}), supply a user-defined prior object
#' (e.g. \code{normal_prior(0,10)}) for that specific parameter.
#'
#' @param df A data frame containing the columns for site-level estimates (\code{tau_hat_col})
#'   and squared SEs (\code{se2_col}), plus optional covariate columns.
#' @param tau_hat_col Name of the column with site-level effect estimates (\(\tau_j\)).
#'   Default is \code{"tau_j_hat"}.
#' @param se2_col Name of the column with squared standard errors. Default is \code{"se2_hat"}.
#' @param covariates Character vector of covariate column names (or \code{NULL} if none).
#' @param prior_hypermean Optional prior object for the hyper-mean \(\mu\).
#'   If \code{NULL} and \code{auto_prior=TRUE}, will be set automatically.
#' @param prior_hypersd Optional prior object for the hyper-SD \(\tau\).
#'   If \code{NULL} and \code{auto_prior=TRUE}, will be set automatically.
#' @param prior_beta Optional prior object for covariate effects \(\beta\).
#'   If \code{NULL} and \code{auto_prior=TRUE}, will be set automatically.
#' @param auto_prior Logical. If \code{TRUE} (default), any \code{NULL} priors
#'   are chosen automatically using \code{\link{auto_prior_rubin_normal}} logic.
#'   If \code{FALSE}, \code{NULL} priors revert to simpler defaults (e.g.,
#'   \code{normal_prior(0,10)}, etc.).
#' @param iter Number of MCMC iterations per chain. Default 2000.
#' @param chains Number of MCMC chains. Default 4.
#' @param cores Number of cores for parallel MCMC. Default 2.
#' @param ... Additional arguments passed on to \code{rstan::sampling}.
#'
#' @details
#' Internally, this function calls \code{\link{compile_rubin_normal_mod}} to compile
#' (or retrieve) the Stan model \code{rubin_normal.stan}. Once compiled, it caches
#' the resulting model object so that subsequent calls do not re-compile.
#'
#' @return A \code{stanfit} object (from the \pkg{rstan} package), containing the posterior draws.
#'
#' @examples
#' \dontrun{
#' # Example data
#' df_with_cov <- data.frame(
#'   tau_j_hat = c(0.1, 0.3, -0.2),
#'   se2_hat   = c(0.01, 0.02, 0.03),
#'   X1        = c(1, 2, 1)
#' )
#'
#' # 1) Default usage => auto_prior=TRUE with all priors NULL => data-dependent priors
#' fit_auto <- fit_rubin_normal(
#'   df          = df_with_cov,
#'   tau_hat_col = "tau_j_hat",
#'   se2_col     = "se2_hat",
#'   covariates  = "X1"
#' )
#'
#' # 2) Manually override just the hyper-SD prior
#' my_hypersd <- half_student_t_prior(df=4, scale=1.5)
#' fit_custom <- fit_rubin_normal(
#'   df             = df_with_cov,
#'   tau_hat_col    = "tau_j_hat",
#'   se2_col        = "se2_hat",
#'   covariates     = "X1",
#'   prior_hypersd  = my_hypersd  # override
#' )
#'
#' # 3) auto_prior=FALSE => anything left as NULL uses simpler defaults (normal(0,10), etc.)
#' fit_simple <- fit_rubin_normal(
#'   df          = df_with_cov,
#'   tau_hat_col = "tau_j_hat",
#'   se2_col     = "se2_hat",
#'   covariates  = "X1",
#'   auto_prior  = FALSE
#' )
#'
#' # Inspect results
#' print(fit_auto, pars=c("tau","sigma_tau"), probs=c(0.025,0.5,0.975))
#' }
#' @export
fit_rubin_normal <- function(
    df,
    tau_hat_col     = "tau_j_hat",
    se2_col         = "se2_hat",
    covariates      = NULL,
    prior_hypermean = NULL,
    prior_hypersd   = NULL,
    prior_beta      = NULL,
    auto_prior      = TRUE,
    iter            = 2000,
    chains          = 4,
    cores           = 2,
    ...
) {
  # 1) Build the Stan data list
  standata <- build_standata_rubin_normal(
    df               = df,
    tau_hat_col      = tau_hat_col,
    se2_col          = se2_col,
    covariates       = covariates,
    prior_hypermean  = prior_hypermean,
    prior_hypersd    = prior_hypersd,
    prior_beta       = prior_beta,
    auto_prior       = auto_prior
  )

  # 2) Compile or retrieve the rubin_normal.stan model
  mod <- compile_rubin_normal_mod()

  # 3) Run MCMC
  fit <- rstan::sampling(
    object  = mod,
    data    = standata,
    iter    = iter,
    chains  = chains,
    cores   = cores,
    ...
  )

  return(fit)
}
