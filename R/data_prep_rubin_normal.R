#' @title Data preparation for the Rubin Normal partial pooling model
#'
#' @name data_prep_rubin_normal
#'
#' @description
#' This file provides two functions for assembling summary-level data
#' for use in a \code{rubin_normal.stan} model:
#' \enumerate{
#'   \item \code{\link{prepare_summary_data_rubin}}: reads columns from a data frame
#'     containing site-level effect estimates (e.g., \code{tau_j_hat} or \code{theta_hat})
#'     and their squared standard errors (\code{se2_col}).
#'   \item \code{\link{build_standata_rubin_normal}}: merges these data with user-chosen
#'     priors (hypermean, hypersd, beta) into a single list ready for \code{rstan::sampling}.
#' }
#'
#' Typically, your data frame might already have site-level effect estimates
#' (a column like \code{"tau_j_hat"}) and squared standard errors
#' (a column like \code{"se2_j"}). This script extracts them, optionally collects
#' covariates, and converts the squared SEs to regular SEs by taking \code{sqrt}.
#'
#' @author
#'   (Your name / or team name)
#'
#' @seealso
#'   \code{\link{prepare_prior_settings_rubin_normal}} for merging prior objects,
#'   \code{\link{fit_rubin_normal}} for the actual fitting.
#'
#' @keywords internal
NULL


########################
# 1) Prepare summary-level data
########################

#' Prepare summary-level data for the rubin_normal model
#'
#' This function expects columns for:
#' - Observed site-level effects: \code{tau_hat_col} (default \code{"tau_j_hat"}).
#' - Squared standard errors: \code{se2_col} (default \code{"se2_hat"}).
#' - Optional covariates (\code{covariates}).
#'
#' If your data has actual SEs (not squared), you can either rename them to \code{"se2_hat"}
#' and square them externally, or just store them in the \code{se2_col} while setting
#' the right name. This function will do \code{sqrt(...)} on whatever \code{se2_col} is,
#' so be sure it truly represents \emph{squared} errors if you follow these defaults.
#'
#' @param df A data frame with summary-level data, typically including
#'   \code{tau_j_hat} (or your chosen column name) and \code{se2_hat}.
#' @param tau_hat_col Character: the column name in \code{df} for the site-level
#'   observed effect estimates (defaults to \code{"tau_j_hat"}).
#' @param se2_col Character: the column name in \code{df} for squared standard errors
#'   (defaults to \code{"se2_hat"}).
#' @param covariates Character vector of covariate column names (or \code{NULL} if none).
#'
#' @return A named list with components:
#' \itemize{
#'   \item \code{J}: number of sites
#'   \item \code{tau_hat_j}: numeric vector of site-level effect estimates
#'         (extracted from \code{df[[tau_hat_col]]})
#'   \item \code{se_j}: numeric vector of standard errors (\code{sqrt(df[[se2_col]])})
#'   \item \code{Nc}: number of covariates
#'   \item \code{X}: matrix of covariates (or a zero-column matrix if \code{covariates=NULL})
#' }
#'
#' @examples
#' # Example data frame might have columns: tau_j_hat, se2_hat, X1, X2, ...
#' df_with_cov <- data.frame(
#'   tau_j_hat = c(0.36, 0.317, 0.765),
#'   se2_hat   = c(0.167, 0.154, 0.222),
#'   X1        = c(0,1,1),
#'   X2        = c(0,1,0)
#' )
#' prep <- prepare_summary_data_rubin(
#'   df          = df_with_cov,
#'   tau_hat_col = "tau_j_hat",
#'   se2_col     = "se2_hat",
#'   covariates  = c("X1","X2")
#' )
#' str(prep)
#'
#' @export
prepare_summary_data_rubin <- function(df,
                                       tau_hat_col = "tau_j_hat",
                                       se2_col     = "se2_hat",
                                       covariates  = NULL) {
  # 1) Check presence of needed columns
  needed_cols <- c(tau_hat_col, se2_col)
  missing_cols <- setdiff(needed_cols, names(df))
  if(length(missing_cols) > 0) {
    stop(
      "Data is missing required column(s): ",
      paste(missing_cols, collapse=", ")
    )
  }

  # 2) Extract vectors
  tau_hat_j <- df[[tau_hat_col]]
  se2_j     <- df[[se2_col]]
  if(any(se2_j < 0, na.rm=TRUE)) {
    stop("Some squared SE values < 0 (invalid). Check column '", se2_col, "'.")
  }
  se_j <- sqrt(se2_j)

  J <- length(tau_hat_j)

  # 3) Covariate matrix
  if(!is.null(covariates) && length(covariates) > 0) {
    if(!all(covariates %in% names(df))) {
      stop("Some covariate columns not found in data: ",
           paste(setdiff(covariates, names(df)), collapse=", "))
    }
    X <- as.matrix(df[covariates])
  } else {
    X <- matrix(0, nrow = J, ncol = 0)  # no covariates
  }

  # 4) Return assembled list
  list(
    J          = J,
    tau_hat_j  = tau_hat_j,
    se_j       = se_j,
    Nc         = ncol(X),
    X          = X
  )
}


########################
# 2) Build Stan data list for rubin_normal model
########################

#' Build Stan data list for rubin_normal.stan (auto priors by default)
#'
#' This function returns a list suitable for \code{rstan::sampling(...)}
#' with the Rubin Normal model. It expects a data frame containing site-level
#' effect estimates (\code{tau_hat_col}) and squared standard errors (\code{se2_col}),
#' plus optional covariates.
#'
#' By default, \code{auto_prior = TRUE}, and all three prior arguments
#' (\code{prior_hypermean}, \code{prior_hypersd}, \code{prior_beta})
#' are \code{NULL}. In that scenario, the function will automatically
#' generate data-dependent priors for each parameter using
#' \code{\link{auto_prior_rubin_normal}()} (baggr-style).
#'
#' If you provide a manual prior for some parameter(s), those override
#' the auto priors for that portion only. If you set \code{auto_prior = FALSE},
#' then any \code{NULL} prior argument reverts to a simple default
#' (e.g. \code{normal_prior(0, 10)}, \code{half_student_t_prior(df=3, scale=2)})
#' as in the earlier behavior.
#'
#' @param df A data frame containing at least the columns \code{tau_hat_col} and
#'   \code{se2_col}, and optionally some covariates.
#' @param tau_hat_col Name of the column with site-level effect estimates
#'   (defaults to \code{"tau_j_hat"}).
#' @param se2_col Name of the column with squared standard errors
#'   (defaults to \code{"se2_hat"}).
#' @param covariates Character vector of covariate column names, or \code{NULL} if none.
#' @param prior_hypermean A prior object for the hyper-mean \(\mu\). If \code{NULL},
#'   it will be set automatically if \code{auto_prior=TRUE}, or to
#'   \code{normal_prior(0, 10)} otherwise.
#' @param prior_hypersd A prior object for the hyper-SD \(\tau\). If \code{NULL},
#'   it will be set automatically if \code{auto_prior=TRUE}, or to
#'   \code{half_student_t_prior(df=3, scale=2)} otherwise.
#' @param prior_beta A prior object for covariate effects \(\beta\). If \code{NULL},
#'   it will be set automatically if \code{auto_prior=TRUE}, or to
#'   \code{normal_prior(0, 10)} otherwise.
#' @param auto_prior Logical indicating whether to use data-dependent priors by default
#'   (baggr style). Defaults to \code{TRUE}.
#'
#' @return A named list ready to pass to Stan, with the elements:
#' \itemize{
#'   \item \code{J, tau_hat_j, se_j, Nc, X} (the data),
#'   \item \code{prior_hypermean_fam}, \code{prior_hypermean_val},
#'         \code{prior_hypersd_fam},   \code{prior_hypersd_val},
#'         \code{prior_beta_fam},      \code{prior_beta_val}.
#' }
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
#' # 1) Default usage (auto_prior=TRUE, all priors=NULL => auto priors for everything)
#' standata_auto <- build_standata_rubin_normal(
#'   df          = df_with_cov,
#'   tau_hat_col = "tau_j_hat",
#'   se2_col     = "se2_hat",
#'   covariates  = "X1"
#' )
#'
#' # 2) If you want to override just one of them, e.g. prior_hypersd:
#' my_prior_sd <- half_student_t_prior(df=4, scale=1.5)
#' standata_mixed <- build_standata_rubin_normal(
#'   df             = df_with_cov,
#'   tau_hat_col    = "tau_j_hat",
#'   se2_col        = "se2_hat",
#'   covariates     = "X1",
#'   prior_hypersd  = my_prior_sd  # user override for hyper-SD
#' )
#'
#' # 3) auto_prior=FALSE => anything left as NULL uses the simpler defaults (normal(0,10), etc.)
#' standata_manual <- build_standata_rubin_normal(
#'   df          = df_with_cov,
#'   tau_hat_col = "tau_j_hat",
#'   se2_col     = "se2_hat",
#'   covariates  = "X1",
#'   auto_prior  = FALSE
#'   # => prior_hypermean=NULL => normal(0,10)
#'   #    prior_hypersd=NULL   => half_student_t_prior(3,2)
#'   #    prior_beta=NULL      => normal(0,10)
#' )
#' }
#'
#' @export
build_standata_rubin_normal <- function(
    df,
    tau_hat_col     = "tau_j_hat",
    se2_col         = "se2_hat",
    covariates      = NULL,
    prior_hypermean = NULL,
    prior_hypersd   = NULL,
    prior_beta      = NULL,
    auto_prior      = TRUE
) {
  #### 1) If auto_prior = TRUE, apply partial or full auto logic
  if (isTRUE(auto_prior)) {
    # Call auto_prior_rubin_normal to get "auto" suggestions
    auto_priors <- auto_prior_rubin_normal(
      df          = df,
      tau_hat_col = tau_hat_col,
      se2_col     = se2_col,
      covariates  = covariates
    )

    # If user didn't supply a specific prior, use auto.
    # Otherwise, use the user-supplied version for that parameter.
    if (is.null(prior_hypermean)) {
      prior_hypermean <- auto_priors$hypermean
      message("Auto prior used for hypermean (mu).")
    }
    if (is.null(prior_hypersd)) {
      prior_hypersd <- auto_priors$hypersd
      message("Auto prior used for hypersd (tau).")
    }
    if (is.null(prior_beta)) {
      prior_beta <- auto_priors$beta
      message("Auto prior used for beta (covariates).")
    }

  } else {
    #### 2) If auto_prior = FALSE, revert to older style defaults:
    if (is.null(prior_hypermean)) {
      prior_hypermean <- normal_prior(0, 10)
      message("Using default normal_prior(0,10) for hypermean.")
    }
    if (is.null(prior_hypersd)) {
      prior_hypersd <- half_student_t_prior(df = 3, scale = 2)
      message("Using default half_student_t_prior(df=3, scale=2) for hypersd.")
    }
    if (is.null(prior_beta)) {
      prior_beta <- normal_prior(0, 10)
      message("Using default normal_prior(0,10) for beta.")
    }
  }

  #### 3) Prepare summary data (tau_j, se_j, covariates, etc.)
  sum_list <- prepare_summary_data_rubin(
    df           = df,
    tau_hat_col  = tau_hat_col,
    se2_col      = se2_col,
    covariates   = covariates
  )

  #### 4) Convert the final set of prior objects into Stan-compatible form
  prior_list <- prepare_prior_settings_rubin_normal(
    hypermean = prior_hypermean,
    hypersd   = prior_hypersd,
    beta      = prior_beta
  )

  #### 5) Return the final list
  list(
    # Data portion
    J            = sum_list$J,
    tau_hat_j    = sum_list$tau_hat_j,
    se_j         = sum_list$se_j,
    Nc           = sum_list$Nc,
    X            = sum_list$X,

    # Priors portion
    prior_hypermean_fam = prior_list$prior_hypermean_fam,
    prior_hypermean_val = prior_list$prior_hypermean_val,
    prior_hypersd_fam   = prior_list$prior_hypersd_fam,
    prior_hypersd_val   = prior_list$prior_hypersd_val,
    prior_beta_fam      = prior_list$prior_beta_fam,
    prior_beta_val      = prior_list$prior_beta_val
  )
}
