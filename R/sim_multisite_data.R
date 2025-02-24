#' Simulate a multi-site dataset with possible rank correlation between tau_j and se2_j
#'
#' This function:
#' \enumerate{
#'   \item Calls \code{\link{gen_priorG2}} to generate \code{J} site-level "true" effects (including
#'         \code{site_index} and possibly covariates).
#'   \item Calls \code{\link{sim_sitesize_withinvar}} to generate \code{se2_j} for each site.
#'   \item Calls \code{\link{sim_observed_effects}}, which \emph{may} reorder \code{se2_j}
#'         internally to achieve the target Spearman rank correlation with \code{tau_j}.
#'         The function then merges these results back into the final data frame, leaving
#'         \code{tau_j} and covariates in their original row order.
#' }
#'
#' The returned data frame has columns:
#' \itemize{
#'   \item \strong{\code{site_index}} (unique ID for each site).
#'   \item \emph{covariates}, if a formula/data were used in \code{\link{gen_priorG2}}.
#'   \item \strong{\code{tau_j}}, the true site effect from \code{gen_priorG2}.
#'   \item \strong{\code{se2_j}}, the (possibly permuted) sampling variance.
#'   \item \strong{\code{tau_j_hat}}, the observed site effect, drawn from \(\mathrm{N}(\tau_j, se2_j)\).
#'   \item \strong{\code{corr_est}}, the final Spearman rank correlation between \code{tau_j}
#'         and \code{se2_j}.
#' }
#'
#' @inheritParams gen_priorG2
#' @inheritParams sim_sitesize_withinvar
#' @inheritParams sim_observed_effects
#'
#' @param set_seed Logical (default \code{FALSE}). If \code{TRUE}, the function begins by
#'   calling \code{set.seed(123)}, thus making the entire pipeline reproducible.
#'
#' @return A \strong{tibble} with \code{J} rows, containing site index, (optional) covariates,
#'   the true effects \code{tau_j}, the sampling variance \code{se2_j}, the observed
#'   \code{tau_j_hat}, and \code{corr_est}.
#'
#' @details
#' \subsection{Reordering of \eqn{se_j^2} for partial correlation}{
#'   If \code{precision_dependence = TRUE}, \code{\link{sim_observed_effects}} reassigns the
#'   \(\mathrm{se2}_j\) values among the sites (in an internal permutation) to achieve the
#'   specified \code{rank_corr} with \(\tau_j\).  \strong{However,} we do \emph{not} reorder
#'   \(\tau_j\) or any covariates. That means, effectively, site \(i\) keeps its row and
#'   \(\tau_j[i]\) but might receive a different \(\mathrm{se2}_j\). This allows us to produce
#'   a rank correlation between the unaltered \(\tau_j\) and a re-labeled vector of \(\mathrm{se2}_j\),
#'   while leaving covariates matched to the same row index.
#' }
#'
#' @examples
#' \dontrun{
#' ##################################################
#' # Example 1: No covariates, partial correlation
#' ##################################################
#' df_no_cov <- sim_multisite_data(
#'   true_dist  = "Gaussian",
#'   J          = 10,
#'   sigma_tau  = 0.25,
#'   variance   = 1,
#'   # site size generator:
#'   nj_mean    = 20,
#'   cv         = 0.3,
#'   # enforce correlation ~ 0.4
#'   precision_dependence = TRUE,
#'   rank_corr  = 0.4,
#'   set_seed   = TRUE
#' )
#' head(df_no_cov)
#'
#' cor(df_no_cov$tau_j, df_no_cov$se2_j, method="spearman")  # near 0.4
#'
#' ##################################################
#' # Example 2: Covariates, binomial random X1,X2,X3
#' ##################################################
#' J <- 100
#' p1 <- 0.5
#' p2 <- 0.3
#' p3 <- 0.7
#'
#' X1 <- rbinom(J, size=1, prob=p1)
#' X2 <- rbinom(J, size=1, prob=p2)
#' X3 <- rbinom(J, size=1, prob=p3)
#'
#' site_data <- data.frame(X1 = X1, X2 = X2, X3 = X3)
#'
#' # formula => ~ X1 + X2 + X3 includes an intercept => total 4 columns
#' beta_vec <- c(2, 1.5, -1, 0.5)
#'
#' df_with_cov <- sim_multisite_data(
#'   true_dist = "Gaussian",
#'   J         = J,
#'   sigma_tau = 0.25,
#'   variance  = 1,
#'   formula   = ~ X1 + X2 + X3,
#'   beta      = beta_vec,
#'   data      = site_data,
#'
#'   nj_mean   = 20,
#'   cv        = 0.5,
#'
#'   precision_dependence = TRUE,
#'   rank_corr = 0.3,
#'
#'   set_seed = TRUE
#' )
#' head(df_with_cov)
#'
#' cor(df_with_cov$tau_j, df_with_cov$se2_j, method="spearman")  # near 0.3
#' }
#'
#' @export
sim_multisite_data <- function(
    # from gen_priorG2
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

  # from sim_sitesize_withinvar
  nj_mean = 10,
  cv = 0.3,
  nj_min = 5,
  p = 0.5,
  R2 = 0.0,

  # from sim_observed_effects
  precision_dependence = FALSE,
  rank_corr = 0.5,
  max_iter = 20000,
  tol = 0.01,

  set_seed = FALSE
) {
  # 1) Possibly set seed
  if (set_seed) {
    set.seed(123)
  }

  # 2) Generate "true" site effects => includes site_index, covariates (if any),
  #    and a column named 'site_effect' which we'll rename to 'tau_j'.
  df_tau <- gen_priorG2(
    true_dist  = true_dist,
    J          = J,
    sigma_tau  = sigma_tau,
    tau        = tau,
    variance   = variance,
    nu         = nu,
    slant      = slant,
    rho        = rho,
    delta      = delta,
    eps        = eps,
    ups        = ups,
    formula    = formula,
    beta       = beta,
    data       = data,
    seed       = NULL,       # do not reset inside
    set_seed   = FALSE
  )
  if (!"site_effect" %in% names(df_tau)) {
    stop("gen_priorG2 should return 'site_effect'; not found.")
  }
  df_tau <- dplyr::rename(df_tau, tau_j = site_effect)

  # 3) Generate site sizes & within-site variances
  df_site <- sim_sitesize_withinvar(
    J       = J,
    nj_mean = nj_mean,
    cv      = cv,
    nj_min  = nj_min,
    p       = p,
    R2      = R2,
    seed    = NULL,
    set_seed= FALSE
  )
  se2_j <- df_site$se2_j

  # 4) Possibly reorder se2_j to produce rank_corr w.r.t. tau_j, then simulate observed tau_j_hat
  df_obs <- sim_observed_effects(
    tau_j               = df_tau$tau_j,
    se2_j               = se2_j,
    precision_dependence= precision_dependence,
    rank_corr           = rank_corr,
    max_iter            = max_iter,
    tol                 = tol,
    seed                = NULL
  )
  # => columns: tau_j, se2_j, tau_j_hat, corr_est
  # The function does NOT reorder df_tau, only reassigns se2_j among sites.

  # 5) Combine results
  #    Because sim_observed_effects() keeps tau_j in the same row i, we can just attach new columns:
  #    'se2_j', 'tau_j_hat', 'corr_est' to df_tau in the same row i.
  df_tau$se2_j     <- df_obs$se2_j
  df_tau$tau_j_hat <- df_obs$tau_j_hat
  # The rank correlation is the same for all rows, so store df_obs$corr_est[1].
  df_tau$corr_est  <- df_obs$corr_est[1]

  # 6) Return as tibble
  dplyr::as_tibble(df_tau)
}
