#' Simulate a complete multi-site dataset (no final n_j), with optional rank correlation
#'
#' This function does the following:
#' \enumerate{
#'   \item Calls \code{\link{gen_priorG2}} to generate \eqn{\tau_j} (true effects).
#'   \item Calls \code{\link{sim_sitesize_withinvar}} to generate \eqn{\text{se2}_j}.
#'         (Internally also generates \eqn{n_j}, but that is not included in the final output.)
#'   \item Calls \code{\link{sim_observed_effects}} to optionally reorder \eqn{\text{se2}_j}
#'         for partial Spearman correlation \eqn{\rho=\text{rank_corr}}, then simulates
#'         \(\hat{\tau}_j \sim N(\tau_j, se2_j)\).
#' }
#'
#' The user can specify \code{set_seed=TRUE} to do \code{set.seed(123)} for reproducible results.
#' By default, \code{set_seed=FALSE} means no seed is set. The final returned tibble has columns
#' \code{tau_j, se2_j, tau_j_hat, corr_est} in that order, omitting any site size column.
#'
#' @inheritParams gen_priorG2
#' @inheritParams sim_sitesize_withinvar
#' @inheritParams sim_observed_effects
#'
#' @param set_seed Logical (default \code{FALSE}). If \code{TRUE}, we do \code{set.seed(123)}
#'   at the start. If \code{FALSE}, no seed is set.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{\code{tau_j}}{The true site effects.}
#'   \item{\code{se2_j}}{The final within-site variances.}
#'   \item{\code{tau_j_hat}}{Observed site effects \(\sim N(\tau_j, se2_j)\).}
#'   \item{\code{corr_est}}{Spearman rank correlation between \code{tau_j} and \code{se2_j}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example:
#' df_sim <- sim_multisite_data(
#'   # For gen_priorG2:
#'   true_dist = "Gaussian", J = 10, variance = 1,
#'   # For sim_sitesize_withinvar:
#'   nj_mean = 20, cv = 0.5,
#'   # partial correlation ~ 0.4
#'   precision_dependence = TRUE, rank_corr = 0.4,
#'   set_seed=TRUE
#' )
#' df_sim
#'
#' cor(df_sim$tau_j, df_sim$se2_j, method="spearman") # near 0.4
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
  # 1) If set_seed=TRUE => set.seed(123)
  if (set_seed) {
    set.seed(123)
  }

  # 2) Generate true effects
  tau_j <- gen_priorG2(
    true_dist = true_dist,
    J         = J,
    sigma_tau = sigma_tau,
    tau       = tau,
    variance  = variance,
    nu        = nu,
    slant     = slant,
    rho       = rho,
    delta     = delta,
    eps       = eps,
    ups       = ups,
    formula   = formula,
    beta      = beta,
    data      = data,
    # don't re-set seed inside
    seed      = NULL,
    set_seed  = FALSE
  )

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
  se2_j <- df_site$se2_j  # ignoring n_j, alpha, beta

  # 4) Simulate observed effects with partial correlation if requested
  df_obs <- sim_observed_effects(
    tau_j    = tau_j,
    se2_j    = se2_j,
    precision_dependence = precision_dependence,
    rank_corr           = rank_corr,
    max_iter            = max_iter,
    tol                 = tol
  )
  # columns: tau_j, se2_j, tau_j_hat, n_j=NA, corr_est

  # 5) Return tibble with tau_j, se2_j, tau_j_hat, corr_est (ignore the n_j col from df_obs)
  dplyr::tibble(
    tau_j     = df_obs$tau_j,
    se2_j     = df_obs$se2_j,
    tau_j_hat = df_obs$tau_j_hat,
    corr_est  = df_obs$corr_est
  )
}
