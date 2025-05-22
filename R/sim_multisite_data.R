#' Simulate a multi-site dataset, optionally with correlation between tau_j and se2_j
#'
#' This function:
#' \enumerate{
#'   \item Calls \code{\link{gen_priorG2}} to generate \code{J} site-level "true" effects (including
#'         \code{site_index} and possibly covariates).
#'   \item Calls \code{\link{sim_sitesize_withinvar}} to generate \code{se2_j} for each site.
#'   \item Depending on user input, calls either \code{\link{sim_observed_effects}} (for Spearman rank correlation)
#'         or \code{\link{sim_observed_effects_copula}} (for Pearson correlation) to potentially reorder \code{se2_j}.
#'   \item The function merges these results back into the final data frame, leaving
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
#'   \item \strong{\code{corr_est}}, the final correlation estimate (Spearman or Pearson) between \code{tau_j}
#'         and \code{se2_j}, depending on which method was used.
#' }
#'
#' @inheritParams gen_priorG2
#' @inheritParams sim_sitesize_withinvar
#'
#' @param precision_dependence Logical (default \code{FALSE}). If \code{TRUE}, will attempt
#'   to induce correlation between \(\tau_j\) and \(\mathrm{se2}_j\). The type of correlation
#'   depends on whether \code{rank_corr} or \code{pearson_corr} is specified (see below).
#' @param rank_corr Numeric in [-1,1], optional. If provided (and \code{precision_dependence=TRUE}),
#'   we enforce a \emph{Spearman} correlation target using \code{\link{sim_observed_effects}}
#'   (hill-climbing approach).
#' @param pearson_corr Numeric in [-1,1], optional. If provided (and \code{precision_dependence=TRUE}),
#'   we enforce a \emph{Pearson} correlation target using \code{\link{sim_observed_effects_copula}}
#'   (Gaussian copula approach).
#'
#' @param max_iter Integer. Maximum number of swap attempts for the hill-climbing approach
#'   (only relevant if \code{rank_corr} is used).
#' @param tol Numeric tolerance for stopping early if the current Spearman correlation is
#'   within \code{tol} of \code{rank_corr} (only relevant if \code{rank_corr} is used).
#'
#' @param set_seed Logical (default \code{FALSE}). If \code{TRUE}, the function begins by
#'   calling \code{set.seed(123)}, thus making the entire pipeline reproducible.
#'
#' @return A \strong{tibble} with \code{J} rows, containing site index, (optional) covariates,
#'   the true effects \code{tau_j}, the sampling variance \code{se2_j}, the observed
#'   \code{tau_j_hat}, and \code{corr_est}.
#'
#' @details
#' \subsection{Choosing correlation type}{
#'   - If \code{precision_dependence=FALSE}, no correlation is enforced (we do a random shuffle
#'     so it's near zero).
#'   - If \code{precision_dependence=TRUE} and \code{rank_corr} is specified (not \code{NULL}),
#'     we use \code{sim_observed_effects} to achieve a Spearman correlation (hill climbing).
#'   - If \code{precision_dependence=TRUE} and \code{pearson_corr} is specified (not \code{NULL}),
#'     we use \code{sim_observed_effects_copula} to achieve a Pearson correlation (Gaussian copula).
#'   - If both are specified or neither is specified, we throw an error or interpret it as zero
#'     correlation. You can customize this logic as you prefer.
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: rank_corr => Spearman correlation ~ 0.4
#' df_sp <- sim_multisite_data(
#'   J = 20,
#'   precision_dependence = TRUE,
#'   rank_corr = 0.4,
#'   set_seed = TRUE
#' )
#' cor(df_sp$tau_j, df_sp$se2_j, method="spearman") # ~0.4
#'
#' # Example 2: pearson_corr => Pearson correlation ~ 0.5
#' df_pe <- sim_multisite_data(
#'   J = 20,
#'   precision_dependence = TRUE,
#'   pearson_corr = 0.5,
#'   set_seed = TRUE
#' )
#' cor(df_pe$tau_j, df_pe$se2_j) # ~0.5
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

  # correlation logic
  precision_dependence = FALSE,
  rank_corr = NULL,
  pearson_corr = NULL,
  max_iter = 20000,
  tol = 0.01,

  set_seed = FALSE
) {
  # 1) Possibly set seed globally
  if (set_seed) {
    set.seed(123)
  }

  # 2) Basic validation: only one correlation type allowed
  if (precision_dependence) {
    # both not null => conflict
    if (!is.null(rank_corr) && !is.null(pearson_corr)) {
      stop("Cannot specify both rank_corr and pearson_corr. Choose one.")
    }
  }

  # 3) Generate "true" site effects => includes site_index, covariates (if any),
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
    seed       = NULL,
    set_seed   = FALSE
  )
  if (!"site_effect" %in% names(df_tau)) {
    stop("gen_priorG2 should return 'site_effect'; not found.")
  }
  df_tau <- dplyr::rename(df_tau, tau_j = site_effect)

  # 4) Generate site sizes & within-site variances
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

  # 5) Build final observed effects with correlation (if requested)
  if (!precision_dependence) {
    # no correlation => near zero
    # => we can simply call sim_observed_effects with rank_corr=0 but not do the hill climb
    #    or call sim_observed_effects_copula with pearson_corr=0, or just random shuffle
    #    For simplicity, let's use sim_observed_effects(...) with precision_dependence=FALSE
    #    => random shuffle => near zero
    df_obs <- sim_observed_effects(
      tau_j   = df_tau$tau_j,
      se2_j   = se2_j,
      precision_dependence = FALSE
    )
  } else {
    # precision_dependence = TRUE => pick which method
    if (!is.null(rank_corr)) {
      # Spearman => sim_observed_effects
      df_obs <- sim_observed_effects(
        tau_j   = df_tau$tau_j,
        se2_j   = se2_j,
        precision_dependence = TRUE,
        rank_corr           = rank_corr,
        max_iter            = max_iter,
        tol                 = tol
      )
    } else if (!is.null(pearson_corr)) {
      # Pearson => sim_observed_effects_copula
      df_obs <- sim_observed_effects_copula(
        tau_j   = df_tau$tau_j,
        se2_j   = se2_j,
        precision_dependence = TRUE,
        pearson_corr         = pearson_corr
      )
    } else {
      stop("precision_dependence=TRUE but neither rank_corr nor pearson_corr was specified.")
    }
  }

  # df_obs => columns: tau_j, se2_j, tau_j_hat, corr_est
  # The function does NOT reorder df_tau, only reassigns se2_j among sites.

  # 6) Combine results
  #    Because sim_observed_effects* keep tau_j in the same row i, we can just attach
  #    new columns: 'se2_j', 'tau_j_hat', 'corr_est' to df_tau in the same row i.
  df_tau$se2_j     <- df_obs$se2_j
  df_tau$tau_j_hat <- df_obs$tau_j_hat
  # The correlation is the same for all rows, so store df_obs$corr_est[1].
  df_tau$corr_est  <- df_obs$corr_est[1]

  # 7) Return as tibble
  dplyr::as_tibble(df_tau)
}
