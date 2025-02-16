#' @keywords internal
#' @importFrom stats integrate pgamma dgamma
trunc_gamma_moments <- function(alpha, beta, m) {
  # For Gamma(shape=alpha, rate=beta), truncated below m => x >= m
  # We want E[x|x>=m], Var[x|x>=m].

  S <- 1 - pgamma(m, shape=alpha, rate=beta)
  if (S < 1e-15) {
    return(c(NA, NA))  # numerical issue if almost all mass < m
  }

  # E[x|x>=m] = (1/S) * ∫(m to ∞) x * dgamma(x;alpha,beta) dx
  int_mean <- integrate(
    f = function(x) x * dgamma(x, shape=alpha, rate=beta),
    lower = m, upper = Inf, rel.tol = 1e-8
  )
  if (int_mean$message != "OK") return(c(NA, NA))
  Ex_trunc <- int_mean$value / S

  # E[x^2|x>=m]
  int_m2 <- integrate(
    f = function(x) x^2 * dgamma(x, shape=alpha, rate=beta),
    lower = m, upper = Inf, rel.tol = 1e-8
  )
  if (int_m2$message != "OK") return(c(NA, NA))
  Ex2_trunc <- int_m2$value / S

  Var_trunc <- Ex2_trunc - (Ex_trunc)^2
  c(Ex_trunc, Var_trunc)
}

#' @keywords internal
#' @importFrom nleqslv nleqslv
solve_trunc_gamma <- function(m, targetMean, targetVar, init_alpha=2, init_rate=1) {

  f_system <- function(par) {
    alpha <- par[1]
    beta  <- par[2]

    # invalid => large penalty
    if (alpha <= 0 || beta <= 0) {
      return(c(1e10, 1e10))
    }

    mm <- trunc_gamma_moments(alpha, beta, m)
    if (any(is.na(mm))) {
      return(c(1e10, 1e10))
    }
    # We want: mean=targetMean, var=targetVar
    c(mm[1] - targetMean,
      mm[2] - targetVar)
  }

  sol <- nleqslv::nleqslv(x = c(init_alpha, init_rate), fn = f_system)

  if (sol$termcd != 1) {
    warning("nleqslv did not converge. termcd=", sol$termcd)
  }
  list(alpha = sol$x[1], beta = sol$x[2], convergence = sol$termcd, details = sol)
}

#' @keywords internal
#' @importFrom stats rgamma
rtrunc_gamma <- function(J, alpha, beta, m, seed=NULL) {
  # Accept-reject sampler for truncated Gamma
  # shape=alpha, rate=beta, truncated below m => X >= m

  if (!is.null(seed)) {
    set.seed(seed)
  }

  out <- numeric(J)
  n_acc <- 0
  while(n_acc < J) {
    # oversample *2 each time
    x_candidate <- rgamma(n = (J - n_acc)*2, shape=alpha, rate=beta)
    x_keep <- x_candidate[x_candidate >= m]
    n_take <- min(length(x_keep), J - n_acc)

    if (n_take > 0) {
      out[(n_acc+1):(n_acc+n_take)] <- x_keep[1:n_take]
      n_acc <- n_acc + n_take
    }
  }
  out
}

#' Generate Site Sizes and Within-Site Variances
#'
#' This function (internally using a truncated Gamma distribution) generates site sizes
#' (\code{n_j}) with a specified mean (\code{nj_mean}) and coefficient of variation
#' (\code{cv}), imposing a lower bound \code{nj_min}. It also computes the **within-site
#' sampling variances** (\code{se2_j}) under a Neyman-type design-based perspective, where
#' the outcome is expressed in *effect size units* (i.e., \eqn{\mathrm{Var}(Y)=1}).
#'
#' @section Connecting sampling errors to site sizes:
#' In effect-size units (\eqn{\mathrm{Var}(Y)=1}), one can approximate the sampling variance
#' of an estimated site-specific effect \eqn{\hat{\tau}_j} as:
#' \deqn{
#'   (\hat{se}_j)^2 = \frac{1}{n_j} \Bigl(\frac{1}{p_j} + \frac{1}{1 - p_j}\Bigr)\,(1 - R^2),
#' }
#' where \eqn{n_j} is the site size, \eqn{p_j} is the proportion assigned to treatment
#' (the propensity score; Rosenbaum & Rubin, 1983), and \eqn{R^2} is the explanatory power
#' of level-1 covariates. This formula follows a Neyman repeated-sampling approach under
#' homoskedasticity (see Miratrix, Weiss, & Henderson, 2021). If we further assume a
#' *constant* proportion \eqn{p} across sites, then
#' \deqn{
#'   (\hat{se}_j)^2 = \frac{1}{n_j}\Bigl(\frac{1}{p} + \frac{1}{1 - p}\Bigr)\,(1 - R^2)
#'   \;=\;\frac{\kappa}{n_j}.
#' }
#' By default, this function sets \eqn{p=0.5} and \eqn{R^2=0}, hence \eqn{\kappa=4}, so
#' that \eqn{(\hat{se}_j)^2=4/n_j}.
#'
#' @param J Integer > 0. Number of sites.
#' @param nj_mean Positive numeric. Desired truncated-gamma mean of site sizes.
#' @param cv Numeric >= 0. Coefficient of variation for the truncated Gamma.
#' @param nj_min Positive numeric. Lower bound for truncation (default 5).
#' @param p Numeric in (0,1). Proportion assigned to treatment (defaults to 0.5).
#' @param R2 Numeric in [0,1). Proportion of variance explained by level-1 covariates
#'           (defaults to 0.0).
#' @param seed Optional integer for random seed (default 123).
#' @param set_seed Logical. If TRUE and \code{seed} != NULL, sets the global seed (default FALSE).
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{\code{alpha}}{Shape parameter for truncated Gamma (found via numeric solver).}
#'   \item{\code{beta}}{Rate parameter for truncated Gamma (found via numeric solver).}
#'   \item{\code{n_j_raw}}{Continuous site sizes drawn from truncated Gamma(\code{alpha}, \code{beta}), \eqn{\ge nj_min}.}
#'   \item{\code{n_j}}{Integer-rounded site sizes.}
#'   \item{\code{se2_j}}{Within-site sampling variance, computed as \eqn{\kappa/n_j}, based on the
#'                      \eqn{\Bigl(\tfrac{1}{p} + \tfrac{1}{1-p}\Bigr) (1 - R^2)} relationship.}
#' }
#'
#' @details
#' \subsection{Truncated Gamma Approach}{
#'   Internally, the function solves for \eqn{\alpha,\beta} so that:
#'   \deqn{
#'     E[X \mid X \ge nj_{\min}] = nj_{\mathrm{mean}}, \quad
#'     Var\bigl(X \mid X \ge nj_{\min}\bigr] = (cv \times nj_{\mathrm{mean}})^2.
#'   }
#'   We then sample from \eqn{\Gamma(\alpha,\beta)} but \emph{reject} any draws below
#'   \eqn{nj_{\min}}. This is known as accept-reject sampling for a truncated distribution.
#'   The final \eqn{n_j} values are \emph{rounded} to integers.
#' }
#'
#' @references
#' \itemize{
#'   \item Rosenbaum, P. R., & Rubin, D. B. (1983). The central role of the propensity score
#'         in observational studies for causal effects. \emph{Biometrika}, 70(1), 41–55.
#'   \item Miratrix, L. W., Weiss, M. J., & Henderson, B. (2021). An applied researcher's
#'         guide to estimating effects from multisite individually randomized trials:
#'         Estimands, estimators, and estimates. \emph{Journal of Research on Educational
#'         Effectiveness}, 14(1), 270–308.
#' }
#'
#' @importFrom dplyr tibble
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' df_test <- sim_sitesize_withinvar(
#'   J = 1000,
#'   nj_mean = 10,
#'   cv = 0.5,
#'   nj_min = 5,
#'   p = 0.5,    # => kappa=4 if R2=0
#'   R2 = 0.0,
#'   seed = 999,
#'   set_seed = TRUE
#' )
#' head(df_test)
#' }
#'
#' @export
sim_sitesize_withinvar <- function(
    J = 25,
    nj_mean = 10,
    cv = 0.3,
    nj_min = 5,
    p = 0.5,
    R2 = 0.0,
    seed = 123,
    set_seed = FALSE
) {
  # ---- Error checks ----
  if (!is.numeric(J) || J <= 0 || floor(J) != J) {
    stop("`J` must be a positive integer.")
  }
  if (!is.numeric(nj_mean) || nj_mean <= 0) {
    stop("`nj_mean` must be a positive number.")
  }
  if (!is.numeric(cv) || cv < 0) {
    stop("`cv` must be >= 0.")
  }
  if (!is.numeric(nj_min) || nj_min <= 0) {
    stop("`nj_min` must be a positive number.")
  }
  if (!is.numeric(p) || p <= 0 || p >= 1) {
    stop("`p` must be in (0,1).")
  }
  if (!is.numeric(R2) || R2 < 0 || R2 >= 1) {
    stop("`R2` must be in [0,1).")
  }

  # Possibly set global seed if requested
  if (set_seed && !is.null(seed)) {
    set.seed(seed)
  }

  # 1) Solve truncated gamma shape=alpha, rate=beta
  targetMean <- nj_mean
  targetVar  <- (cv * nj_mean)^2
  sol <- solve_trunc_gamma(
    m           = nj_min,
    targetMean  = targetMean,
    targetVar   = targetVar,
    init_alpha  = 2,
    init_rate   = 0.5
  )
  alpha_hat <- sol$alpha
  beta_hat  <- sol$beta

  # 2) Sample from truncated gamma (no re-seed in rtrunc_gamma)
  n_j_raw <- rtrunc_gamma(
    J     = J,
    alpha = alpha_hat,
    beta  = beta_hat,
    m     = nj_min,
    seed  = NULL
  )

  # 3) Round to integer for final site sizes
  n_j <- round(n_j_raw)

  # 4) Compute sampling variance se2_j
  #    varY=1 => kappa = (1/p + 1/(1-p))*(1 - R2)
  varY <- 1
  kappa <- varY * (1/p + 1/(1 - p)) * (1 - R2)
  se2_j <- kappa * (1 / n_j)

  # 5) Build final tibble
  df <- dplyr::tibble(
    alpha   = alpha_hat,
    beta    = beta_hat,
    n_j_raw = n_j_raw,
    n_j     = n_j,
    se2_j   = se2_j
  )

  return(df)
}

