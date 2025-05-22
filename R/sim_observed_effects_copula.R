#' Copula-based approach to simulate observed effects with optional Pearson correlation
#'
#' This function takes user-supplied vectors of \emph{true} site effects (\(\tau_j\))
#' and site-level sampling variances (\(\mathrm{se2}_j\)), then optionally reorders
#' or re-generates \(\mathrm{se2}_j\) via a Gaussian copula approach so that
#' the final (\(\tau_j\), \(\mathrm{se2}_j\)) has a specified Pearson correlation
#' (\code{pearson_corr}) in the \emph{standard-normal (copula)} transform space.
#'
#' @param tau_j Numeric vector of length J (the "true" site effects).
#' @param se2_j Numeric vector of length J (the within-site sampling variances).
#' @param precision_dependence Logical: if \code{TRUE}, enforce a correlation structure
#'   via a Gaussian copula. If \code{FALSE}, we simply random-shuffle \code{se2_j},
#'   aiming for near-zero correlation.
#' @param pearson_corr Numeric in [-1,1], the desired correlation in the \emph{Gaussian copula}
#'   space. Defaults to 0.5. Note that for heavily skewed or discrete-like data,
#'   the actual Pearson correlation \code{cor(tau_j, se2_j)} in the final dataset
#'   may not be \emph{exactly} \code{pearson_corr}, but typically it will be close.
#' @param seed Optional integer for reproducibility. If not \code{NULL}, \code{set.seed(seed)} is called
#'   at the start of this function.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{\code{tau_j}}{The original \(\tau_j\), unchanged.}
#'   \item{\code{se2_j}}{The updated \(\mathrm{se2}_j\), either random-shuffled
#'         (if \code{precision_dependence=FALSE}) or derived via copula
#'         (if \code{precision_dependence=TRUE}). The marginal distribution
#'         remains the same as the original \(\mathrm{se2}_j\) (empirical).}
#'   \item{\code{tau_j_hat}}{Simulated observed effects, drawn from
#'         \(\mathcal{N}(\tau_j, \mathrm{se2}_j)\).}
#'   \item{\code{corr_est}}{The final Pearson correlation \code{cor(tau_j, se2_j)}
#'         in the output dataset.}
#' }
#'
#' @details
#' \itemize{
#'   \item If \code{precision_dependence = FALSE}, we shuffle \code{se2_j} randomly
#'         (via \code{sample()}) to break most existing correlation. The resulting
#'         correlation is typically near zero, but not guaranteed to be exactly zero.
#'   \item If \code{precision_dependence = TRUE}, we apply the standard copula logic:
#'         (1) convert \(\tau_j\) and \(\mathrm{se2}_j\) to normal \code{z}-scores via
#'         rank \(\to\) percentile \(\to\) \code{qnorm}; (2) impose a linear correlation
#'         \code{pearson_corr} in \code{z}-space; (3) map back to the empirical distribution
#'         of \(\mathrm{se2}_j\). This preserves the marginal distribution of each vector
#'         while creating an approximate Pearson correlation in the final dataset.
#'   \item Ties are broken randomly (\code{ties.method="random"}) for rank calculations.
#'         For data with many identical values, results might vary slightly.
#' }
#'
#' @examples
#' \dontrun{
#'   set.seed(101)
#'   J <- 50
#'   tau_vec <- rnorm(J, 2, 1)
#'   se2_vec <- rgamma(J, shape=2, rate=4)
#'
#'   # 1) No correlation (just random shuffle)
#'   df_nocorr <- sim_observed_effects_copula(tau_vec, se2_vec,
#'       precision_dependence = FALSE)
#'   cor(df_nocorr$tau_j, df_nocorr$se2_j)  # near zero (not exact)
#'
#'   # 2) Copula approach => target ~0.5 correlation in normal space
#'   df_corr <- sim_observed_effects_copula(tau_vec, se2_vec,
#'       precision_dependence = TRUE, pearson_corr = 0.5)
#'   cor(df_corr$tau_j, df_corr$se2_j)      # check actual correlation
#' }
#'
#' @export
sim_observed_effects_copula <- function(
    tau_j,
    se2_j,
    precision_dependence = FALSE,
    pearson_corr = 0.5,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  J <- length(tau_j)
  if (length(se2_j) != J) {
    stop("tau_j and se2_j must have the same length.")
  }

  #------------------------------------------------------
  # case 1: precision_dependence=FALSE => random shuffle
  #------------------------------------------------------
  if (!precision_dependence) {
    se2_new <- sample(se2_j)

    #------------------------------------------------------
    # case 2: precision_dependence=TRUE => Gaussian copula
    #------------------------------------------------------
  } else {
    # 1) Convert each vector to "z-space" via rank -> percentile -> qnorm
    convert_to_z <- function(x) {
      r <- rank(x, ties.method="random")
      # percentile in (0,1)
      u <- (r - 0.5) / length(x)
      z <- qnorm(u)  # standard normal
      z
    }
    z_tau  <- convert_to_z(tau_j)
    z_se2  <- convert_to_z(se2_j)

    # 2) We want new z_se2_new that has correlation = pearson_corr with z_tau
    rho <- pearson_corr
    if (abs(rho) > 1) {
      stop("pearson_corr must be in [-1,1].")
    }

    # z_se2_new[i] = rho * z_tau[i] + sqrt(1-rho^2)* eps_i
    eps <- rnorm(J, mean=0, sd=1)
    z_se2_new <- rho * z_tau + sqrt(1 - rho^2) * eps

    # 3) Convert z_se2_new back to percentile => map to the empirical distribution of se2_j
    r_se2_new <- pnorm(z_se2_new)
    empDist <- sort(se2_j)

    idx <- ceiling(r_se2_new * J)
    idx[idx < 1]  <- 1
    idx[idx > J]  <- J
    se2_new <- empDist[idx]
  }

  # final correlation (Pearson) in the generated data
  corr_est <- suppressWarnings(cor(tau_j, se2_new))

  # generate observed tau_j_hat
  tau_j_hat <- rnorm(J, mean = tau_j, sd = sqrt(se2_new))

  dplyr::tibble(
    tau_j     = tau_j,
    se2_j     = se2_new,
    tau_j_hat = tau_j_hat,
    corr_est  = corr_est
  )
}
