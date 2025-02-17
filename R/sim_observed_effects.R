#' @keywords internal
reorder_for_spearman <- function(
    tau_j,
    se2_j,
    target = 0.5,
    max_iter = 20000,
    tol = 0.01,
    seed = NULL
) {
  # Internal function: iterative swapping to reorder `se2_j` so that
  # Spearman rank correlation with `tau_j` is near `target`.
  #
  # Hill-climbing approach:
  #   1) Start with a random permutation of se2_j
  #   2) Randomly pick two indices, swap them if it brings correlation
  #      closer to `target`
  #   3) Stop if we reach `max_iter` or if within `tol` of `target`

  if (!is.null(seed)) set.seed(seed)

  J <- length(tau_j)
  if (length(se2_j) != J) {
    stop("Lengths of tau_j and se2_j must match.")
  }
  if (target < -1 || target > 1) {
    stop("Target must be between -1 and 1.")
  }

  # Start with random permutation
  best_perm <- sample.int(J)
  best_se2  <- se2_j[best_perm]
  best_corr <- suppressWarnings(cor(tau_j, best_se2, method="spearman"))

  for (iter in seq_len(max_iter)) {
    swap_idx <- sample.int(J, 2)
    new_perm <- best_perm

    # swap
    tmp <- new_perm[swap_idx[1]]
    new_perm[swap_idx[1]] <- new_perm[swap_idx[2]]
    new_perm[swap_idx[2]] <- tmp

    test_se2 <- se2_j[new_perm]
    test_corr <- suppressWarnings(cor(tau_j, test_se2, method="spearman"))

    if (abs(test_corr - target) < abs(best_corr - target)) {
      best_corr <- test_corr
      best_perm <- new_perm
      best_se2  <- test_se2

      if (abs(best_corr - target) < tol) {
        break
      }
    }
  }
  best_se2
}

#' Simulate observed site-specific effects with optional rank-based partial correlation
#'
#' This function takes user-supplied vectors of \emph{true} site effects (\eqn{\tau_j})
#' and site-level sampling variances (\eqn{\text{se2}_j}), then optionally reorders
#' \eqn{\text{se2}_j} to achieve a given Spearman rank correlation with \eqn{\tau_j}.
#' Finally, it simulates observed estimates \eqn{\hat{\tau}_j} from a normal distribution:
#' \deqn{\hat{\tau}_j \sim N(\tau_j, \text{se2}_j).}
#'
#' By default, \code{precision_dependence = FALSE} means we forcibly random-shuffle
#' \eqn{\text{se2}_j} so that its Spearman correlation with \eqn{\tau_j} is near zero.
#' If \code{precision_dependence = TRUE}, an \strong{iterative swap approach} is used
#' to reorder \eqn{\text{se2}_j} until its Spearman correlation with \eqn{\tau_j} is
#' near \code{rank_corr} (default \code{0.50}). This procedure preserves the marginal
#' distributions of both vectors but does not guarantee an exact match to \code{rank_corr}
#' or a global optimum. It typically suffices for simulation contexts where a moderate
#' rank correlation is desired without distorting the underlying values.
#'
#' @param tau_j Numeric vector. True site effects of length \eqn{J}, e.g. from \code{gen_priorG2}.
#' @param se2_j Numeric vector of the same length as \code{tau_j}, e.g. from
#'  \code{sim_sitesize_withinvar}, specifying the within-site sampling variance
#'  \eqn{\text{se2}_j}.
#' @param precision_dependence Logical indicating whether to enforce a partial correlation
#'  (default \code{FALSE}). If \code{FALSE}, we random-shuffle \code{se2_j} to break
#'  any existing correlation (near-zero result).
#'  If \code{TRUE}, we reorder \code{se2_j} to achieve \code{rank_corr}.
#' @param rank_corr Numeric in \eqn{[-1,1]}, the target Spearman rank correlation to achieve
#'  between \eqn{\tau_j} and \eqn{\text{se2}_j} (default 0.50). Only used if
#'  \code{precision_dependence = TRUE}.
#' @param max_iter Integer. Maximum number of swap attempts for the iterative
#'  correlation injection (default 20000).
#' @param tol Numeric tolerance for stopping early if the current correlation is
#'  within \code{tol} of \code{rank_corr} (default 0.01).
#' @param seed Optional integer for reproducibility (default \code{NULL}). If not \code{NULL},
#'  it is used to set the random seed before reordering (and before simulating observed effects).
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{\code{tau_j}}{The original true site effects (unchanged).}
#'   \item{\code{se2_j}}{The within-site variances after correlation injection
#'                       (\emph{or random shuffle} if \code{precision_dependence=FALSE}).}
#'   \item{\code{tau_j_hat}}{Simulated observed effects, drawn from \eqn{N(\tau_j, se2_j)}.}
#'   \item{\code{corr_est}}{The final Spearman rank correlation between \eqn{\tau_j} and
#'                          the new \eqn{se2_j}.}
#' }
#'
#' @details
#' \subsection{Iterative Swap Approach}{
#' If \code{precision_dependence = TRUE}, we reorder \eqn{\text{se2}_j} using a
#' hill-climbing swap algorithm (see \code{\link[=reorder_for_spearman]{reorder_for_spearman()}} internally).
#' We initialize from a random permutation, then repeatedly pick two indices at random
#' and swap them if it brings the Spearman correlation closer to \code{rank_corr}.
#' This preserves the original set of values in \eqn{\text{se2}_j}, guaranteeing
#' no marginal distortion, but does not ensure a global optimum. Typically,
#' it yields a final correlation within \code{tol} of \code{rank_corr}.
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose you have 50 sites
#' set.seed(123)
#' tau_vec <- rnorm(50, 0, 1)
#' se2_vec <- rexp(50, rate=1) + 0.1
#'
#' # 1) precision_dependence=FALSE => forcibly random-shuffle se2_j => near-zero correlation
#' df_no_corr <- sim_observed_effects(
#'   tau_j = tau_vec,
#'   se2_j = se2_vec,
#'   precision_dependence = FALSE
#' )
#' df_no_corr$corr_est
#'
#' # 2) precision_dependence=TRUE => iterative approach => ~0.50 rank correlation
#' df_corr <- sim_observed_effects(
#'   tau_j = tau_vec,
#'   se2_j = se2_vec,
#'   precision_dependence = TRUE,
#'   rank_corr = 0.5,
#'   max_iter = 20000,
#'   tol = 0.01,
#'   seed = 999
#' )
#' df_corr$corr_est
#' }
#'
#' @export
sim_observed_effects <- function(
    tau_j,
    se2_j,
    precision_dependence = FALSE,
    rank_corr = 0.5,
    max_iter = 20000,
    tol = 0.01,
    seed = NULL
) {
  J <- length(tau_j)
  if (length(se2_j) != J) {
    stop("tau_j and se2_j must have the same length.")
  }

  if (!is.null(seed)) set.seed(seed)

  if (!precision_dependence) {
    # forcibly random-shuffle se2_j to near-zero correlation
    se2_new <- sample(se2_j)
  } else {
    # iterative approach to get rank_corr
    se2_new <- reorder_for_spearman(
      tau_j = tau_j,
      se2_j = se2_j,
      target = rank_corr,
      max_iter = max_iter,
      tol = tol,
      seed = NULL  # already set globally
    )
  }

  # final correlation
  corr_est <- suppressWarnings(cor(tau_j, se2_new, method="spearman"))
  # observed estimate
  tau_j_hat <- rnorm(J, mean = tau_j, sd = sqrt(se2_new))

  dplyr::tibble(
    tau_j     = tau_j,
    se2_j     = se2_new,
    tau_j_hat = tau_j_hat,
    corr_est  = corr_est
  )
}
