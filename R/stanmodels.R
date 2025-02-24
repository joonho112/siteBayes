# Create an environment to store compiled Stan models
stanmodels <- new.env(parent = emptyenv())

#' Compile the rubin_normal Stan model (internal)
#'
#' This function checks if the rubin_normal Stan model is already compiled
#' (stored in the \code{stanmodels} environment). If not, it compiles
#' \code{rubin_normal.stan} using \code{rstan::stan_model()} and caches
#' the resulting model in \code{stanmodels$rubin_normal}.
#'
#' @details
#' The \code{rubin_normal.stan} file is assumed to reside in the
#' \code{inst/stan/} folder of this package. The compiled model is stored
#' in \code{stanmodels$rubin_normal}, so subsequent calls will not re-compile.
#'
#' This function is primarily called internally by \code{fit_rubin_normal()}.
#' Most users will not need to call it directly.
#'
#' @return A \code{stanmodel} object corresponding to \code{rubin_normal.stan}.
#' @keywords internal
compile_rubin_normal_mod <- function() {
  # If the model is already in stanmodels, reuse it
  if (!exists("rubin_normal", envir = stanmodels)) {
    # Path to the .stan file (in inst/stan/)
    stan_file <- system.file("stan", "rubin_normal.stan", package = "siteBayes")

    # Compile via rstan
    mod <- rstan::stan_model(file = stan_file)

    # Store in environment
    stanmodels$rubin_normal <- mod
  }

  return(stanmodels$rubin_normal)
}
