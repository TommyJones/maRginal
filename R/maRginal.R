#' maRginal
#'
#' Marginal Effects for Arbitrary Predictive Models
#'
#' Calculate marginal effects, similar to linear coefficients, for any
#' predictive model. Also calculates and plots individual conditional expectation
#' (ICE) curves. As of this writing, the concept is experimental. Caveat emptor!
#'
#' @importFrom textmineR TmParallelApply
#' @importFrom grDevices rgb
#' @importFrom graphics lines par plot points
#' @importFrom stats model.matrix predict sd

#' @docType package
#' @name maRginal
NULL
