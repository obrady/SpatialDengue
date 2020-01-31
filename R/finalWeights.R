#' Posterior distribution of fitted model parameters
#'
#' A list of data.frames detailing the likelihood of each parameter value
#' '
#' @format A list of length four (each parameter), within each contains a data.frame composed of 10000 rows and 2 columns
#' \describe{
#'   \item{names}{parameter each posterior is related to}
#'   \item{val}{value of the parameter at that sample}
#'   \item{prob}{Likelihood of the data given the model with that parameter value}
#'   }
#' @source {contact author for more detail on fitting}
"finalWeights"