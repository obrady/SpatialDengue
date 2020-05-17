#' Plots a map of a given model compartment state
#'
#' Takes a vector or matrix of a model state, e.g. S, I or R, and plots it back on the original supplied population map
#' @param state Matrix or vector of a model state, e.g. S, I or R
#' @param examras Population raster used to determin the landscape, e.g. sgpop
#' @param unipix Universal pixel lookup table, see ?make.unipix
#' @param main Optional character title
#' @keywords plot map state
#' @export
#' @examples
#' data(sgdat)
#' data(sgpop)
#' sgpop <- pop.process(sgpop, agg = 10)
#' unipix <- make.unipix(sgpop)
#' sgdat <- data.frame(sgdat, patchID = apply(cbind(sgdat[, 3:2]), 1, pix.id.find, unipix))
#' sero <- c(0.577, 0.659, 0.814)
#' startweek <- 40
#' stim <- stim.generate(sgdat, sero, startweek, sgpop, unipix)
#' plot.state(stim)


plot.state <- function(state, examras, unipix, main = NULL){
  # it state is a matrix collapse into a vector
  if(class(state) == "matrix"){
    state = apply(state, 2, sum)
  }
  examras2 = as.vector(examras)
  examras2[unipix$pixID] = state
  values(examras) = examras2
  plot(examras, main = main)
}
