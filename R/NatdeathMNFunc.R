#' Nested mosquito death sampler
#'
#' Only for use within multinom()
#' @param mos_Es Matrix, exposed but not yet infectious mosquitoes that have not yet completed their EIP 
#' @param mospdeath single value between 0 and 1, the daily probability of mosquito mortality
#' @details Returns a matrix of dead mosquitoes
#' @export
#' @examples
#' mos_Es = matrix(round(runif(100, 0, 10), 0), nrow = 10, ncol = 10)
#' probs <- function(x) rbeta(x, shape1 = 26.84987, shape2 = 107.3995)
#' natdeath.MN.func(mos_Es, probs)




natdeath.MN.func <- function(mos_Es, probs){
  newmosD = rep(0, length(mos_Es))
  pindex <- mos_Es > 0
  newmosD[pindex] = apply(cbind(mos_Es[pindex], rep(probs, length(mos_Es[pindex]))), 1, function(u) sum(rbinom(u[1], 1, u[2])))
  return(matrix(newmosD, ncol = ncol(mos_Es)))
}