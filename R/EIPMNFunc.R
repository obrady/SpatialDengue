#' Nested EIP completion sampler
#'
#' Only for use within multinom()
#' @param mos_Es Matrix, exposed but not yet infectious mosquitoes that have not yet died 
#' @param probs Function, giving the probability of completing EIP on a given day
#' @details Returns a matrix of newly infectious mosquitoes
#' @export
#' @examples
#' mos_Es = matrix(round(runif(100, 0, 10), 0), nrow = 10, ncol = 10)
#' probs <- function(times) plnorm(times, meanlog = log(7), sdlog = 0.21)
#' EIP.MN.func(mos_Es, probs)



EIP.MN.func <- function(mos_Es, probs){
  newmosI = matrix(0, ncol(mos_Es), nrow = nrow(mos_Es))
  for(k in 1:nrow(mos_Es)){
    pindex <- mos_Es[k, ] > 0
    if(sum(pindex) > 0){
      newmosI[k, pindex] = apply(cbind(mos_Es[k, ][pindex], probs(k)), 1, function(u) sum(rbinom(u[1], 1, u[2])))
    }
  }
  return(newmosI)
}