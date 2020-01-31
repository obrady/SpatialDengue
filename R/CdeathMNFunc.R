#' Nested mosquito death sampler (due to vector control)
#'
#' Only for use within multinom()
#' @param mos_Es Matrix or vector, matrix = exposed but not yet infectious mosquitoes that have not yet completed their EIP, vector = infectious mosquitoes 
#' @param probs vector, patch-specific daily probability mosqutio dies due to reactive vector control
#' @details Returns a matrix of dead mosquitoes
#' @export
#' @examples
#' mos_Es = matrix(round(runif(100, 0, 10), 0), nrow = 10, ncol = 10)
#' probs <- c(rep(0.5, 5), rep(0, 5))
#' Cdeath.MN.func(mos_Es, probs)
#' 
#' # infectious mosquitoes
#' mos_Is = round(runif(10, 0, 10), 0)
#' Cdeath.MN.func(mos_Is, probs)




Cdeath.MN.func <- function(mos_Es, probs){
  # probability of dying dure to vector control doesn't depend on how old a mosquito is
  if(class(mos_Es) == "matrix"){
    probs = matrix(rep(probs, nrow(mos_Es)), ncol = ncol(mos_Es), byrow = T)
  }
  newmosD2 = rep(0, length(mos_Es))
  if(sum(probs) > 0){
    pindex <- mos_Es > 0
    newmosD2[pindex] = apply(cbind(mos_Es[pindex], probs[pindex]), 
                             1, function(u) sum(rbinom(u[1], 1, u[2])))
  }
  if(class(mos_Es) == "matrix"){
    return(matrix(newmosD2, ncol = ncol(mos_Es)))
  }else{return(newmosD2)}
}