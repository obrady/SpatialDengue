#' Nested Human drug treatment sampler
#'
#' Only for use within multinom()
#' @param sing_I Matrix, Infectious humans, each row is days since becoming infectious
#' @param probs Function, giving the probability of completing EIP on a given day
#' @details Returns a matrix of newly prophylactic drug treated humans
#' @export
#' @examples
#' sing_I = matrix(round(runif(100, 0, 10), 0), nrow = 10, ncol = 10)
#' probs <- c(rep(0.5, 5), rep(0, 5))
#' Hum.Drug.Treat.MN(sing_I, probs)

Hum.Drug.Treat.MN <- function(sing_I, probs){
  probs = matrix(rep(probs, nrow(sing_I)), ncol = ncol(sing_I), byrow = T)
  if(sum(probs) > 0){
    return(matrix(apply(cbind(as.vector(sing_I), as.vector(probs)), 1, function(u) rbinom(1, u[1], u[2])), ncol = ncol(sing_I)))
  }else{return(matrix(0, ncol = ncol(sing_I), nrow = nrow(sing_I)))}
}