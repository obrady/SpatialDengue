#' Matrix time transition function
#'
#' As rows indicate time in each state, at the end of the timestep, each of the matrix-based states need to be moved down a row
#' @param base matrix, the model state, e.g. mos_E
#' @param plus vector, new additions to the state, e.g. newmosE
#' @param minus matrix, new subtractions from the state, e.g. newmosD
#' @details Removes those who graduate from the state, shuffles down the rows by one then adds new additions to the first row of the matrix
#' @export
#' @examples
#' mos_Es = matrix(round(runif(100, 5, 10), 0), nrow = 10, ncol = 10)
#' newmosE = round(runif(10, 0, 2), 0)
#' newmosD = matrix(round(runif(100, 0, 1), 0), nrow = 10, ncol = 10)
#' mos_EsT2 <- transition(mos_Es, newmosE, newmosD)

transition <- function(base, plus, minus){
  base2 = base - minus
  base2[2:nrow(base2), ] = base2[1:(nrow(base2) - 1), ]
  base2[1, ] = plus
  return(base2)
}