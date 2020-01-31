#' Radiation model
#'
#' Human movement model based on radiation, for details see Simini et al. (2012) Nature. 484, 96-100
#' @param m human population size in patch m (home patch)
#' @param n human population size in patch n
#' @param s human population size in all pixels within the radius m -> n, excluding m and n
#' @details returns the predicted flux between pixels m and n
#' @keywords movement
#' @export
#' @examples
#' popmat <- matrix(round(runif(16, 5, 20), 0), nrow = 4, ncol = 4)
#' m <- popmat[1, 1]
#' n <- popmat[3, 2]
#' distm <- as.matrix(dist(popmat)) # creates a distance matrix
#' r <- distm[3, 2]
#' s <- sum(popmat[distm <= r]) - m - n
#' rad.model(m, n, s)




rad.model<- function(m, n, s){if(m*n == 0){return(0)}else{return(m * n / ((m + s) * (m + n + s)))}}
