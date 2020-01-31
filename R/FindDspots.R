#' Find detected spots
#'
#' Based on the location of detected (symptomatic) infections and a specified radius, identify which patches aroudn their residence will be subject to control activities
#' @param newsingD Matrix of newly detected dengue infections, row = day of infection when they were detected, col = the patch in which they were detected
#' @param pixdistmat A distance matrix between patches
#' @param radius The radius in meters around the index case's residence that is to be searched
#' @details Returns a numeric vector with patch IDs that fall within the specified radius of an index case
#' @keywords index radius
#' @export
#' @examples
#' data(sgpop)
#' sgpop <- pop.process(sgpop, agg = 10)
#' unipix <- make.unipix(sgpop)
#' pixdistmat <- distm(cbind(unipix$x, unipix$y))
#' newsingD <- matrix(sample(c(0, 1), nrow(unipix) * 2, prob = c(0.99, 0.01), T), nrow = 2)
#' radius = 1000
#' find.Dspots(newsingD, pixdistmat, radius)


#find.Dspots <- function(newsingD, pixdistmat, radius){return(newsingD)}


find.Dspots <- function(newsingD, pixdistmat, radius){
  if(radius > 0){
    activespots <- (1:ncol(newsingD))[colSums(newsingD) > 0]
    if(sum(activespots) > 0){
      pixdistmatF <- pixdistmat[activespots, ]
      pixdistmatF <- pixdistmatF <= radius
      if(is.vector(pixdistmatF)){pixdistmatF = matrix(pixdistmatF, nrow = 1)}
      chosenspots <- (1:ncol(newsingD))[colSums(pixdistmatF) > 0]
      return(chosenspots)
    }else{return(NA)}
  }else{
    nsdcols <- colSums(newsingD) > 0
    if(sum(nsdcols) > 0){
      return((1:ncol(newsingD))[nsdcols])
    }else{return(NA)}
    }
}