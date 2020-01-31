#' Match Longitude and Latitude of clusters to pixel IDs in grid
#'
#' This function takes latitude and longitude of identified cases and matched them to a specific identifier from the unipix dataframe
#' @param longlat A two element vector of Longitude followed by Latitude
#' @param unipix A universal lookup table detailing all pixels within the grid. Must contain a pixel ID in the first column and the fields "x" (Longitude) and "y" (Latitude)
#' @details uses pixel stepwise distance (Longitude + Latitude) to identify nearest pixel to the given LatLong.
#' @keywords pixel
#' @export
#' @examples
#' data(sgdat)
#' data(sgpop)
#' unipix <- make.unipix(sgpop)
#' pix.id.find(c(sgdat$Longitude[1], sgdat$Latitude[1]), unipix)



# function for pixel finding given a lat long
pix.id.find <- function(longlat, unipix){
  # trim to only those pixels with resident population
  unipix2 = unipix[unipix$pop > 0, ]
  xclose = (unipix2$x - longlat[1])^2
  yclose = (unipix2$y - longlat[2])^2
  #far = distm(longlat, cbind(unipix$x, unipix$y))
  return(as.numeric(unipix2[which.min(xclose + yclose), 1]))
}