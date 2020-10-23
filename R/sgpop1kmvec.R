#' Population of Singapore in vector format
#'
#' The same as "sgpop" but aggregated to 1km x 1km (instead of 100m x 100m) using the pop.process function and converted to vector format. For use on serveres that do not support "rgdal" or "raster" packages
#'
#' @format A vector of population values per pixel
#' \describe{
#'   \item{value}{number of people residing in pixel}
#' }
#' @source \url{https://www.openstreetmap.org/#map=6/54.910/-3.432}
"sgpop1kmvec"