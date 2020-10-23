#' Proportion immune in Singapore in vector format
#'
#' The same as "sero" but aggregated to 1km x 1km (instead of 100m x 100m) using mean valeus over the area and converted to vector format. For use on servers that do not support "rgdal" or "raster" packages
#'
#' @format A vector of proportion immune per pixel
#' \describe{
#'   \item{value}{proportion immune in pixel}
#' }
#' @source \url{https://doi.org/10.1093/aje/kwz110}
#' @source \url{https://www.singstat.gov.sg/publications/population-trends}
"sero1kmvec"