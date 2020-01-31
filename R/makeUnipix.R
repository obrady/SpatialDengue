#' Create a universal lookup table from a population raster
#'
#' This function takes takes a raster of population and creates a data.frame with key details including unique identifiers of the pixels and their characteristics including Longitude and Latitude and population.
#' This is used later on by the model for reference
#' @param popras A population raster, see ?sgpop
#' @keywords unipix
#' @export
#' @examples
#' data(sgpop)
#' unipix <- make.unipix(sgpop)

make.unipix <- function(popras){
  unipix = data.frame(pixID = (1:length(as.vector(popras)))[!is.na(as.vector(popras))],
                      pop = as.vector(popras)[!is.na(as.vector(popras))])
  unipix = unipix[unipix$pop > 0, ]
  unipix = data.frame(patchID = 1:nrow(unipix), unipix, coordinates(popras)[unipix$pixID, ])
}