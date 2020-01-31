#' Population processing
#'
#' This function processes a raster of human population (counts) to be compatible with the dengue model
#' @param sgpop A raster of human population counts
#' @param agg Multiplier for level of aggregation of the raster, defaults to FALSE (0)
#' @details Processes performed involve 
#' i) making all pixels integer value (so each pixel contains a discrete number of people)
#' ii) if supplied "agg" determines the scale of pixel aggregation, e.g. agg = 10 will convert a raster with pixels sized 100m x 100m to 1000m x 1000m
#' @keywords population
#' @export
#' @examples
#' data(sgpop)
#' sgpop10 <- pop.process(sgpop, agg = 10)
#' par(mfrow = c(2, 1))
#' plot(sgpop)
#' plot(sgpop10)



# processes the population surfaces so they are compatible with the model
pop.process <- function(sgpop, agg = FALSE){
  # make integer
  sgpop = round(sgpop, 0)
  # aggregate if supplied
  if(agg != FALSE){sgpop = aggregate(sgpop, agg, fun = sum)}
  return(sgpop)
}