#' Seroprevalence of dengue in Singapore
#'
#' A dataset containing the number of people with IgG prevalence to at lease one dengue virus serotype in Singapore. Uses data from a 2013 seroprevalence survey to calculate age-specific
#' seroprevalence which is then extrapolated across singapore using age data from the census
#'
#' @format A georeferenced raster file with pixels aligned in a 184 row, 271 column grid:
#' \describe{
#'   \item{value}{The proportion of people in each pixel that are IgG positive}
#' }
#' @source \url{https://doi.org/10.1093/aje/kwz110}
#' @source \url{https://www.singstat.gov.sg/publications/population-trends}
"sero"