#' A universal lookup table for pixels and patches in the focus area
#'
#' A lookup table where each row is a patch in the model, with a corresponding pixel ID in the "sgpop" raster, a popultion count and a latitude and longitude. 
#' Comes from the make.unipix() function but a precalcualted version for Singapore (at 1km x 1km patch resolution) is included for use on the server that avoids having to deal with rasters directly
#'
#' @format A data frame where each row is a patch of the model
"unipix"