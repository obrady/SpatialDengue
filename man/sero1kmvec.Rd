% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sero1kmvec.R
\docType{data}
\name{sero1kmvec}
\alias{sero1kmvec}
\title{Proportion immune in Singapore in vector format}
\format{A vector of proportion immune per pixel
\describe{
  \item{value}{proportion immune in pixel}
}}
\source{
\url{https://doi.org/10.1093/aje/kwz110}

\url{https://www.singstat.gov.sg/publications/population-trends}
}
\usage{
sero1kmvec
}
\description{
The same as "sero" but aggregated to 1km x 1km (instead of 100m x 100m) using mean valeus over the area and converted to vector format. For use on servers that do not support "rgdal" or "raster" packages
}
\keyword{datasets}
