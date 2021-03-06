% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sgpop1kmvec.R
\docType{data}
\name{sgpop1kmvec}
\alias{sgpop1kmvec}
\title{Population of Singapore in vector format}
\format{A vector of population values per pixel
\describe{
  \item{value}{number of people residing in pixel}
}}
\source{
\url{https://www.openstreetmap.org/#map=6/54.910/-3.432}
}
\usage{
sgpop1kmvec
}
\description{
The same as "sgpop" but aggregated to 1km x 1km (instead of 100m x 100m) using the pop.process function and converted to vector format. For use on serveres that do not support "rgdal" or "raster" packages
}
\keyword{datasets}
