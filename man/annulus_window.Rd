% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BPI.R
\name{annulus_window}
\alias{annulus_window}
\title{Creates annulus focal window}
\usage{
annulus_window(radius, unit = "cell", resolution, return_dismat = FALSE)
}
\arguments{
\item{radius}{radius of inner annulus c(inner,outer)}

\item{unit}{unit for radius. Either "cell" (number of cells, the default) or "map" for map units (e.g. meters).}

\item{resolution}{resolution of intended raster layer (one number or a vector of length 2). Only necessary if unit= "map"}

\item{return_dismat}{logical, if TRUE return a matrix of distances from focal cell instead of a matrix to pass to terra::focal (default FALSE)}
}
\value{
if a matrix of 1's and NA's showing which cells to include and exclude respectively in focal calculations, or if return_dismat=TRUE, a matrix indicating the distance from the focal cell.
}
\description{
Creates annulus focal window around central pixel.
}