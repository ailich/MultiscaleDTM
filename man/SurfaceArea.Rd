% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SurfaceArea.R
\name{SurfaceArea}
\alias{SurfaceArea}
\title{Calculates surface area of a DTM}
\usage{
SurfaceArea(
  r,
  na.rm = FALSE,
  filename = NULL,
  overwrite = FALSE,
  wopt = list()
)
}
\arguments{
\item{r}{DTM as a SpatRaster or RasterLayer in a projected coordinate system where map units match elevation/depth units}

\item{na.rm}{Logical indicating whether to remove NAs from calculations. When FALSE, the sum of the eight triangles is calculated. When TRUE, the mean of the created triangles is calculated and multiplied by 8 to scale it to the proper area.}

\item{filename}{character Output filename.}

\item{overwrite}{logical. If TRUE, filename is overwritten (default is FALSE).}

\item{wopt}{list with named options for writing files as in writeRaster}
}
\value{
a SpatRaster or RasterLayer
}
\description{
Calculates surface area on a per cell basis of a DTM based on Jenness, 2004.
}
\examples{
r<- erupt()
sa<- SurfaceArea(r)
plot(sa)
}
\references{
Jenness, J.S., 2004. Calculating landscape surface area from digital elevation models. Wildlife Society Bulletin 32, 829-839.
}
