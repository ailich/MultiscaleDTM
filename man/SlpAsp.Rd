% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SlpAsp.R
\name{SlpAsp}
\alias{SlpAsp}
\title{Multiscale Slope and Aspect}
\usage{
SlpAsp(
  r,
  w = c(3, 3),
  unit = "degrees",
  method = "queen",
  metrics = c("slope", "aspect", "eastness", "northness"),
  na.rm = FALSE,
  include_scale = FALSE,
  mask_aspect = TRUE,
  filename = NULL,
  overwrite = FALSE,
  wopt = list()
)
}
\arguments{
\item{r}{DTM as a SpatRaster or RasterLayer in a projected coordinate system where map units match elevation/depth units}

\item{w}{A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.}

\item{unit}{"degrees" or "radians"}

\item{method}{"queen" or "rook", indicating how many neighboring cells to use to compute slope for any cell. queen uses 8 neighbors (up, down, left, right, and diagonals) and rook uses 4 (up, down, left, right). Alternatively, instead of "queen" or "rook", method can be specified as 8 and 4 respectively.}

\item{metrics}{a character string or vector of character strings of which terrain attributes to return. Default is to return all available attributes which are c("slope", "aspect", "eastness", "northness").}

\item{na.rm}{Logical indicating whether or not to remove NA values before calculations. Only applicable if method is "queen" or "8".}

\item{include_scale}{logical indicating whether to append window size to the layer names (default = FALSE)}

\item{mask_aspect}{A logical. When mask_aspect is TRUE (the default), if slope evaluates to 0, aspect will be set to NA and both eastness and northness will be set to 0. When mask_aspect is FALSE, when slope is 0 aspect will be pi/2 radians or 90 degrees which is the behavior of raster::terrain, and northness and eastness will be calculated from that.}

\item{filename}{character Output filename. Can be a single filename, or as many filenames as there are layers to write a file for each layer}

\item{overwrite}{logical. If TRUE, filename is overwritten (default is FALSE).}

\item{wopt}{list with named options for writing files as in writeRaster}
}
\value{
a SpatRaster or RasterStack of slope and/or aspect (and components of aspect)
}
\description{
Calculates multiscale slope and aspect based on the slope.k/aspect.k algorithm from Misiuk et al (2021) which extends classical formulations of slope restricted to a 3x3 window to multiple scales. The code from Misiuk et al (2021) was modified to allow for rectangular rather than only square windows.
}
\details{
When method="rook", slope and aspect are computed according to Fleming and Hoffer (1979) and Ritter (1987). When method="queen", slope and aspect are computed according to Horn (1981). These are the standard slope algorithms found in many GIS packages but are traditionally restricted to a 3 x 3 window size. Misiuk et al (2021) extended these classical formulations  to multiple window sizes. This function modifies the code from Misiuk et al (2021) to allow for rectangular rather than only square windows and also added aspect.
}
\examples{
r<- rast(volcano, extent= ext(2667400, 2667400 + 
ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), 
crs = "EPSG:27200")
slp_asp<- SlpAsp(r = r, w = c(5,5), unit = "degrees", 
method = "queen", metrics = c("slope", "aspect", 
"eastness", "northness"))
plot(slp_asp)
}
\references{
Fleming, M.D., Hoffer, R.M., 1979. Machine processing of landsat MSS data and DMA topographic data for forest cover type mapping (No. LARS Technical Report 062879). Laboratory for Applications of Remote Sensing, Purdue University, West Lafayette, Indiana.

Horn, B.K., 1981. Hill Shading and the Reflectance Map. Proceedings of the IEEE 69, 14-47.

Misiuk, B., Lecours, V., Dolan, M.F.J., Robert, K., 2021. Evaluating the Suitability of Multi-Scale Terrain Attribute Calculation Approaches for Seabed Mapping Applications. Marine Geodesy 44, 327-385. https://doi.org/10.1080/01490419.2021.1925789

Ritter, P., 1987. A vector-based slope and aspect generation algorithm. Photogrammetric Engineering and Remote Sensing 53, 1109-1111.
}
