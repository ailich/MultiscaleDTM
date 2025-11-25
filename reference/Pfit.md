# Calculates multiscale slope and aspect using a local planar fit.

Calculates multiscale slope and aspect of a DTM over a sliding
rectangular window using a using a local planar fit to the surface
(Sharpnack and Akin 1969).

## Usage

``` r
Pfit(
  r,
  w = c(3, 3),
  unit = "degrees",
  metrics = c("pslope", "paspect", "peastness", "pnorthness"),
  na.rm = FALSE,
  include_scale = FALSE,
  mask_aspect = TRUE,
  filename = NULL,
  overwrite = FALSE,
  wopt = list()
)
```

## Arguments

- r:

  DTM as a SpatRaster (terra) or RasterLayer (raster) in a projected
  coordinate system where map units match elevation/depth units (up is
  assumed to be north for calculations of aspect, northness, and
  eastness).

- w:

  Vector of length 2 specifying the dimensions of the rectangular window
  to use where the first number is the number of rows and the second
  number is the number of columns. Window size must be an odd number.
  Default is 3x3.

- unit:

  "degrees" or "radians".

- metrics:

  Character vector specifying which terrain attributes to return. The
  default is to return c("pslope", "paspect", "peastness", and
  "pnorthness"). These are preceded with a 'p' to differentiate them
  from the measures calculated by SlpAsp() and 'Qfit' where the 'p'
  indicates that a planar surface was used for the calculation.
  Additional measures available include "dz.dx" and "dz.dy" which are
  the x and y components of slope respectively.

- na.rm:

  Logical indicating whether or not to remove NA values before
  calculations.

- include_scale:

  logical indicating whether to append window size to the layer names
  (default = FALSE)

- mask_aspect:

  Logical. If TRUE (default), aspect will be set to NA and northness and
  eastness will be set to 0 when slope = 0. If FALSE, aspect is set to
  270 degrees or 3*pi/2 radians ((-pi/2)- atan2(0,0)+2*pi) and northness
  and eastness will be calculated from this.

- filename:

  character Output filename. Can be a single filename, or as many
  filenames as there are layers to write a file for each layer

- overwrite:

  logical. If TRUE, filename is overwritten (default is FALSE).

- wopt:

  list with named options for writing files as in writeRaster

## Value

a SpatRaster (terra) or RasterStack/RasterLayer (raster)

## Details

If only first order derivatives are needed, Pfit is faster than Qfit and
should provide equivalent results to Qfit for first order derivatives
(Jones, 1998) when na.rm=FALSE and approximately the same results
otherwise.

## References

Evans, I.S., 1980. An integrated system of terrain analysis and slope
mapping. Zeitschrift f¨ur Geomorphologic Suppl-Bd 36, 274–295. Jones, K.
H. (1998). A comparison of algorithms used to compute hill slope as a
property of the DEM. Computers & Geosciences, 24(4), 315–323.
https://doi.org/10.1016/S0098-3004(98)00032-6 Wood, J., 1996. The
geomorphological characterisation of digital elevation models (Ph.D.).
University of Leicester.

## Examples

``` r
r<- erupt()
pmetrics<- Pfit(r, w = c(5,5), unit = "degrees", na.rm = TRUE)
plot(pmetrics)
```
