# Implementation of the Sappington et al., (2007) vector ruggedness measure

Implementation of the Sappington et al., (2007) vector ruggedness
measure, modified from Evans (2021).

## Usage

``` r
VRM(
  r,
  w = c(3, 3),
  na.rm = FALSE,
  include_scale = FALSE,
  filename = NULL,
  overwrite = FALSE,
  wopt = list()
)
```

## Arguments

- r:

  DTM as a SpatRaster or RasterLayer

- w:

  A vector of length 2 specifying the dimensions of the rectangular
  window to use where the first number is the number of rows and the
  second number is the number of columns. Window size must be an odd
  number. Default is 3x3.

- na.rm:

  A logical indicating whether or not to remove NA values before
  calculations. See details for more information.

- include_scale:

  logical indicating whether to append window size to the layer names
  (default = FALSE)

- filename:

  character Output filename.

- overwrite:

  logical. If TRUE, filename is overwritten (default is FALSE).

- wopt:

  list with named options for writing files as in writeRaster

## Value

a RasterLayer

## Details

If the crs is cartesian, when na.rm=TRUE, NA's will be removed from the
slope/aspect calculations. When the crs is lat/lon, na.rm=TRUE will not
affect the calculation of slope/aspect as terra::terrain will be used
since it can calculate slope and aspect for spherical geometry but it
does not support na.rm. In both cases when na.rm=TRUE, the x, y, and z
components will be summed with na.rm=TRUE, and the N used in the
denominator of the VRM equation will be the number of non-NA cells in
the window rather than the total number of cells.

## References

Evans JS (2021). spatialEco. R package version 1.3-6,
https://github.com/jeffreyevans/spatialEco.

Sappington, J.M., Longshore, K.M., Thompson, D.B., 2007. Quantifying
Landscape Ruggedness for Animal Habitat Analysis: A Case Study Using
Bighorn Sheep in the Mojave Desert. The Journal of Wildlife Management
71, 1419-1426. https://doi.org/10.2193/2005-723

## Examples

``` r
r<- erupt()
vrm<- VRM(r, w=c(5,5), na.rm = TRUE)
plot(vrm)
```
