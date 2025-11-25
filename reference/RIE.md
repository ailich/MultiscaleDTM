# Calculates Roughness Index-Elevation

Calculates Roughness Index-Elevation. This is the standard deviation of
residual topography in a focal window where residual topography is
calculated as the focal pixel minus the focal mean.

## Usage

``` r
RIE(
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
  calculation of SD

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

a SpatRaster or RasterLayer

## Details

Note the original paper by Cavalli et al (2008) uses a fixed 5x5 window
and uses 25 as the denominator indicating use of the population standard
deviation. This implementation provides a flexible window size and
istead calculates the sample standard deviation which uses a denominator
of n-1.

## References

Cavalli, M., Tarolli, P., Marchi, L., Dalla Fontana, G., 2008. The
effectiveness of airborne LiDAR data in the recognition of channel-bed
morphology. CATENA 73, 249–260.
https://doi.org/10.1016/j.catena.2007.11.001

## Examples

``` r
r<- erupt()
rie<- RIE(r, w=c(5,5), na.rm = TRUE)
plot(rie)
```
