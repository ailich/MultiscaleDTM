# Calculates standard deviation of bathymetry (a measure of rugosity) adjusted for slope

Calculates standard deviation of bathymetry (a measure of rugosity).
Using a sliding rectangular window a plane is fit to the data and the
standard deviation of the residuals is calculated (Ilich et al., 2023)

## Usage

``` r
AdjSD(
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

  DTM as a SpatRaster or RasterLayer in a projected coordinate system
  where map units match elevation/depth units

- w:

  A vector of length 2 specifying the dimensions of the rectangular
  window to use where the first number is the number of rows and the
  second number is the number of columns. Window size must be an odd
  number.

- na.rm:

  A logical indicating whether or not to remove NA values before
  calculations

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

a SpatRaster or RasterLayer of adjusted rugosity

## References

Ilich, A. R., Misiuk, B., Lecours, V., & Murawski, S. A. (2023).
MultiscaleDTM: An open-source R package for multiscale geomorphometric
analysis. Transactions in GIS, 27(4). https://doi.org/10.1111/tgis.13067

## Examples

``` r
r<- erupt()
adjsd<- AdjSD(r, w=c(5,5), na.rm = TRUE)
plot(adjsd)
```
