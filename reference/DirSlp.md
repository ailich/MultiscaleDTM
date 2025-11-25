# Directional Slope

Calculates the slope along a specified direction. Upslope values are
positive and downslope values are negative.

## Usage

``` r
DirSlp(
  alpha,
  dz.dx,
  dz.dy,
  unit = "degrees",
  abs = FALSE,
  include_dir = FALSE,
  filename = NULL,
  overwrite = FALSE,
  wopt = list()
)
```

## Arguments

- alpha:

  Angle (in specified 'unit') at which you would like to calculate
  slope. 0 represents up in map direction (usually North) and it
  increases clockwise. This can be a single number or it can be a raster
  of cell values.

- dz.dx:

  The change in elevation per unit distance in the x direction as a
  SpatRaster, RasterLayer, or a single number. Positive is to the right.
  See details for more.

- dz.dy:

  The change in elevation per unit distance in the y direction as a
  SpatRaster or RasterLayer,or a single number Positive is up. See
  details for more.

- unit:

  "degrees" or "radians" (default is "degrees")

- abs:

  logical indicating whether or not to return the absolute value of
  slope (default is FALSE)

- include_dir:

  logical indicating whether to append direction to layer name (default
  is FALSE)

- filename:

  character Output filename. Can be a single filename, or as many
  filenames as there are layers to write a file for each layer

- overwrite:

  logical. If TRUE, filename is overwritten (default is FALSE).

- wopt:

  list with named options for writing files as in writeRaster

## Value

a SpatRaster or RasterStack of slope and/or aspect (and components of
aspect)

## Details

dz.dx and dz.dy can be calculated at a specified scale via `SlpAsp`,
`Pfit`, `Qfit` (zx and zy), or from an existing layer calculated by
another program.

## References

Neteler, M., & Mitasova, H. (2008). Open source GIS: A GRASS GIS
approach (3rd ed.). Springer.

## Examples

``` r
r<- erupt()
dz1<- SlpAsp(r, metrics = c("dz.dx", "dz.dy"))
dz2<- Qfit(r, metrics = c(), return_params = TRUE, as_derivs=TRUE)
dz3<- Pfit(r, metrics = c("dz.dx", "dz.dy"))
dirslp1<- DirSlp(alpha = 45, dz.dx= dz1$dz.dx, dz.dy= dz1$dz.dy)
dirslp2<- DirSlp(alpha = 45, dz.dx= dz2$zx, dz.dy= dz2$zy)
dirslp3<- DirSlp(alpha = 45, dz.dx= dz3$dz.dx, dz.dy= dz3$dz.dy)
```
