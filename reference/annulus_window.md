# Creates annulus focal window

Creates annulus focal window around central pixel.

## Usage

``` r
annulus_window(radius, unit, resolution)
```

## Arguments

- radius:

  radius of inner annulus c(inner,outer). Inner radius must be less than
  or equal to outer radius.

- unit:

  unit for radius. Either "cell" (number of cells, the default) or "map"
  for map units (e.g. meters).

- resolution:

  resolution of intended raster layer (one number or a vector of length
  2). Only necessary if unit= "map"

## Value

a matrix of 1's and NA's showing which cells to include and exclude
respectively in focal calculations. It also contains attributes
attributes 'unit', 'scale', and 'shape'.
