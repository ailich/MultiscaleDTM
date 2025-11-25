# Creates circular focal window

Creates circular focal window around central pixel.

## Usage

``` r
circle_window(radius, unit, resolution, return_dismat = FALSE)
```

## Arguments

- radius:

  radius of circular window

- unit:

  unit for radius. Either "cell" (number of cells) or "map" for map
  units (e.g. meters).

- resolution:

  resolution of intended raster layer (one number or a vector of length
  2). Only necessary if unit= "map"

- return_dismat:

  logical, if TRUE return a matrix of distances from focal cell instead
  of a matrix to pass to terra::focal.

## Value

a matrix of 1's and NA's showing which cells to include and exclude
respectively in focal calculations, or if return_dismat=TRUE, a matrix
indicating the distance from the focal cell. It also contains attributes
attributes 'unit', 'scale', and 'shape' if return_dismat=FALSE, and if
return_dismat=TRUE the attribute 'unit'.
