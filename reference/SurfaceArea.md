# Calculates surface area of a DTM

Calculates surface area on a per cell basis of a DTM based on Jenness,
2004.

## Usage

``` r
SurfaceArea(
  r,
  na.rm = FALSE,
  filename = NULL,
  overwrite = FALSE,
  wopt = list()
)
```

## Arguments

- r:

  DTM as a SpatRaster or RasterLayer in a projected coordinate system
  where map units match elevation/depth units

- na.rm:

  Logical indicating whether to remove NAs from calculations. When
  FALSE, the sum of the eight triangles is calculated. When TRUE, the
  mean of the created triangles is calculated and multiplied by 8 to
  scale it to the proper area.

- filename:

  character Output filename.

- overwrite:

  logical. If TRUE, filename is overwritten (default is FALSE).

- wopt:

  list with named options for writing files as in writeRaster

## Value

a SpatRaster or RasterLayer

## References

Jenness, J.S., 2004. Calculating landscape surface area from digital
elevation models. Wildlife Society Bulletin 32, 829-839.

## Examples

``` r
r<- erupt()
sa<- SurfaceArea(r)
plot(sa)
```
