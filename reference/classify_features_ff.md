# Helper function factory to classify morphometric features

Helper function factory to classify morphometric features according to a
modified version of Wood 1996 page 120

## Usage

``` r
classify_features_ff(slope_tolerance = 1, curvature_tolerance = 1e-04)
```

## Arguments

- slope_tolerance:

  Slope tolerance that defines a 'flat' surface (degrees; default is
  1.0). Relevant for the features layer.

- curvature_tolerance:

  Curvature tolerance that defines 'planar' surface (default is 0.0001).
  Relevant for the features layer.

## Value

A function that can be passed to raster::overlay to classify
morphometric features

## References

Wood, J., 1996. The geomorphological characterisation of digital
elevation models (Ph.D.). University of Leicester.
