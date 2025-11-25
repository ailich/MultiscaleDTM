# Helper function to filter outliers from regression parameters using interquartile range

Helper function to filter outliers from regression parameters using
interquartile range

## Usage

``` r
outlier_filter(params, outlier_quantile, wopt = list())
```

## Arguments

- params:

  regression parameters for fitted surface

- outlier_quantile:

  A numeric vector of length two or three. If two numbers are used it
  specifies the lower (Q1) and upper (Q2) quantiles used for determining
  the interquantile range (IQR). These values should be between 0 and 1
  with Q2 \> Q1. An optional third number can be used to specify a the
  size of a regular sample to be taken which can be useful if the full
  dataset is too large to fit in memory. Values are considered outliers
  if they are less than Q1-(100*IQR) or greater than Q2+(100*IQR), where
  IQR=Q2-Q1.

- wopt:

  list with named options for writing files as in writeRaster
