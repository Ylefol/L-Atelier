# Apply SCAD Penalty Thresholding

Smoothly Clipped Absolute Deviation (SCAD) penalty. Treats small,
medium, and large coefficients differently:

- Small coefficients: shrunk like LASSO

- Medium coefficients: partially shrunk

- Large coefficients: no penalty

## Usage

``` r
MINERVA_soft_threshold_scad(x, tau, a = 3.7)
```

## Arguments

- x:

  Input value or vector

- tau:

  Threshold parameter

- a:

  SCAD shape parameter (default = 3.7, from Fan & Li 2001)

## Value

SCAD-thresholded value(s)
