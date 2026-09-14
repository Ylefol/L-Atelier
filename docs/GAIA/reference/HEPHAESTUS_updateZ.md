# Update Consensus Variable (z-update)

Projects onto the unit ball for consensus constraint in ADMM. Ensures
\|\|z\|\| \<= 1 by normalization if needed.

## Usage

``` r
HEPHAESTUS_updateZ(x, old)
```

## Arguments

- x:

  Current value before projection

- old:

  Previous value (fallback if NA encountered)

## Value

Projected consensus variable z
