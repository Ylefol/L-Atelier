# Build a hephaestus_mofa result object from a trained/loaded MOFA model

Build a hephaestus_mofa result object from a trained/loaded MOFA model

## Usage

``` r
.hephaestus_mofa_build_result(trained, params = NULL)
```

## Details

MOFA2::get_factors()/get_weights()/calculate_variance_explained() shapes
used below (list-per-group of sample x factor / feature x factor / view
x factor matrices) match MOFA2's own documentation and vignettes as of
this writing – not independently verified against every MOFA2 release,
so confirm dimnames/orientation on first real run if a future MOFA2
version changes them.
