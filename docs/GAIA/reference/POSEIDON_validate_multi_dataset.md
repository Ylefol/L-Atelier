# Validate Multi-Dataset Structure

Checks that all datasets in a list have:

- Same number of samples (rows)

- No missing values

- Valid dimensions

## Usage

``` r
POSEIDON_validate_multi_dataset(X_list, sample_names = NULL)
```

## Arguments

- X_list:

  List of data matrices

- sample_names:

  Optional vector of expected sample names

## Value

Invisible TRUE if valid, stops with error otherwise
