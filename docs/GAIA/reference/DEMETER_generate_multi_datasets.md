# Generate Multi-Dataset Simulations with Shared Latent Structure

Generates multiple synthetic datasets that share a common latent
variable structure. Used for testing and validating multi-omics
integration methods like ConvCCA and RelPMDCCA.

## Usage

``` r
DEMETER_generate_multi_datasets(
  n_datasets = 3,
  n_samples = 100,
  n_features = c(500, 1000, 750),
  signal_samples = 25,
  signal_features = c(50, 50, 40),
  seed = 1234
)
```

## Arguments

- n_datasets:

  Integer. Number of datasets to generate (default = 3).

- n_samples:

  Integer. Number of samples in each dataset (default = 100).

- n_features:

  Integer vector. Number of features for each dataset. Must have length
  equal to n_datasets (default = c(500, 1000, 750)).

- signal_samples:

  Integer. Number of samples containing true signal in the shared latent
  variable (default = 25). Remaining samples are noise.

- signal_features:

  Integer vector. Number of features containing true signal in each
  dataset. Must have length equal to n_datasets (default = c(50, 50,
  40)).

- seed:

  Integer. Random seed for reproducibility (default = 1234).

## Value

A list with the following components:

- datasets:

  List of length n_datasets containing simulated data matrices. Each
  matrix has dimensions n_samples x n_features.

- true_canonical:

  List of true canonical vectors (one per dataset). Signal features have
  weight 1, noise features have weight 0.

- shared_latent:

  Matrix (n_samples x 1) containing the shared latent variable. First
  signal_samples rows are 1, remaining rows are 0.

- signal_samples:

  Integer indicating number of signal samples.

- signal_features:

  Integer vector indicating number of signal features per dataset.

## Details

The simulation model follows the sparse CCA framework:

For each dataset i: \$\$X_i = u \cdot v_i^T + E_i\$\$

Where:

- u is the shared latent variable (n_samples x 1)

- v_i is the canonical vector for dataset i (n_features x 1)

- E_i is Gaussian noise ~ N(0, 1)

The true signal is sparse:

- Only the first signal_samples samples have non-zero latent values

- Only the first `signal_features[i]` features in each dataset
  participate in the shared structure

This simulation allows controlled testing of:

- Feature selection accuracy (can methods recover signal features?)

- Sample selection accuracy (can methods identify signal samples?)

- Robustness to noise

- Scalability with varying dimensions

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate 3 datasets with default parameters
sim <- DEMETER_generate_multi_datasets()

# Generate 2 high-dimensional datasets
sim <- DEMETER_generate_multi_datasets(
  n_datasets = 2,
  n_samples = 50,
  n_features = c(2000, 3000),
  signal_samples = 10,
  signal_features = c(100, 150),
  seed = 42
)

# Access simulated data
X_1 <- sim$datasets[[1]]
X_2 <- sim$datasets[[2]]
true_weights <- sim$true_canonical[[1]]

} # }
```
