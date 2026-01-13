# GAIA Setup Guide

## Prerequisites

## Function Naming Convention

All functions in GAIA follow a strict naming convention:

```
MODULE_function_name
```

Where `MODULE` is the subordinate module name in uppercase:

- **DEMETER** - Utilities (matrix operations)
- **POSEIDON** - Preprocessing & Normalization
- **MINERVA** - Machine Learning, Validation, & Penalties
- **HEPHAESTUS** - Integration Methods
- **AETHER** - Visualization
- **ARTEMIS** - Analysis & Statistics
- **APOLLO** - Annotation & Enrichment
- **ELEUTHIA** - Import, Export & Reporting
- **HADES** - Quality Control & Filtering

### Examples:

```r
# Demeter (utilities)
DEMETER_matrix_power(Sigma, -0.5)  # Function form
Sigma %^% (-0.5)                   # Operator form (also available)
DEMETER_gettingInverse(X)

# Poseidon (preprocessing)
POSEIDON_scale_data(X)
POSEIDON_filter_zero_variance(X_list)
POSEIDON_validate_multi_dataset(X_list)

# Minerva (penalties and CV)
MINERVA_soft_threshold_lasso(x, tau)
MINERVA_cv_scca(X_1, X_2, method = "ConvCCA", tau_grid)

# Hephaestus (integration)
HEPHAESTUS_convCCA(X_1, X_2, tauW_1 = 0.3, tauW_2 = 0.3)
HEPHAESTUS_multi_relPMDCCA(X_list, lambda = 10, tau = list(0.8, 0.8, 0.8))

# Aether (visualization)
AETHER_plot_feature_weights(weights, method_name = "ConvCCA")
```

## Usage Example

### Basic Two-Dataset Integration with Cross-Validation

```r
# Functions imported via source

# Prepare data
X_1 <- your_rnaseq_data  # n samples x p1 features
X_2 <- your_atacseq_data # n samples x p2 features

# Preprocess
X_1 <- POSEIDON_filter_zero_variance(X_1)
X_2 <- POSEIDON_filter_zero_variance(X_2)

# Run cross-validation
cv_result <- MINERVA_cv_scca(
    X_1 = X_1,
    X_2 = X_2,
    method = "ConvCCA",
    tau_grid = expand.grid(
        tau1 = seq(0.1, 0.5, 0.1),
        tau2 = seq(0.1, 0.5, 0.1)
    ),
    k = 5,
    nIter = 100,
    penalty = "LASSO"
)

# Visualize results
p <- AETHER_plot_feature_weights(
    cv_result$best_model$W[[1]],
    method_name = "ConvCCA - RNA-seq"
)
print(p)
```

### Multi-Dataset Integration

```r
# Functions imported via source

# Prepare multi-omics data
X_list <- list(rnaseq, atacseq, cutandtag)

# Preprocess - filter features with zero variance in ANY dataset
X_list <- POSEIDON_filter_shared_variance(X_list, by_row = FALSE)
POSEIDON_validate_multi_dataset(X_list)

# Run multi-dataset sCCA
result <- HEPHAESTUS_multi_convCCA(
    X = X_list,
    tau = list(0.3, 0.3, 0.3),
    nIter = 100,
    penalty = "LASSO"
)

# Visualize
p <- AETHER_plot_multi_feature_weights(
    weights_list = result$W,
    method_name = "multi.convCCA",
    dataset_names = c("RNA-seq", "ATAC-seq", "CUT&TAG")
)
print(p)
```

## Module Dependencies

The dependency structure is:

```
Demeter (no dependencies)
    ↓
Minerva/penalty_functions, Poseidon (use Demeter)
    ↓
Hephaestus (uses Demeter + Minerva/penalty_functions)
    ↓
Minerva/sCCA_CV (uses Hephaestus)

Aether (no dependencies)
```

## Benefits of This Approach

1. **Clear Origin**: Every function name indicates which module it comes from
2. **No Name Conflicts**: PREFIX ensures unique function names across modules
3. **Selective Imports**: NOT IMPLEMENTED YET - each subordinate will eventually be a package
4. **Namespace Isolation**: No pollution of global namespace
5. **Easy Refactoring**: Moving functions between modules is straightforward
6. **Self-Documenting**: Code clearly shows dependencies

## Notes

- The `%^%` operator (matrix power) is available without prefix for convenience
- All internal helper functions also use the prefix convention
- When in doubt about a function's location, the prefix tells you exactly where to find it
