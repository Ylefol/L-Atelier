# Load a Previously Trained MOFA2 Model

Reloads a MOFA2 model saved by
[`HEPHAESTUS_run_mofa()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HEPHAESTUS_run_mofa.md)
(or by
[`MOFA2::run_mofa()`](https://rdrr.io/pkg/MOFA2/man/run_mofa.html)
directly) from its `.hdf5` file, without retraining. Training is slow
and stochastic, so keeping a loader separate from the run step matters
in practice – mirrors the run/parse split already used for HOMER
([`.homer_run()`](https://ylefol.github.io/L-Atelier/GAIA/reference/dot-homer_run.md)
vs
[`APOLLO_load_homer_results()`](https://ylefol.github.io/L-Atelier/GAIA/reference/APOLLO_load_homer_results.md)).

## Usage

``` r
HEPHAESTUS_load_mofa_model(file)
```

## Arguments

- file:

  Character. Path to a trained MOFA2 `.hdf5` model file.

## Value

A classed `hephaestus_mofa` list, identical in shape to the return value
of
[`HEPHAESTUS_run_mofa()`](https://ylefol.github.io/L-Atelier/GAIA/reference/HEPHAESTUS_run_mofa.md)
– `params` will be `NULL` since training options aren't recoverable from
the saved model file alone.

## Examples

``` r
if (FALSE) { # \dontrun{
result   <- HEPHAESTUS_run_mofa(X = list(...), outfile = "my_model.hdf5")
reloaded <- HEPHAESTUS_load_mofa_model("my_model.hdf5")

} # }
```
