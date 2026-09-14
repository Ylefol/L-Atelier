# Convert GAIA's samples x features list-of-matrices into MOFA2's long format

MOFA2's list-of-matrices input path (`create_mofa_from_matrix`) requires
identical, identically-ordered colnames (samples) across every view and
does not pad partial overlap. The long-format path
(`sample, feature, view, value[, group]`) has no such requirement – rows
for missing sample/view combinations are simply omitted, which is
exactly what's needed for GAIA's list-of-matrices inputs where views may
only partially overlap in samples.

## Usage

``` r
.hephaestus_mofa_to_long(X, groups = NULL)
```
