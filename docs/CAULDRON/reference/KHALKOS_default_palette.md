# Default CAULDRON colour palettes

Returns a vector of `n` colours from the CAULDRON palette suite.
Categorical palettes are based on a hand-curated 20-colour set extended
via `colorRampPalette` for `n > 20`. Sequential and diverging palettes
are interpolated between key anchor colours.

## Usage

``` r
KHALKOS_default_palette(n, type = c("categorical", "sequential", "diverging"))
```

## Arguments

- n:

  Integer. Number of colours to return.

- type:

  Character. Palette type: `"categorical"` (default), `"sequential"`, or
  `"diverging"`.

## Value

Character vector of `n` hex colour codes.
