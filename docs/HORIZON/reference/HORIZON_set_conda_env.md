# Set the conda environment used by HORIZON for all external tools

Registers the path to a conda environment containing `samtools`,
`bedtools`, `macs3`, and `deeptools`. Must be called once per session
before any pipeline function that invokes these tools.

## Usage

``` r
HORIZON_set_conda_env(path)
```

## Arguments

- path:

  Character. Path to the conda environment directory (the directory that
  contains a `bin/` subdirectory).

## Value

The resolved path, invisibly.

## Details

Recommended setup (run once in a terminal):

    conda create -n horizon_cli -c bioconda -c conda-forge \
      samtools=1.23.1 bedtools=2.31.1 macs3=3.0.2 deeptools=3.5.5
