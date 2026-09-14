# CYAN: Niche Genomic Analyses

CYAN handles specialised genomic analyses — RNA editing detection via
REDItools2 and somatic variant calling via GATK MuTect2 — positioned
between upstream BAM processing (HORIZON) and expression-level analysis
(GAIA/CAULDRON).

## Details

Named after CYAN from Horizon Zero Dawn: a focused, independent
analytical unit rather than a large orchestrating system.

## Package design

- All public functions are prefixed `CYAN_`

- One sample per call; cohort looping is left to the user

- Large intermediates written to `work_dir`; final outputs to
  `output_dir`

- Python tools (REDItools2) run through a managed basilisk environment

- Java tools (GATK/MuTect2) invoked via `system2("gatk", ...)`

## See also

Useful links:

- <https://ylefol.github.io/L-Atelier/CYAN/>

- <https://github.com/Ylefol/L-Atelier/tree/master/CYAN>

- Report bugs at <https://github.com/Ylefol/L-Atelier/issues>

## Author

**Maintainer**: Yohan Lefol <yohanlefol@proton.me>

Authors:

- Yohan Lefol <yohanlefol@proton.me>
