# Convert Consensus Peaks to SAF Format

Converts a consensus peak data.frame to SAF (Simplified Annotation
Format) for use with Rsubread::featureCounts().

## Usage

``` r
ELEUTHIA_peaks_to_saf(consensus_peaks)
```

## Arguments

- consensus_peaks:

  A consensus peak data.frame from ELEUTHIA_create_consensus_peaks().

## Value

A data.frame in SAF format with columns: GeneID, Chr, Start, End,
Strand.

## Details

SAF format is required by Rsubread::featureCounts() for counting reads
in custom genomic regions. The "GeneID" column contains the peak_id.

## Examples

``` r
if (FALSE) { # \dontrun{
consensus <- ELEUTHIA_create_consensus_peaks(peak_list)
saf <- ELEUTHIA_peaks_to_saf(consensus)
counts <- Rsubread::featureCounts(bam_files, annot.ext = saf, isGTFAnnotationFile = FALSE)

} # }
```
