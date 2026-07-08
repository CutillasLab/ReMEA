README
================
Luke Higgins
2026-07-08

# ReMEA

## Response Marker Enrichment Analysis

`ReMEA` is an enrichment analysis method, whereby protein perturbation
response signatures are used to score perturbations according to the
inferred sensitivity profile.

## Installation

### R console

#### One-liner

Run: `devtools::install_github("CutillasLab/ReMEA")`

- You may first need to install `devtools`:
  `install.packages("devtools")`

### Manual

1.  Download the Package Archive File (ReMEAx.x.x.tar.gz) of the latest
    release.
    - Not the Source code
2.  Run: devtools::install_local(path = “C:/path/to/ReMEAx.x.x.tar.gz”)
    - Replace the string argument to path with the actual location
    - You may first need to install `devtools`:
      `install.packages("devtools")`

## Usage

For a complete analysis of your data, the
`ReMEA::complete_ReMEA_analysis` function runs the analysis for all
perturbation modalities and signatures. Default settings apply the
analysis using the “CellLines” signatures. This returns a list of all
the resulting scores. This includes the individual db scores and main
ReMEA score (av.effect).

``` r
results <- ReMEA::complete_ReMEA_analysis(protein_data = proteomics_data)
```

We to analyse data using signatures derived from the beat AML data,
enter select the “BeatAML” signatures in the function. In addition, the
haem tumour type must be specified:

``` r
results <- ReMEA::complete_ReMEA_analysis(protein_data = proteomics_data,
                                          signature_version = "BeatAML",
                                          tumour_type = c("haem"))
```

The core component of ReMEA is the
`ReMEA::response_marker_enrichment_analysis` function. If you wish use
use specific signature sets, these can be specified here:

``` r
results <- ReMEA::response_marker_enrichment_analysis(protein_data = proteomics_data,
                                                      marker_type = "DRUG",
                                                      tumour_type = "pan",
                                                      signature_version = "CellLines")
```

If you wish to use your own signatures, this can be achieved using the
base `ReMEA::signature_enrichment` function.

## Interpretation of results

The primary score for interpretation is the av.effect score, this is the
principal ReMEA score which incorporates the output from the different
signature sets (constructed from the different cell line proteomics
datasets).

A positive ReMEA score (av.effect) is interpreted as an increase in
sensitivity to that perturbation. A negative score suggests resistance
to that perturbation. The magnitude of scores can differ between
datasets (similar to other enrichment methods). Accompanying p-values
and overall ranking of scores should be used to interpret the returned
scores.

## Interpretation of results

For reference, an explanation of the fields from the output can be
reviewd using the `ReMEA::output_description` function.

``` r
ReMEA::output_description(
  combined_scores = TRUE,
  individual_scores = FALSE
)
```

    ## 
    ## Combined scores
    ## +---------------------+-----------+--------------------------------------------+
    ## | Column              | Type      | Description                                |
    ## +---------------------+-----------+--------------------------------------------+
    ## | perturbagen         | character | The perturbation in question.              |
    ## | av.effect           | numeric   | The composite score across signature       |
    ## |                     |           | selections.                                |
    ## | combined_ks_pvalue  | numeric   | The combined Kolmogorov-Smirnov p-value    |
    ## |                     |           | from the different signature selections.   |
    ## |                     |           | P-values are combined using Stouffer's     |
    ## |                     |           | method.                                    |
    ## | combined_bws_pvalue | numeric   | The combined BWS p-value from the          |
    ## |                     |           | different signature selections. P-values   |
    ## |                     |           | are combined using Stouffer's method.      |
    ## | alpha               | numeric   | -log10(combined_ks_pvalue) * av.effect.    |
    ## | ks_qvalue           | numeric   | Benjamini-Hochberg corrected               |
    ## |                     |           | combined_ks_pvalue.                        |
    ## | ks_padj_bonferroni  | numeric   | Bonferroni corrected combined_ks_pvalue.   |
    ## | bws_qvalue          | numeric   | Benjamini-Hochberg corrected               |
    ## |                     |           | combined_bws_pvalue.                       |
    ## | bws_padj_bonferroni | numeric   | Bonferroni corrected combined_bws_pvalue.  |
    ## +---------------------+-----------+--------------------------------------------+

To further investigate a particular ReMEA score, the individual dataset
signature metrics can be reviewed. See below for a description of the
data fields.

``` r
ReMEA::output_description(
  combined_scores = FALSE,
  individual_scores = TRUE
)
```

    ## 
    ## Individual database scores
    ## +---------------------------+-----------+--------------------------------------+
    ## | Column                    | Type      | Description                          |
    ## +---------------------------+-----------+--------------------------------------+
    ## | perturbagen               | character | The perturbation in question. For    |
    ## |                           |           | compounds, this includes the         |
    ## |                           |           | site-specific code from the GDSC2    |
    ## |                           |           | drug annotations.                    |
    ## | zscore.resistance         | numeric   | The z-score for the resistance       |
    ## |                           |           | signature.                           |
    ## | zscore.sensitivity        | numeric   | The z-score for the sensitivity      |
    ## |                           |           | signature.                           |
    ## | pvalue.ks                 | numeric   | The p-value from the                 |
    ## |                           |           | Kolmogorov-Smirnov test comparing    |
    ## |                           |           | the resistance and sensitivity       |
    ## |                           |           | signatures.                          |
    ## | stat.ks                   | numeric   | The Kolmogorov-Smirnov test          |
    ## |                           |           | statistic comparing the sensitivity  |
    ## |                           |           | and resistance signatures.           |
    ## | pvalue.bws                | numeric   | The BWS p-value comparing the        |
    ## |                           |           | resistance and sensitivity           |
    ## |                           |           | signatures.                          |
    ## | resistance.markers        | character | Proteins used in the resistance      |
    ## |                           |           | signature.                           |
    ## | sensitivity.markers       | character | Proteins used in the sensitivity     |
    ## |                           |           | signature.                           |
    ## | av.rank.resistance        | numeric   | Average rank of the resistance       |
    ## |                           |           | signature.                           |
    ## | av.rank.sensitivity       | numeric   | Average rank of the sensitivity      |
    ## |                           |           | signature.                           |
    ## | geommean.rank.resistance  | numeric   | Geometric mean of the rank of the    |
    ## |                           |           | resistance signature.                |
    ## | geommean.rank.sensitivity | numeric   | Geometric mean of the rank of the    |
    ## |                           |           | sensitivity signature.               |
    ## | med.rank.resistance       | numeric   | Median rank of the resistance        |
    ## |                           |           | signature.                           |
    ## | med.rank.sensitivity      | numeric   | Median rank of the sensitivity       |
    ## |                           |           | signature.                           |
    ## | n.resistance.markers      | integer   | Number of resistance markers.        |
    ## | n.sensitivity.markers     | integer   | Number of sensitivity markers.       |
    ## | signature.type            | character | Origin of the signature, structured  |
    ## |                           |           | as dataset_cancertype_modality.      |
    ## | delta.zscore              | numeric   | zscore.sensitivity -                 |
    ## |                           |           | zscore.resistance.                   |
    ## | total_counts              | integer   | n.resistance.markers +               |
    ## |                           |           | n.sensitivity.markers.               |
    ## | delta_rank_mean           | numeric   | log2(av.rank.sensitivity) -          |
    ## |                           |           | log2(av.rank.resistance).            |
    ## | delta_rank_median         | numeric   | log2(med.rank.sensitivity) -         |
    ## |                           |           | log2(med.rank.resistance).           |
    ## | delta_rank_geomean        | numeric   | log2(geommean.rank.sensitivity) -    |
    ## |                           |           | log2(geommean.rank.resistance).      |
    ## | max_delta_score           | numeric   | Maximum absolute delta score across  |
    ## |                           |           | delta.zscore, delta_rank_geomean,    |
    ## |                           |           | delta_rank_mean, and                 |
    ## |                           |           | delta_rank_median.                   |
    ## | ks_qvalue                 | numeric   | Benjamini-Hochberg corrected         |
    ## |                           |           | Kolmogorov-Smirnov p-value.          |
    ## |                           |           | Adjusted within each                 |
    ## |                           |           | signature.type.                      |
    ## | ks_p_padj_bonf            | numeric   | Bonferroni corrected                 |
    ## |                           |           | Kolmogorov-Smirnov p-value.          |
    ## |                           |           | Adjusted within each                 |
    ## |                           |           | signature.type.                      |
    ## | bws_qvalue                | numeric   | Benjamini-Hochberg corrected BWS     |
    ## |                           |           | p-value. Adjusted within each        |
    ## |                           |           | signature.type.                      |
    ## | bws_padj_bonf             | numeric   | Bonferroni corrected BWS p-value.    |
    ## |                           |           | Adjusted within each                 |
    ## |                           |           | signature.type.                      |
    ## +---------------------------+-----------+--------------------------------------+
