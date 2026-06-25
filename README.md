# MethylRegion

**MethylRegion: A unified framework for differential methylation regions across designs, traits, and covariates**

MethylRegion consolidates the AMRfinder development branches into one R package
with a small, design-oriented API. The package separates DMR functions for
discrete group comparisons from AMR functions for continuous phenotype
associations.

## Functions

| Function | Design | Trait | Covariates | Source branch |
| --- | --- | --- | --- | --- |
| `dmr_case_control()` | Case-control | Binary/discrete | No | `dmr.no.cov.2dks-mwu` |
| `dmr_case_control_cov()` | Case-control | Binary/discrete | Yes | `dmr.with.cov.lm` |
| `dmr_paired()` | Pre/post treatment or paired case-control | Binary/discrete | No | `dmr.no.cov.2dks-paired-wilcox` |
| `amr_continuous()` | EWAS / continuous phenotype association | Continuous | Yes/No | `main` |
| `amr_longitudinal()` | Longitudinal or repeated-measure AMR | Continuous | Yes/No | `dmr.with.cov.lm.phenotype-model` |

`amr_longitudinal()` currently supports continuous phenotypes only. It does not
support discrete or binary phenotypes.

## Installation

```r
remotes::install_github("godot4k/MethylRegion")
```

For local development:

```sh
R CMD INSTALL MethylRegion
```

## Input

The methylation matrix is supplied as a data frame with genomic columns followed
by one methylation column per sample:

```r
input_dat <- data.frame(
  chr = rep("chr1", 6),
  pos = seq(100, 600, by = 100),
  sample1 = c(0.1, 0.2, 0.2, 0.7, 0.8, 0.8),
  sample2 = c(0.2, 0.2, 0.3, 0.6, 0.7, 0.7),
  sample3 = c(0.8, 0.7, 0.7, 0.2, 0.2, 0.1),
  sample4 = c(0.7, 0.8, 0.8, 0.1, 0.2, 0.2)
)
```

## Examples

Independent case-control DMR analysis:

```r
library(MethylRegion)

y <- data.frame(group = c(0, 0, 1, 1))

dmr <- dmr_case_control(
  input_dat,
  y,
  controlist = list(mincpgs = 2, trend = 0)
)
```

Case-control DMR analysis with covariates:

```r
covariates <- data.frame(age = c(43, 51, 39, 58))

dmr_cov <- dmr_case_control_cov(
  input_dat,
  y,
  cov.mod = covariates,
  controlist = list(mincpgs = 2, trend = 0)
)
```

Paired pre/post DMR analysis:

```r
y_paired <- data.frame(
  state = c(0, 1, 0, 1),
  pair_id = c("p1", "p1", "p2", "p2")
)

dmr_prepost <- dmr_paired(
  input_dat,
  y_paired,
  controlist = list(mincpgs = 2, trend = 0)
)
```

Continuous phenotype AMR analysis:

```r
y_continuous <- data.frame(phenotype = c(1.2, 2.0, 4.2, 5.1))

amr <- amr_continuous(
  input_dat,
  y_continuous,
  controlist = list(mincpgs = 2, trend = 0)
)
```

Longitudinal continuous-phenotype AMR analysis:

```r
y_long <- data.frame(
  phenotype = c(1.2, 1.8, 3.9, 4.5),
  subject_id = c("s1", "s1", "s2", "s2"),
  visit = c("pre", "post", "pre", "post")
)

amr_long <- amr_longitudinal(
  input_dat,
  y_long,
  controlist = list(mincpgs = 2, trend = 0)
)
```

## Development notes

The five exported functions intentionally preserve branch-specific statistics
and output columns. This keeps the original prototype behavior available while
placing the implementations behind one package namespace.
