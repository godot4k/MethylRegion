# Bulk TestData Multi-Group `mr_bi()` Check

This records the script used to test multi-group pairwise `mr_bi()` on
`tests/TestData/bulk.sub.txt.20.Rds`, plus the observed summary and first 10
result rows.

## Script

```r
.libPaths(c("/tmp/methylregion-r-lib", .libPaths()))
library(MethylRegion)

dat <- readRDS("tests/TestData/bulk.sub.txt.20.Rds")
y <- data.frame(group = rep(c("A", "B", "C"), c(7, 7, 6)))

set.seed(1)
res <- mr_bi(dat, y, data.type = "independent")

cat("comparisons:", paste(attr(res, "comparisons"), collapse = ", "), "\n")
cat("dim:", paste(dim(res), collapse = " x "), "\n")
print(table(res$comparison))

sig <- subset(res, FDR < 0.05 & abs(mean_diff) >= 0.1)
cat("significant_FDR_lt_0.05_absdiff_ge_0.1:", nrow(sig), "\n")
print(table(sig$comparison))

head(res[order(res$FDR, res$p_value), ], 10)
```

## Summary Output

```text
comparisons: A_vs_B, A_vs_C, B_vs_C
dim: 827 x 15

A_vs_B A_vs_C B_vs_C
   288    271    268

significant_FDR_lt_0.05_absdiff_ge_0.1: 22

A_vs_B A_vs_C B_vs_C
     8      9      5
```

## First 10 Result Rows

Rows are ordered by `FDR`, then `p_value`.

| comparison | groupA | groupB | chr | start | end | N.CpGs | ks_stat | mean_diff | p_value | methX | methY | e_value | FDR | e_adjust |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| A_vs_B | A | B | chr21 | 46902889 | 46902982 | 8 | 0.6160714 | -0.116247115 | 2.001988e-11 | 0.065137175 | 0.5000000 | 3.735190e+14 | 5.765725e-09 | 1.556329e+13 |
| A_vs_C | A | C | chr21 | 38070968 | 38071211 | 31 | 0.2565284 | -0.009417608 | 1.128579e-07 | 0.013237281 | 0.4615385 | 8.218612e+02 | 3.058450e-05 | 1.304060e+02 |
| A_vs_B | A | B | chr21 | 36258812 | 36259286 | 18 | 0.3134921 | 0.108920584 | 3.069506e-07 | 0.165617804 | 0.5000000 | 2.532848e+02 | 4.420089e-05 | 5.892390e+01 |
| A_vs_B | A | B | chr21 | 47401831 | 47401992 | 10 | 0.4142857 | 0.048359949 | 4.778519e-07 | 0.041675305 | 0.5000000 | 1.586203e+09 | 4.587378e-05 | 1.046454e+08 |
| A_vs_C | A | C | chr21 | 15352003 | 15352137 | 13 | 0.4532967 | -0.141409599 | 3.457615e-07 | 0.292349690 | 0.4615385 | 3.370361e+02 | 4.685069e-05 | 5.969643e+01 |
| A_vs_C | A | C | chr21 | 47062451 | 47062652 | 16 | 0.2931548 | -0.018271833 | 2.522667e-06 | 0.015514598 | 0.4615385 | 6.067383e+04 | 2.278809e-04 | 4.701662e+03 |
| B_vs_C | B | C | chr21 | 46677187 | 46677379 | 9 | 0.4365079 | 0.219347988 | 9.769686e-07 | 0.393591945 | 0.4615385 | 6.033667e+04 | 2.618276e-04 | 3.827326e+03 |
| A_vs_B | A | B | chr21 | 33104622 | 33104719 | 14 | 0.2244898 | -0.003993775 | 6.480199e-06 | 0.002599976 | 0.5000000 | 2.144933e+10 | 4.548230e-04 | 1.266106e+09 |
| A_vs_B | A | B | chr21 | 35831957 | 35832096 | 10 | 0.3642857 | 0.081631731 | 7.896233e-06 | 0.092660533 | 0.5000000 | 1.526039e+05 | 4.548230e-04 | 1.589624e+04 |
| A_vs_B | A | B | chr21 | 9826610 | 9826644 | 10 | 0.3428571 | 0.036299707 | 8.133278e-05 | 0.184916941 | 0.5000000 | 6.699443e+02 | 3.903974e-03 | 1.307468e+02 |

Note: the `A/B/C` labels are simulated from the 20 sample columns in order
(`A = S1-S7`, `B = S8-S14`, `C = S15-S20`), so this check validates the
multi-group pairwise behavior rather than a biological grouping.
