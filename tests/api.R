library(MethylRegion)

stopifnot(is.function(dmr_case_control))
stopifnot(is.function(dmr_case_control_cov))
stopifnot(is.function(dmr_paired))
stopifnot(is.function(dmr_longitudinal))
stopifnot(is.function(amr_continuous))
stopifnot(is.function(amr_longitudinal))
stopifnot(is.function(mr_bi))
stopifnot(is.function(mr_continuous))

bad_dat <- data.frame(sample1 = 0.1, sample2 = 0.2)
err <- try(dmr_case_control(bad_dat, data.frame(group = c(0, 1))), silent = TRUE)
stopifnot(inherits(err, "try-error"))

dat <- data.frame(
  chr = rep("chr1", 8),
  pos = seq(100, 800, by = 100),
  sample1 = c(0.1, 0.1, 0.2, 0.2, 0.8, 0.8, 0.7, 0.7),
  sample2 = c(0.2, 0.1, 0.2, 0.3, 0.7, 0.8, 0.8, 0.7),
  sample3 = c(0.8, 0.8, 0.7, 0.7, 0.2, 0.1, 0.2, 0.1),
  sample4 = c(0.7, 0.8, 0.8, 0.7, 0.1, 0.2, 0.2, 0.1)
)
control <- list(mincpgs = 20)

binary_y <- data.frame(group = c(0, 0, 1, 1))
paired_y <- data.frame(group = c(0, 1, 0, 1), pair_id = c("p1", "p1", "p2", "p2"))
binary_long_y <- data.frame(group = c(0, 1, 0, 1), family = c("f1", "f1", "f2", "f2"))
continuous_y <- data.frame(phenotype = c(1, 2, 4, 5))
long_y <- data.frame(
  phenotype = c(1.0, 1.6, 3.8, 4.4),
  subject_id = c("s1", "s1", "s2", "s2"),
  visit = c("pre", "post", "pre", "post")
)
covariates <- data.frame(age = c(40, 42, 57, 59))

results <- list(
  dmr_case_control = dmr_case_control(dat, binary_y, controlist = control),
  dmr_case_control_cov = dmr_case_control_cov(dat, binary_y, cov.mod = covariates, controlist = control),
  dmr_paired = dmr_paired(dat, paired_y, controlist = control),
  dmr_longitudinal = dmr_longitudinal(dat, binary_long_y, controlist = control),
  amr_continuous = amr_continuous(dat, continuous_y, controlist = control),
  amr_longitudinal = amr_longitudinal(dat, long_y, controlist = control),
  mr_bi_independent = mr_bi(dat, binary_y, data.type = "independent", controlist = control),
  mr_bi_independent_cov = mr_bi(dat, binary_y, data.type = "independent", cov.mod = covariates, controlist = control),
  mr_bi_paired = mr_bi(dat, paired_y, data.type = "paired", controlist = control),
  mr_bi_paired_cov = mr_bi(dat, paired_y, data.type = "paired", cov.mod = covariates, controlist = control),
  mr_bi_longitudinal = mr_bi(dat, binary_long_y, data.type = "longitudinal", controlist = control),
  mr_continuous_independent = mr_continuous(dat, continuous_y, data.type = "independent", controlist = control),
  mr_continuous_longitudinal = mr_continuous(dat, long_y, data.type = "longitudinal", controlist = control)
)

stopifnot(all(vapply(results, function(x) is.null(x) || is.data.frame(x), logical(1))))

e_value_methods <- c("dmr_case_control", "mr_bi_independent")
for (name in setdiff(names(results), e_value_methods)) {
  stopifnot(!any(c("e_value", "e_adjust", "e_bh_significant") %in% names(results[[name]])))
}

mixed_control <- list(mincpgs = 2, trend = 0)
mixed_res <- mr_bi(dat, binary_long_y, data.type = "longitudinal", controlist = mixed_control)
stopifnot(is.data.frame(mixed_res))
if (nrow(mixed_res) > 0) {
  stopifnot(all(c("coef_lmm", "p_value", "FDR") %in% names(mixed_res)))
  stopifnot(!any(c("e_value", "e_adjust", "e_bh_significant") %in% names(mixed_res)))
}

multi_dat <- data.frame(
  chr = rep("chr1", 12),
  pos = seq(100, 1200, by = 100),
  sample1 = rep(0.05, 12),
  sample2 = rep(0.06, 12),
  sample3 = rep(0.95, 12),
  sample4 = rep(0.94, 12),
  sample5 = rep(c(0.4, 0.6), 6),
  sample6 = rep(c(0.5, 0.7), 6)
)
multi_y <- data.frame(group = c("A", "A", "B", "B", "C", "C"))
multi_res <- mr_bi(
  multi_dat,
  multi_y,
  data.type = "independent",
  controlist = list(mincpgs = 2, trend = 0, valley = 0)
)
stopifnot(is.data.frame(multi_res))
stopifnot(all(c("comparison", "groupA", "groupB") %in% names(multi_res)))
stopifnot(identical(attr(multi_res, "comparisons"), c("A_vs_B", "A_vs_C", "B_vs_C")))
stopifnot(nrow(multi_res) == 3)
stopifnot(identical(as.integer(table(multi_res$comparison)[attr(multi_res, "comparisons")]), c(1L, 1L, 1L)))
