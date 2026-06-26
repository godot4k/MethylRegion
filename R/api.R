.methylregion_default_control <- function() {
  list(
    maxdist = 300,
    method = "pearson",
    maxseg = -1,
    mincpgs = 5,
    threads = 1,
    mode = 1,
    mtc = 1,
    name = "sample",
    trend = 0.6,
    minNo = -1,
    minFactor = 0.8,
    valley = 0.7,
    minMethDist = 0.1,
    randomseed = 26061981
  )
}

.methylregion_control <- function(controlist) {
  if (is.null(controlist)) {
    controlist <- list()
  }
  if (!is.list(controlist)) {
    stop("controlist must be a list.", call. = FALSE)
  }
  defaults <- .methylregion_default_control()
  for (name in names(controlist)) {
    defaults[[name]] <- controlist[[name]]
  }
  defaults
}

.methylregion_input <- function(input_dat, intput_dat = NULL) {
  if (missing(input_dat)) {
    if (is.null(intput_dat)) {
      stop("input_dat is required.", call. = FALSE)
    }
    input_dat <- intput_dat
  }
  input_dat <- as.data.frame(input_dat)
  required <- c("chr", "pos")
  missing_cols <- setdiff(required, names(input_dat))
  if (length(missing_cols) > 0) {
    stop(
      "input_dat must contain columns: ",
      paste(required, collapse = ", "),
      call. = FALSE
    )
  }
  input_dat
}

.methylregion_result <- function(result) {
  if (is.null(result) || (!is.data.frame(result) && length(result) == 0)) {
    return(data.frame())
  }
  result
}

#' Detect methylation regions for a binary phenotype
#'
#' `mr_bi()` is the main entry point for binary phenotype designs. The `data`
#' argument selects the study structure, while `cov.mod = NULL` versus a
#' non-`NULL` covariate data frame determines whether a covariate-adjusted
#' branch is used when that branch exists.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Data frame, matrix, or vector whose first column is the binary group
#'   or state indicator coded as 0 and 1.
#' @param data Study structure. One of `"independent"`, `"paired"`, or
#'   `"longitudinal"`.
#' @param cov.mod Optional data frame of sample-level covariates. For
#'   `data = "independent"`, `NULL` dispatches to `dmr_case_control()` and
#'   non-`NULL` dispatches to `dmr_case_control_cov()`.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
mr_bi <- function(input_dat, y, data = c("independent", "paired", "longitudinal"),
                  cov.mod = NULL, controlist = list(), intput_dat = NULL) {
  data <- match.arg(data)

  if (data == "independent") {
    if (is.null(cov.mod)) {
      return(dmr_case_control(input_dat, y, controlist = controlist, intput_dat = intput_dat))
    }
    return(dmr_case_control_cov(input_dat, y, cov.mod = cov.mod, controlist = controlist, intput_dat = intput_dat))
  }

  if (data == "paired") {
    if (!is.null(cov.mod)) {
      stop("Binary paired designs with cov.mod are not implemented yet.", call. = FALSE)
    }
    return(dmr_paired(input_dat, y, controlist = controlist, intput_dat = intput_dat))
  }

  stop("Binary longitudinal designs are not implemented yet.", call. = FALSE)
}

#' Detect methylation regions for a continuous phenotype
#'
#' `mr_continuous()` is the main entry point for continuous phenotype designs.
#' The `data` argument selects independent versus longitudinal/repeated-measure
#' analyses. `cov.mod` is passed through and may be `NULL`.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Data frame, matrix, or vector whose first column is a numeric
#'   continuous phenotype.
#' @param data Study structure. One of `"independent"` or `"longitudinal"`.
#' @param cov.mod Optional data frame of sample-level covariates.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
mr_continuous <- function(input_dat, y, data = c("independent", "longitudinal"),
                          cov.mod = NULL, controlist = list(), intput_dat = NULL) {
  data <- match.arg(data)

  if (data == "independent") {
    return(amr_continuous(input_dat, y, cov.mod = cov.mod, controlist = controlist, intput_dat = intput_dat))
  }

  amr_longitudinal(input_dat, y, cov.mod = cov.mod, controlist = controlist, intput_dat = intput_dat)
}

#' Detect DMRs for a binary case-control design without covariates
#'
#' `dmr_case_control()` compares two independent sample groups coded as 0 and 1.
#' It uses the 2D-KS/Mann-Whitney branch implementation.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Data frame, matrix, or vector whose first column is the binary group
#'   indicator. Controls must be coded 0 and cases/tests must be coded 1.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
dmr_case_control <- function(input_dat, y, controlist = list(), intput_dat = NULL) {
  input_dat <- .methylregion_input(input_dat, intput_dat)
  .methylregion_result(
    .dmr_case_control_impl(input_dat, y, cov.mod = NULL, controlist = .methylregion_control(controlist))
  )
}

#' Detect DMRs for a binary case-control design with covariates
#'
#' `dmr_case_control_cov()` tests the group effect while adjusting for
#' sample-level covariates using the covariate linear-model branch.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Binary group variable; controls are coded 0 and cases/tests are coded 1.
#' @param cov.mod Optional data frame of sample-level covariates.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
dmr_case_control_cov <- function(input_dat, y, cov.mod = NULL, controlist = list(), intput_dat = NULL) {
  input_dat <- .methylregion_input(input_dat, intput_dat)
  .methylregion_result(
    .dmr_case_control_cov_impl(input_dat, y, cov.mod = cov.mod, controlist = .methylregion_control(controlist))
  )
}

#' Detect paired pre/post treatment DMRs
#'
#' `dmr_paired()` compares paired binary states coded as 0 and 1. If `y`
#' contains a `pair_id` column, the region p value is calculated with a paired
#' Wilcoxon test; otherwise it falls back to the unpaired rank-sum test.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Data frame whose first column is the binary pre/post state. Include
#'   `pair_id` for paired testing.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
dmr_paired <- function(input_dat, y, controlist = list(), intput_dat = NULL) {
  input_dat <- .methylregion_input(input_dat, intput_dat)
  .methylregion_result(
    .dmr_paired_impl(input_dat, y, cov.mod = NULL, controlist = .methylregion_control(controlist))
  )
}

#' Detect AMRs associated with a continuous phenotype
#'
#' `amr_continuous()` identifies methylation regions associated with a numeric
#' phenotype, optionally adjusting for sample-level covariates.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Numeric continuous phenotype in the first column.
#' @param cov.mod Optional data frame of sample-level covariates.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
amr_continuous <- function(input_dat, y, cov.mod = NULL, controlist = list(), intput_dat = NULL) {
  input_dat <- .methylregion_input(input_dat, intput_dat)
  .methylregion_result(
    .amr_continuous_impl(input_dat, y, cov.mod = cov.mod, controlist = .methylregion_control(controlist))
  )
}

#' Detect longitudinal AMRs associated with a continuous phenotype
#'
#' `amr_longitudinal()` identifies regions associated with a numeric continuous
#' phenotype while allowing adjustment variables and repeated-measure identifiers.
#' If a subject/id column is available and `lmerTest` is installed, the branch
#' implementation attempts a random-intercept mixed model; otherwise it falls
#' back to a fixed-effect linear model. Discrete phenotypes are not supported.
#'
#' @param input_dat Data frame with `chr`, `pos`, and one methylation column per sample.
#' @param y Data frame whose first column is a numeric continuous phenotype.
#'   Additional columns are adjustment variables or repeated-measure identifiers.
#' @param cov.mod Optional data frame of additional sample-level covariates.
#' @param controlist Optional list of segmentation controls.
#' @param intput_dat Deprecated spelling retained for compatibility.
#' @return A data frame of candidate methylation regions.
amr_longitudinal <- function(input_dat, y, cov.mod = NULL, controlist = list(), intput_dat = NULL) {
  input_dat <- .methylregion_input(input_dat, intput_dat)
  .methylregion_result(
    .amr_longitudinal_impl(input_dat, y, cov.mod = cov.mod, controlist = .methylregion_control(controlist))
  )
}
