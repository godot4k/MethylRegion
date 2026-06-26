\name{dmr_case_control}
\alias{mr_bi}
\alias{mr_continuous}
\alias{dmr_case_control}
\alias{dmr_case_control_cov}
\alias{dmr_paired}
\alias{dmr_longitudinal}
\alias{amr_continuous}
\alias{amr_longitudinal}
\title{MethylRegion study-design functions}
\description{
MethylRegion provides separate entry points for methylation-region analyses
across binary DMR designs and continuous-phenotype AMR designs.
}
\usage{
mr_bi(input_dat, y, data = c("independent", "paired", "longitudinal"),
  cov.mod = NULL, controlist = list(), intput_dat = NULL)

mr_continuous(input_dat, y, data = c("independent", "longitudinal"),
  cov.mod = NULL, controlist = list(), intput_dat = NULL)

dmr_case_control(input_dat, y, controlist = list(), intput_dat = NULL)

dmr_case_control_cov(input_dat, y, cov.mod = NULL, controlist = list(),
  intput_dat = NULL)

dmr_paired(input_dat, y, controlist = list(), intput_dat = NULL)

dmr_longitudinal(input_dat, y, cov.mod = NULL, controlist = list(),
  intput_dat = NULL)

amr_continuous(input_dat, y, cov.mod = NULL, controlist = list(),
  intput_dat = NULL)

amr_longitudinal(input_dat, y, cov.mod = NULL, controlist = list(),
  intput_dat = NULL)
}
\arguments{
\item{input_dat}{A data frame with columns \code{chr}, \code{pos}, and one
methylation column per sample.}
\item{y}{Study-design response data. For DMR functions, the first column is a
binary group or state coded 0 and 1. For AMR functions, the first column is a
numeric continuous phenotype.}
\item{data}{Study structure for the main entry points. For \code{mr_bi()},
choose \code{"independent"}, \code{"paired"}, or \code{"longitudinal"}.
For \code{mr_continuous()}, choose \code{"independent"} or
\code{"longitudinal"}.}
\item{controlist}{Optional list of segmentation controls. Missing entries are
filled from MethylRegion defaults.}
\item{intput_dat}{Deprecated misspelling retained for compatibility with the
original AMRfinder prototypes.}
\item{cov.mod}{Optional data frame of sample-level covariates. Main entry
points use \code{is.null(cov.mod)} to select no-covariate versus
covariate-aware branches where both branches exist.}
}
\details{
\code{mr_bi()} is the main entry point for binary phenotypes. With
\code{data = "independent"}, \code{cov.mod = NULL} dispatches to
\code{dmr_case_control()} and non-\code{NULL} \code{cov.mod} dispatches to
\code{dmr_case_control_cov()}. With \code{data = "paired"},
\code{cov.mod = NULL} dispatches to \code{dmr_paired()}; paired binary designs
with covariates dispatch to \code{dmr_longitudinal()}. With
\code{data = "longitudinal"}, \code{mr_bi()} dispatches to
\code{dmr_longitudinal()} whether or not \code{cov.mod} is supplied.

\code{mr_continuous()} is the main entry point for continuous phenotypes. With
\code{data = "independent"}, it dispatches to \code{amr_continuous()}. With
\code{data = "longitudinal"}, it dispatches to \code{amr_longitudinal()}.
\code{cov.mod} may be \code{NULL} or a covariate data frame in both continuous
branches.

\code{dmr_case_control()} handles independent binary case-control comparisons
without covariates using the 2D-KS/Mann-Whitney branch.

\code{dmr_case_control_cov()} handles binary case-control comparisons with
sample-level covariates using the covariate linear-model branch.

\code{dmr_paired()} handles pre/post or paired binary designs. If \code{y}
contains \code{pair_id}, the region p value is calculated with a paired
Wilcoxon test.

\code{dmr_longitudinal()} handles binary longitudinal or family-structured
designs using a linear mixed model. The region-level model is
\code{meth ~ phenotype + (1 | family)} when a family-like random-effect column
is available. Additional non-random columns in \code{y} or \code{cov.mod} are
included as fixed effects, matching the model-building behavior of
\code{amr_longitudinal()}. If mixed-model fitting is unavailable, it falls back
to a fixed-effect linear model.

\code{amr_continuous()} handles continuous phenotype-associated methylation
regions.

\code{amr_longitudinal()} handles longitudinal or repeated-measure AMRs for
continuous phenotypes only. If an id-like column is present and \pkg{lmerTest}
is installed, the implementation attempts a random-intercept mixed model;
otherwise it uses a fixed-effect linear model.
}
\value{
A data frame of candidate methylation regions. Output columns vary by design
because each function preserves the corresponding branch-specific statistics.
}
\examples{
dat <- data.frame(
  chr = rep("chr1", 6),
  pos = seq(100, 600, by = 100),
  sample1 = c(0.1, 0.2, 0.2, 0.7, 0.8, 0.8),
  sample2 = c(0.2, 0.2, 0.3, 0.6, 0.7, 0.7),
  sample3 = c(0.8, 0.7, 0.7, 0.2, 0.2, 0.1),
  sample4 = c(0.7, 0.8, 0.8, 0.1, 0.2, 0.2)
)
y_bin <- data.frame(group = c(0, 0, 1, 1))
res <- mr_bi(dat, y_bin, data = "independent")

y_cont <- data.frame(phenotype = c(1, 2, 4, 5))
res2 <- mr_continuous(dat, y_cont, data = "independent")
}
\keyword{models}
\keyword{htest}
