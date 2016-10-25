\name{QDM}
\alias{QDM}
\title{
Univariate bias correction via quantile delta mapping
}
\description{
Univariate bias correction based on the quantile delta mapping \code{QDM}
version of the quantile mapping algorithm from Cannon et al. (2015).
\code{QDM} constrains model-projected changes in quantiles to
be preserved following bias correction by quantile mapping.
}
\usage{
QDM(o.c, m.c, m.p, ratio=FALSE, trace=0.05, trace.calc=0.5*trace,
    jitter.factor=0, n.tau=NULL, ratio.max=2,
    ratio.max.trace=10*trace, ECBC=FALSE, ties='first',
    subsample=NULL)
}
\arguments{
  \item{o.c}{
      vector of observed samples during the calibration period.
  }
  \item{m.c}{
      vector of model outputs during the calibration period.
  }
  \item{m.p}{
      vector of model outputs during the projected period.
  }
  \item{ratio}{
      logical value indicating if samples are of a ratio quantity (e.g., precipitation).
  }
  \item{trace}{
      numeric value indicating the threshold below which values of a ratio quantity (e.g., \code{ratio=TRUE}) should be considered exact zeros.
  }
  \item{trace.calc}{
      numeric value of a threshold used internally when handling of exact zeros; defaults to one half of \code{trace}.
  }
  \item{jitter.factor}{
      optional strength of jittering to be applied when quantities are quantized.
  }
  \item{n.tau}{
      number of quantiles used in the quantile mapping; \code{NULL} equals the length of the \code{m.p} series.
  }
  \item{ratio.max}{
      numeric value indicating the maximum proportional change allowed for ratio quantities below the \code{ratio.max.trace} threshold.
  }
  \item{ratio.max.trace}{
      numeric value of a trace threshold used to constrain the proportional change in ratio quantities to \code{ratio.max}; defaults to ten times \code{trace}.
  }
  \item{ECBC}{
      logical value indicating whether \code{mhat.p} outputs should be ordered according to \code{o.c} ranks, i.e., as in the empirical copula-bias correction (ECBC) algorithm.
  }
  \item{ties}{
      method used to handle ties when calculating ordinal ranks.
  }
  \item{subsample}{
      use \code{subsample} draws of size \code{n.tau} to calculate empirical
      quantiles; if \code{NULL}, calculate normally.
  }
}
\value{
     a list of with elements consisting of:
      \item{mhat.c}{vector of bias corrected \code{m.c} values for the calibration period.}
      \item{mhat.p}{vector of bias corrected \code{m.p} values for the projection period.}
}

\seealso{
\code{\link{MBCp}, \link{MBCr}, \link{MRS}, \link{escore}}
}
\references{
Cannon, A.J., S.R. Sobie, and T.Q. Murdock, 2015. Bias correction
of simulated precipitation by quantile mapping: How well do methods
preserve relative changes in quantiles and extremes? Journal of
Climate, 28:6938-6959. doi:10.1175/JCLI-D-14-00754.1
}
