\name{MBCp}
\alias{MBCp}
\title{
Multivariate bias correction (Pearson correlation)
}
\description{
Multivariate bias correction that matches marginal
distributions using \code{\link{QDM}} and the Pearson correlation
dependence structure following Cannon (2016).
}
\usage{
MBCp(o.c, m.c, m.p, iter=20, cor.thresh=1e-4,
     ratio.seq=rep(FALSE, ncol(o.c)), trace=0.05,
     trace.calc=0.5*trace, jitter.factor=0, n.tau=NULL,
     ratio.max=2, ratio.max.trace=10*trace, ties='first',
     qmap.precalc=FALSE, silent=FALSE)
}
\arguments{
  \item{o.c}{
      matrix of observed samples during the calibration period.
  }
  \item{m.c}{
      matrix of model outputs during the calibration period.
  }
  \item{m.p}{
      matrix of model outputs during the projected period.
  }
  \item{iter}{
      maximum number of algorithm iterations.
  }
  \item{cor.thresh}{
      if greater than zero, a threshold indicating the change in magnitude of Pearson correlations required for convergence.
  }
  \item{ratio.seq}{
      vector of logical values indicating if samples are of a ratio quantity (e.g., precipitation).
  }
  \item{trace}{
      numeric values indicating thresholds below which values of a ratio quantity (e.g., \code{ratio=TRUE}) should be considered exact zeros.
  }
  \item{trace.calc}{
      numeric values of thresholds used internally when handling of exact zeros; defaults to one half of \code{trace}.
  }
  \item{jitter.factor}{
      optional strength of jittering to be applied when quantities are quantized.
  }
  \item{n.tau}{
      number of quantiles used in the quantile mapping; \code{NULL} equals the length of the \code{m.p} series.
  }
  \item{ratio.max}{
      numeric values indicating the maximum proportional changes allowed for ratio quantities below the \code{ratio.max.trace} threshold.
  }
  \item{ratio.max.trace}{
      numeric values of trace thresholds used to constrain the proportional change in ratio quantities to \code{ratio.max}; defaults to ten times \code{trace}.
  }
  \item{ties}{
      method used to handle ties when calculating ordinal ranks.
  }
  \item{qmap.precalc}{
      logical value indicating if \code{m.c} and \code{m.p} are outputs from \code{QDM}. 
  }
  \item{silent}{
      logical value indicating if algorithm progress should be reported.
  }
}
\value{
     a list of with elements consisting of:
      \item{mhat.c}{matrix of bias corrected \code{m.c} values for the calibration period.}
      \item{mhat.p}{matrix of bias corrected \code{m.p} values for the projection period.}
}

\seealso{
\code{\link{QDM}, \link{MBCr}, \link{MRS}, \link{MBCn} \link{escore}}
}
\references{
Cannon, A.J., 2016. Multivariate bias correction of climate model outputs:
Matching marginal distributions and inter-variable dependence structure.
Journal of Climate, 29:7045-7064. doi:10.1175/JCLI-D-15-0679.1

Cannon, A.J., S.R. Sobie, and T.Q. Murdock, 2015. Bias correction
of simulated precipitation by quantile mapping: How well do methods
preserve relative changes in quantiles and extremes? Journal of
Climate, 28:6938-6959. doi:10.1175/JCLI-D-14-00754.1
}
