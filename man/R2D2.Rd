\name{R2D2}
\alias{R2D2}
\title{
Multivariate bias correction (R2D2)
}
\description{
Multivariate bias correction that matches the multivariate
distribution using \code{\link{QDM}} and the R2D2 
algorithm following Vrac (2018).
}
\usage{
R2D2(o.c, m.c, m.p, ref.column=1, ratio.seq=rep(FALSE, ncol(o.c)),
     trace=0.05, trace.calc=0.5*trace, jitter.factor=0,
     n.tau=NULL, ratio.max=2, ratio.max.trace=10*trace,
     ties='first', qmap.precalc=FALSE, subsample=NULL, pp.type=7)
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
  \item{ref.column}{
      index of the reference column used for the 1D nearest neighbour matching
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
  \item{subsample}{
      use \code{subsample} draws of size \code{n.tau} to calculate initial
      empirical quantiles; if \code{NULL}, calculate normally.
  }
  \item{pp.type}{
      type of plotting position used in \code{quantile}.
  }
}
\value{
     a list of with elements consisting of:
      \item{mhat.c}{matrix of bias corrected \code{m.c} values for the calibration period.}
      \item{mhat.p}{matrix of bias corrected \code{m.p} values for the projection period.}
}

\seealso{
\code{\link{QDM}, \link{MBCp}, \link{MBCr}, \link{MRS}, \link{MBCn}}
}
\references{
Cannon, A.J., S.R. Sobie, and T.Q. Murdock, 2015. Bias correction
of simulated precipitation by quantile mapping: How well do methods
preserve relative changes in quantiles and extremes? Journal of
Climate, 28:6938-6959. doi:10.1175/JCLI-D-14-00754.1

Vrac, M., 2018. Multivariate bias adjustment of high-dimensional
climate simulations: the Rank Resampling for Distributions and
Dependences (R2D2) bias correction. Hydrology and Earth System
Sciences, 22:3175-3196. doi:10.5194/hess-22-3175-2018
}
