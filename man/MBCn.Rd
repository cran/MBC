\name{MBCn}
\alias{MBCn}
\title{
Multivariate bias correction (N-pdft)
}
\description{
Multivariate bias correction that matches the multivariate
distribution using \code{\link{QDM}} and the N-dimensional
probability density function transform (N-pdft) following
Cannon (2017).
}
\usage{
MBCn(o.c, m.c, m.p, iter=30, ratio.seq=rep(FALSE, ncol(o.c)),
     trace=0.05, trace.calc=0.5*trace, jitter.factor=0,
     n.tau=NULL, ratio.max=2, ratio.max.trace=10*trace,
     ties='first', qmap.precalc=FALSE, rot.seq=NULL,
     silent=FALSE, n.escore=0, return.all=FALSE, subsample=NULL,
     pp.type=7)
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
  \item{rot.seq}{
      use a supplied list of random rotation matrices. \code{NULL} generates on
      the fly.
  }
  \item{silent}{
      logical value indicating if algorithm progress should be reported.
  }
  \item{n.escore}{
      number of cases used to calculate the energy distance when monitoring 
      convergence.
  }
  \item{return.all}{
      logical value indicating whether results from all iterations are returned.
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
\code{\link{QDM}, \link{MBCp}, \link{MBCr}, \link{MRS}, \link{escore}, \link{rot.random}}
}
\references{
Cannon, A.J., 2017. Multivariate quantile mapping bias correction: An 
N-dimensional probability density function transform for climate model
simulations of multiple variables.
Climate Dynamics. doi:10.1007/s00382-017-3580-6

Cannon, A.J., S.R. Sobie, and T.Q. Murdock, 2015. Bias correction
of simulated precipitation by quantile mapping: How well do methods
preserve relative changes in quantiles and extremes? Journal of
Climate, 28:6938-6959. doi:10.1175/JCLI-D-14-00754.1

Pitié, F., A.C. Kokaram, and R. Dahyot, 2005. N-dimensional probability
density function transfer and its application to color transfer.
In Tenth IEEE International Conference on Computer Vision, 2005.
ICCV 2005. (Vol. 2, pp. 1434-1439). IEEE.

Pitié, F., A.C. Kokaram, and R. Dahyot, 2007. Automated colour grading
using colour distribution transfer. Computer Vision and Image
Understanding, 107(1):123-137.
}
