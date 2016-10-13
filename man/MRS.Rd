\name{MRS}
\alias{MRS}
\title{
Multivariate linear rescaling using Cholesky decomposition
}
\description{
Multivariate linear bias correction based on Cholesky decomposition of the
covariance matrix following Scheuer and Stoller (1962) and Bürger
et al. (2011). Bias correction matches the multivariate mean and covariance
structure.
}
\usage{
MRS(o.c, m.c, m.p, o.c.chol=NULL, o.p.chol=NULL, m.c.chol=NULL,
    m.p.chol=NULL)
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
  \item{o.c.chol}{
      precalculated Cholesky decomposition of the \code{o.c} covariance matrix; \code{NULL} calculates internally.
  }
  \item{o.p.chol}{
      precalculated Cholesky decomposition of the target \code{o.p} covariance matrix; \code{NULL} defaults to \code{o.c.chol}.
  }
  \item{m.c.chol}{
      precalculated Cholesky decomposition of the \code{m.c} covariance matrix; \code{NULL} calculates internally.
  }
  \item{m.p.chol}{
      precalculated Cholesky decomposition of the \code{m.p} covariance matrix; \code{NULL} calculates internally.
  }
}
\value{
     a list of with elements consisting of:
      \item{mhat.c}{matrix of bias corrected \code{m.c} values for the calibration period.}
      \item{mhat.p}{matrix of bias corrected \code{m.p} values for the projection period.}
}

\seealso{
\code{\link{MBCp}, \link{MBCr}}
}
\references{
Scheuer, E.M. and D.S. Stoller, 1962. On the generation of normal random
vectors. Technometrics, 4(2):278-281.

Bürger, G., J. Schulla, and A.T. Werner, 2011. Estimates of future flow,
including extremes, of the Columbia River headwaters. Water Resources
Research, 47(10):W10520. doi:10.1029/2010WR009716
}
