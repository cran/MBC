\name{escore}
\alias{escore}
\title{
Energy distance score
}
\description{
Calculate the energy distance score measuring the statistical
discrepancy between samples \code{x} and \code{y} from two multivariate
distributions.
}
\usage{
escore(x, y, scale.x=FALSE, n.cases=NULL, alpha=1, method='cluster')
}
\arguments{
  \item{x}{
      numeric matrix.
  }
  \item{y}{
      numeric matrix.
  }
  \item{scale.x}{
      logical indicating whether data should be standardized based on \code{x}.
  }
  \item{n.cases}{
      the number of sub-sampled cases; \code{NULL} uses all data.
  }
  \item{alpha}{
      distance exponent in (0,2]
  }
  \item{method}{
      method used to weight the statistics
  }
}

\references{
Székely, G.J. and M.L. Rizzo, 2004. Testing for equal
distributions in high dimension, InterStat, November (5).

Székely, G.J. and M.L. Rizzo, 2013. Energy statistics: statistics
based on distances. Journal of Statistical Planning and Inference,
143(8):1249-1272. doi:10.1016/j.jspi.2013.03.018

Rizzo, M.L. and G.L. Székely, 2016. Energy distance. Wiley
Interdisciplinary Reviews: Computational Statistics, 8(1):27-38.
}
