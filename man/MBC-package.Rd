\name{MBC-package}
\alias{MBC-package}
\alias{MBC}
\docType{package}
\title{
Multivariate Bias Correction of Climate Model Outputs
}
\description{
Calibrate and apply multivariate bias correction
algorithms for climate model simulations of multiple
climate variables. Three iterative methods are supported: (i) MBC
Pearson correlation (MBCp), (ii) MBC rank correlation (MBCr),
and (iii) MBC N-dimensional probability density function
transform (MBCn). The first two, MBCp and MBCr
(Cannon, 2016), match marginal distributions and inter-variable
dependence structure. Dependence structure can be measured
either by the Pearson correlation (\code{\link{MBCp}})
or by the Spearman rank correlation (\code{\link{MBCr}}).
The energy distance score (\code{\link{escore}}) is recommended
for model selection. The third, \code{\link{MBCn}} (Cannon, 2018),
which operates on the full multivariate distribution, is more
flexible and can be considered to be a multivariate analogue
of univariate quantile mapping. All aspects of the observed
distribution are transferred to the climate model simulations.
In each of the three methods, marginal distributions are corrected
by the change-preserving quantile delta mapping (\code{\link{QDM}})
algorithm (Cannon et al., 2015). Finally, an implementation
of the Rank Resampling for Distributions and Dependences (R2D2)
method introduced by Vrac (2018) is also included.

An example application of the three MBC methods using the
\code{cccma} dataset can be run via:

\code{example(MBC, run.dontrun=TRUE)}

Note: because empirical quantiles and their changes are used by
\code{\link{QDM}}, sample sizes of the observed, model calibration,
and model projection datasets should be approximately equal.

}
\details{
\tabular{ll}{
Package: \tab MBC\cr
Type: \tab Package\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}

\references{
Cannon, A.J., 2018. Multivariate quantile mapping bias correction: An 
N-dimensional probability density function transform for climate model
simulations of multiple variables. Climate Dynamics, 50(1-2):31-49.
doi:10.1007/s00382-017-3580-6

Cannon, A.J., 2016. Multivariate bias correction of climate model output:
Matching marginal distributions and inter-variable dependence structure.
Journal of Climate, 29:7045-7064. doi:10.1175/JCLI-D-15-0679.1

Cannon, A.J., S.R. Sobie, and T.Q. Murdock, 2015. Bias correction
of simulated precipitation by quantile mapping: How well do methods
preserve relative changes in quantiles and extremes? Journal of
Climate, 28:6938-6959. doi:10.1175/JCLI-D-14-00754.1

Francois, B., M. Vrac, A.J. Cannon, Y. Robin, and D. Allard, 2020.
Multivariate bias corrections of climate simulations: Which benefits
for which losses? Earth System Dynamics, 11:537-562.
doi:10.5194/esd-11-537-2020

Vrac, M., 2018. Multivariate bias adjustment of high-dimensional
climate simulations: the Rank Resampling for Distributions and
Dependences (R2D2) bias correction. Hydrology and Earth System
Sciences, 22:3175-3196. doi:10.5194/hess-22-3175-2018

}
\keyword{ package }

\seealso{
\code{\link{QDM}, \link{MBCp}, \link{MBCr}, \link{MBCn}, \link{R2D2}, \link{escore}, \link{rot.random}, \link{cccma}}
}


\examples{
\dontrun{
data(cccma)
set.seed(1)

# Univariate quantile mapping
qdm.c <- cccma$gcm.c*0
qdm.p <- cccma$gcm.p*0
for(i in seq(ncol(cccma$gcm.c))){
    fit.qdm <- QDM(o.c=cccma$rcm.c[,i], m.c=cccma$gcm.c[,i],
                   m.p=cccma$gcm.p[,i], ratio=cccma$ratio.seq[i],
                   trace=cccma$trace[i])
    qdm.c[,i] <- fit.qdm$mhat.c
    qdm.p[,i] <- fit.qdm$mhat.p
}

# Multivariate MBCp bias correction
fit.mbcp <- MBCp(o.c=cccma$rcm.c, m.c=cccma$gcm.c,
                 m.p=cccma$gcm.p, ratio.seq=cccma$ratio.seq,
                 trace=cccma$trace)
mbcp.c <- fit.mbcp$mhat.c
mbcp.p <- fit.mbcp$mhat.p

# Multivariate MBCr bias correction
fit.mbcr <- MBCr(o.c=cccma$rcm.c, m.c=cccma$gcm.c,
                 m.p=cccma$gcm.p, ratio.seq=cccma$ratio.seq,
                 trace=cccma$trace)
mbcr.c <- fit.mbcr$mhat.c
mbcr.p <- fit.mbcr$mhat.p

# Multivariate MBCn bias correction
fit.mbcn <- MBCn(o.c=cccma$rcm.c, m.c=cccma$gcm.c,
                 m.p=cccma$gcm.p, ratio.seq=cccma$ratio.seq,
                 trace=cccma$trace)
mbcn.c <- fit.mbcn$mhat.c
mbcn.p <- fit.mbcn$mhat.p
colnames(mbcn.c) <- colnames(mbcn.p) <-
    colnames(cccma$rcm.c)

# Correlation matrices (Pearson and Spearman)
# MBCp
dev.new()
par(mfrow=c(2, 2))
plot(c(cor(cccma$rcm.c)), c(cor(qdm.c)), col='black',
     pch=19, xlim=c(-1, 1), ylim=c(-1, 1),
     xlab='CanRCM4', ylab='CanESM2 MBCp',
     main='Pearson correlation\nMBCp calibration')
abline(0, 1)
grid()
points(c(cor(cccma$rcm.c)), c(cor(mbcp.c)), col='red')
plot(c(cor(cccma$rcm.p)), c(cor(qdm.p)),
     col='black', pch=19, xlim=c(-1, 1), ylim=c(-1, 1),
     xlab='CanRCM4', ylab='CanESM2 MBCp',
     main='Pearson correlation\nMBCp evaluation')
abline(0, 1)
grid()
points(c(cor(cccma$rcm.p)), c(cor(mbcp.p)), col='red')
plot(c(cor(cccma$rcm.c, m='s')), c(cor(qdm.c, m='s')),
     col='black', pch=19, xlim=c(-1, 1), ylim=c(-1, 1),
     xlab='CanRCM4', ylab='CanESM2 MBCp',
     main='Spearman correlation\nMBCp calibration')
abline(0, 1)
grid()
points(c(cor(cccma$rcm.c, m='s')), c(cor(mbcp.c, m='s')),
       col='red')
plot(c(cor(cccma$rcm.p, m='s')), c(cor(qdm.p, m='s')),
     col='black', pch=19, xlim=c(-1, 1), ylim=c(-1, 1),
     xlab='CanRCM4', ylab='CanESM2 MBCp',
     main='Spearman correlation\nMBCp evaluation')
abline(0, 1)
grid()
points(c(cor(cccma$rcm.p, m='s')), c(cor(mbcp.p, m='s')),
       col='red')

# MBCr
dev.new()
par(mfrow=c(2, 2))
plot(c(cor(cccma$rcm.c)), c(cor(qdm.c)), col='black',
     pch=19, xlim=c(-1, 1), ylim=c(-1, 1),
     xlab='CanRCM4', ylab='CanESM2 MBCr',
     main='Pearson correlation\nMBCr calibration')
abline(0, 1)
grid()
points(c(cor(cccma$rcm.c)), c(cor(mbcr.c)), col='blue')
plot(c(cor(cccma$rcm.p)), c(cor(qdm.p)),
     col='black', pch=19, xlim=c(-1, 1), ylim=c(-1, 1),
     xlab='CanRCM4', ylab='CanESM2 MBCr',
     main='Pearson correlation\nMBCr evaluation')
abline(0, 1)
grid()
points(c(cor(cccma$rcm.p)), c(cor(mbcr.p)), col='blue')
plot(c(cor(cccma$rcm.c, m='s')), c(cor(qdm.c, m='s')),
     col='black', pch=19, xlim=c(-1, 1), ylim=c(-1, 1),
     xlab='CanRCM4', ylab='CanESM2 MBCr',
     main='Spearman correlation\nMBCr calibration')
abline(0, 1)
grid()
points(c(cor(cccma$rcm.c, m='s')), c(cor(mbcr.c, m='s')),
       col='blue')
plot(c(cor(cccma$rcm.p, m='s')), c(cor(qdm.p, m='s')),
     col='black', pch=19, xlim=c(-1, 1), ylim=c(-1, 1),
     xlab='CanRCM4', ylab='CanESM2 MBCr',
     main='Spearman correlation\nMBCr evaluation')
abline(0, 1)
grid()
points(c(cor(cccma$rcm.p, m='s')), c(cor(mbcr.p, m='s')),
       col='blue')

# MBCn
dev.new()
par(mfrow=c(2, 2))
plot(c(cor(cccma$rcm.c)), c(cor(qdm.c)), col='black',
     pch=19, xlim=c(-1, 1), ylim=c(-1, 1),
     xlab='CanRCM4', ylab='CanESM2 MBCn',
     main='Pearson correlation\nMBCn calibration')
abline(0, 1)
grid()
points(c(cor(cccma$rcm.c)), c(cor(mbcn.c)), col='orange')
plot(c(cor(cccma$rcm.p)), c(cor(qdm.p)),
     col='black', pch=19, xlim=c(-1, 1), ylim=c(-1, 1),
     xlab='CanRCM4', ylab='CanESM2 MBCn',
     main='Pearson correlation\nMBCn evaluation')
abline(0, 1)
grid()
points(c(cor(cccma$rcm.p)), c(cor(mbcn.p)), col='orange')
plot(c(cor(cccma$rcm.c, m='s')), c(cor(qdm.c, m='s')),
     col='black', pch=19, xlim=c(-1, 1), ylim=c(-1, 1),
     xlab='CanRCM4', ylab='CanESM2 MBCn',
     main='Spearman correlation\nMBCn calibration')
abline(0, 1)
grid()
points(c(cor(cccma$rcm.c, m='s')), c(cor(mbcn.c, m='s')),
       col='orange')
plot(c(cor(cccma$rcm.p, m='s')), c(cor(qdm.p, m='s')),
     col='black', pch=19, xlim=c(-1, 1), ylim=c(-1, 1),
     xlab='CanRCM4', ylab='CanESM2 MBCn',
     main='Spearman correlation\nMBCn evaluation')
abline(0, 1)
grid()
points(c(cor(cccma$rcm.p, m='s')), c(cor(mbcn.p, m='s')),
       col='orange')

# Pairwise scatterplots
dev.new()
pairs(cccma$gcm.c, main='CanESM2 calibration', col='#0000001A')
dev.new()
pairs(cccma$rcm.c, main='CanRCM4 calibration', col='#0000001A')
dev.new()
pairs(qdm.c, main='QDM calibration', col='#0000001A')
dev.new()
pairs(mbcp.c, main='MBCp calibration', col='#FF00001A')
dev.new()
pairs(mbcr.c, main='MBCr calibration', col='#0000FF1A')
dev.new()
pairs(mbcn.c, main='MBCn calibration', col='#FFA5001A')

# Energy distance skill score relative to univariate QDM
escore.qdm <- escore(cccma$rcm.p, qdm.p, scale.x=TRUE)
escore.mbcp <- escore(cccma$rcm.p, mbcp.p, scale.x=TRUE)
escore.mbcr <- escore(cccma$rcm.p, mbcr.p, scale.x=TRUE)
escore.mbcn <- escore(cccma$rcm.p, mbcn.p, scale.x=TRUE)

cat('ESS (MBCp):', 1-escore.mbcp/escore.qdm, '\n')
cat('ESS (MBCr):', 1-escore.mbcr/escore.qdm, '\n')
cat('ESS (MBCn):', 1-escore.mbcn/escore.qdm, '\n')
}
}
