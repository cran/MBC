################################################################################
# MBC-QDM.R - Multivariate bias correction based on quantile delta mapping
# and iterative application of Cholesky decomposition rescaling (MBCp and MBCr)
# Multivariate bias correction based on quantile delta mapping and the 
# N-dimensional pdf transform (MBCn)
# Alex J. Cannon (alex.cannon@canada.ca)
################################################################################

library(Matrix)
library(energy)

QDM <-
# Quantile delta mapping bias correction for preserving changes in quantiles
# Note: QDM is equivalent to the equidistant and equiratio forms of quantile
# mapping (Cannon et al., 2015).
# Cannon, A.J., Sobie, S.R., and Murdock, T.Q. 2015. Bias correction of
#  simulated precipitation by quantile mapping: How well do methods preserve
#  relative changes in quantiles and extremes? Journal of Climate,
#  28: 6938-6959. doi:10.1175/JCLI-D-14-00754.1
function(o.c, m.c, m.p, ratio=FALSE, trace=0.05, trace.calc=0.5*trace,
         jitter.factor=0, n.tau=NULL, ratio.max=2, ratio.max.trace=10*trace,
         ECBC=FALSE, ties='first', subsample=NULL, pp.type=7){
    # o = vector of observed values; m = vector of modelled values
    # c = current period;  p = projected period
    # ratio = TRUE --> preserve relative trends in a ratio variable
    # trace = 0.05 --> replace values less than trace with exact zeros
    # trace.calc = 0.5*trace --> treat values below trace.calc as censored
    # jitter.factor = 0.01 --> jitter to accommodate ties
    # n.tau = NULL --> number of empirical quantiles (NULL=sample length)
    # ratio.max = 2 --> maximum delta when values are less than ratio.max.trace
    # ratio.max.trace = 10*trace --> values below which ratio.max is applied
    # ECBC = TRUE --> apply Schaake shuffle to enforce o.c temporal sequencing
    # subsample = NULL --> use this number of repeated subsamples of size n.tau
    #  to calculate empirical quantiles (e.g., when o.c, m.c, and m.p are of
    #  very different size)
    # pp.type = 7 --> plotting position type used in quantile
    # tau.m-p = F.m-p(x.m-p)
    # delta.m = x.m-p {/,-} F.m-c^-1(tau.m-p)
    # xhat.m-p = F.o-c^-1(tau.m-p) {*,+} delta.m
    #
    # If jitter.factor > 0, apply a small amount of jitter to accommodate ties
    # due to limited measurement precision
    if(jitter.factor==0 && 
      (length(unique(o.c))==1 ||
       length(unique(m.c))==1 ||
       length(unique(m.p))==1)){
        jitter.factor <- sqrt(.Machine$double.eps)
    }
    if(jitter.factor > 0){
        o.c <- jitter(o.c, jitter.factor)
        m.c <- jitter(m.c, jitter.factor)
        m.p <- jitter(m.p, jitter.factor)
    }
    # For ratio data, treat exact zeros as left censored values less than
    # trace.calc
    if(ratio){
        epsilon <- .Machine$double.eps
        o.c[o.c < trace.calc] <- runif(sum(o.c < trace.calc), min=epsilon,
                                       max=trace.calc)
        m.c[m.c < trace.calc] <- runif(sum(m.c < trace.calc), min=epsilon,
                                       max=trace.calc)
        m.p[m.p < trace.calc] <- runif(sum(m.p < trace.calc), min=epsilon,
                                       max=trace.calc)
    }
    # Calculate empirical quantiles
    n <- length(m.p)
    if(is.null(n.tau)) n.tau <- n
    tau <- seq(0, 1, length=n.tau)
    if(!is.null(subsample)){
        quant.o.c <- rowMeans(apply(replicate(subsample,
                              sample(o.c, size=length(tau))),
                              2, quantile, probs=tau, type=pp.type))
        quant.m.c <- rowMeans(apply(replicate(subsample,
                              sample(m.c, size=length(tau))),
                              2, quantile, probs=tau, type=pp.type))
        quant.m.p <- rowMeans(apply(replicate(subsample,
                              sample(m.p, size=length(tau))),
                              2, quantile, probs=tau, type=pp.type))
    } else{
        quant.o.c <- quantile(o.c, tau, type=pp.type)
        quant.m.c <- quantile(m.c, tau, type=pp.type)
        quant.m.p <- quantile(m.p, tau, type=pp.type)
    }
    # Apply quantile delta mapping bias correction
    tau.m.p <- approx(quant.m.p, tau, m.p, rule=2)$y    
    if(ratio){
        approx.t.qmc.tmp <- approx(tau, quant.m.c, tau.m.p, rule=2)$y
        delta.m <- m.p/approx.t.qmc.tmp
        delta.m[(delta.m > ratio.max) &
                (approx.t.qmc.tmp < ratio.max.trace)] <- ratio.max
        mhat.p <- approx(tau, quant.o.c, tau.m.p, rule=2)$y*delta.m
    } else{
        delta.m <- m.p - approx(tau, quant.m.c, tau.m.p, rule=2)$y
        mhat.p <- approx(tau, quant.o.c, tau.m.p, rule=2)$y + delta.m
    }
    mhat.c <- approx(quant.m.c, quant.o.c, m.c, rule=2)$y
    # For ratio data, set values less than trace to zero
    if(ratio){
        mhat.c[mhat.c < trace] <- 0
        mhat.p[mhat.p < trace] <- 0
    }
    if(ECBC){
        # empirical copula coupling/Schaake shuffle
        if(length(mhat.p)==length(o.c)){
            mhat.p <- sort(mhat.p)[rank(o.c, ties.method=ties)]
        } else{
            stop('Schaake shuffle failed due to incompatible lengths')
        }
    }
    list(mhat.c=mhat.c, mhat.p=mhat.p)
}

################################################################################
# Multivariate goodness-of-fit scoring function

escore <-
# Energy score for assessing the equality of two multivariate samples
# Székely, G.J. and Rizzo, M.L. 2013. Energy statistics: A class of statistics
#  based on distances. Journal of Statistical Planning and Inference, 143(8),
#  1249-1272. doi:10.1016/j.jspi.2013.03.018
# Baringhaus, L. and Franz, C. 2004. On a new multivariate two-sample test.
#  Journal of Multivariate Analysis, 88(1), 190-206.
#  doi:10.1016/S0047-259X(03)00079-4
function (x, y, scale.x = FALSE, n.cases = NULL, alpha = 1, method = "cluster")
{
    n.x <- nrow(x)
    n.y <- nrow(y)
    if (scale.x) {
        x <- scale(x)
        y <- scale(y, center = attr(x, "scaled:center"), scale = attr(x,
            "scaled:scale"))
    }
    if (!is.null(n.cases)) {
        n.cases <- min(n.x, n.y, n.cases)
        x <- x[sample(n.x, size = n.cases), , drop = FALSE]
        y <- y[sample(n.y, size = n.cases), , drop = FALSE]
        n.x <- n.cases
        n.y <- n.cases
    }
    edist(rbind(x, y), sizes = c(n.x, n.y), distance = FALSE,
        alpha = alpha, method = method)[1]/2
}

################################################################################
# Multivariate bias correction based on iterative application of quantile
# mapping (or ranking) and multivariate rescaling via Cholesky
# decomposition of the covariance matrix. Results in simulated marginal 
# distributions and Pearson or Spearman rank correlations that match
# observations
# Cannon, A.J. 2016. Multivariate Bias Correction of Climate Model Outputs:
#  Matching Marginal Distributions and Inter-variable Dependence Structure.
#  Journal of Climate, doi:10.1175/JCLI-D-15-0679.1

MRS <-
# Multivariate rescaling based on Cholesky decomposition of the covariance
#  matrix
# Scheuer, E.M., Stoller, D.S., 1962. On the generation of normal random
#  vectors. Technometrics, 4(2), 278-281.
# Bürger, G., Schulla, J., & Werner, A.T. 2011. Estimates of future flow,
#  including extremes, of the Columbia River headwaters. Wat Resour Res 47(10),
#  W10520, doi:10.1029/2010WR009716
function(o.c, m.c, m.p, o.c.chol=NULL, o.p.chol=NULL, m.c.chol=NULL,
         m.p.chol=NULL){
    # Center based on multivariate means
    o.c.mean <- colMeans(o.c)
    m.c.mean <- colMeans(m.c)
    m.p.mean <- colMeans(m.p)
    o.c <- sweep(o.c, 2, o.c.mean, '-')
    m.c <- sweep(m.c, 2, m.c.mean, '-')
    m.p <- sweep(m.p, 2, m.p.mean, '-')
    # Cholesky decomposition of covariance matrix
    # If !is.null(o.p.chol) --> projected target
    if(is.null(o.c.chol)) o.c.chol <- chol(cov(o.c))
    if(is.null(o.p.chol)) o.p.chol <- chol(cov(o.c))
    if(is.null(m.c.chol)) m.c.chol <- chol(cov(m.c))
    if(is.null(m.p.chol)) m.p.chol <- chol(cov(m.c))
    # Bias correction factors
    mbcfactor <- solve(m.c.chol) %*% o.c.chol
    mbpfactor <- solve(m.p.chol) %*% o.p.chol
    # Multivariate bias correction
    mbc.c <- m.c %*% mbcfactor
    mbc.p <- m.p %*% mbpfactor
    # Recenter and account for change in means
    mbc.c <- sweep(mbc.c, 2, o.c.mean, '+')
    mbc.p <- sweep(mbc.p, 2, o.c.mean, '+')
    mbc.p <- sweep(mbc.p, 2, m.p.mean-m.c.mean, '+')
    list(mhat.c=mbc.c, mhat.p=mbc.p)
}

MBCr <-
# Multivariate quantile mapping bias correction (Spearman correlation)
# Cannon, A.J., 2016. Multivariate bias correction of climate model outputs: 
#  matching marginal distributions and inter-variable dependence structure.
#  Journal of Climate, 29(19):7045–7064. doi:10.1175/JCLI-D-15-0679.1
function(o.c, m.c, m.p, iter=20, cor.thresh=1e-4,
         ratio.seq=rep(FALSE, ncol(o.c)), trace=0.05,
         trace.calc=0.5*trace, jitter.factor=0, n.tau=NULL, ratio.max=2,
         ratio.max.trace=10*trace, ties='first', qmap.precalc=FALSE,
         silent=FALSE, subsample=NULL, pp.type=7){
    if(length(trace.calc)==1)
        trace.calc <- rep(trace.calc, ncol(o.c))
    if(length(trace)==1)
        trace <- rep(trace, ncol(o.c))
    if(length(jitter.factor)==1)
        jitter.factor <- rep(jitter.factor, ncol(o.c))
    if(length(ratio.max))
        ratio.max <- rep(ratio.max, ncol(o.c))
    if(length(ratio.max.trace)==1)
        ratio.max.trace <- rep(ratio.max.trace, ncol(o.c))
    m.c.qmap <- m.c
    m.p.qmap <- m.p
    if(!qmap.precalc){
        # Quantile delta mapping bias correction
        for(i in seq(ncol(o.c))){
            fit.qmap <- QDM(o.c=o.c[,i], m.c=m.c[,i], m.p=m.p[,i],
                            ratio=ratio.seq[i], trace.calc=trace.calc[i],
                            trace=trace[i], jitter.factor=jitter.factor[i],
                            n.tau=n.tau, ratio.max=ratio.max[i],
                            ratio.max.trace=ratio.max.trace[i],
                            subsample=subsample, pp.type=pp.type)
            m.c.qmap[,i] <- fit.qmap$mhat.c
            m.p.qmap[,i] <- fit.qmap$mhat.p
        }
    }
    # Ordinal ranks of observed and modelled data
    o.c.r <- apply(o.c, 2, rank, ties.method=ties)
    m.c.r <- apply(m.c, 2, rank, ties.method=ties)
    m.p.r <- apply(m.p, 2, rank, ties.method=ties)
    m.c.i <- m.c.r
    if(cor.thresh > 0){
        # Spearman correlation to assess convergence
        cor.i <- cor(m.c.r)
        cor.i[is.na(cor.i)] <- 0
    } else{
        cor.diff <- 0
    }
    # Iterative MBC/reranking
    o.c.chol <- o.p.chol <- as.matrix(chol(nearPD(cov(o.c.r))$mat))
    for(i in seq(iter)){
        m.c.chol <- m.p.chol <- as.matrix(chol(nearPD(cov(m.c.r))$mat))
        fit.mbc <- MRS(o.c=o.c.r, m.c=m.c.r, m.p=m.p.r, o.c.chol=o.c.chol,
                       o.p.chol=o.p.chol, m.c.chol=m.c.chol, m.p.chol=m.p.chol)
        m.c.r <- apply(fit.mbc$mhat.c, 2, rank, ties.method=ties)
        m.p.r <- apply(fit.mbc$mhat.p, 2, rank, ties.method=ties)
        if(cor.thresh > 0){
            # Check on Spearman correlation convergence
            cor.j <- cor(m.c.r)
            cor.j[is.na(cor.j)] <- 0
            cor.diff <- mean(abs(cor.j-cor.i))
            cor.i <- cor.j
        }
        if(!silent){
            cat(i, mean(m.c.r==m.c.i), cor.diff, '')
        }
        if(cor.diff < cor.thresh) break
        if(identical(m.c.r, m.c.i)) break
        m.c.i <- m.c.r
    }
    if(!silent) cat('\n')
    for(i in seq(ncol(o.c))){
        # Replace ordinal ranks with QDM outputs
        m.c.r[,i] <- sort(m.c.qmap[,i])[m.c.r[,i]]
        m.p.r[,i] <- sort(m.p.qmap[,i])[m.p.r[,i]]
    }
    list(mhat.c=m.c.r, mhat.p=m.p.r)
}

MBCp <-
# Multivariate quantile mapping bias correction (Pearson correlation)
# Cannon, A.J., 2016. Multivariate bias correction of climate model outputs: 
#  matching marginal distributions and inter-variable dependence structure.
#  Journal of Climate, 29(19):7045–7064. doi:10.1175/JCLI-D-15-0679.1
function(o.c, m.c, m.p, iter=20, cor.thresh=1e-4,
         ratio.seq=rep(FALSE, ncol(o.c)), trace=0.05, trace.calc=0.5*trace,
         jitter.factor=0, n.tau=NULL, ratio.max=2, ratio.max.trace=10*trace,
         ties='first', qmap.precalc=FALSE, silent=FALSE, subsample=NULL,
         pp.type=7){
    if(length(trace.calc)==1)
        trace.calc <- rep(trace.calc, ncol(o.c))
    if(length(trace)==1)
        trace <- rep(trace, ncol(o.c))
    if(length(jitter.factor)==1)
        jitter.factor <- rep(jitter.factor, ncol(o.c))
    if(length(ratio.max))
        ratio.max <- rep(ratio.max, ncol(o.c))
    if(length(ratio.max.trace)==1)
        ratio.max.trace <- rep(ratio.max.trace, ncol(o.c))
    m.c.qmap <- m.c
    m.p.qmap <- m.p
    if(!qmap.precalc){
        # Quantile delta mapping bias correction
        for(i in seq(ncol(o.c))){
            fit.qmap <- QDM(o.c=o.c[,i], m.c=m.c[,i], m.p=m.p[,i],
                            ratio=ratio.seq[i], trace.calc=trace.calc[i],
                            trace=trace[i], jitter.factor=jitter.factor[i],
                            n.tau=n.tau, ratio.max=ratio.max[i],
                            ratio.max.trace=ratio.max.trace[i],
                            subsample=subsample, pp.type=pp.type)
            m.c.qmap[,i] <- fit.qmap$mhat.c
            m.p.qmap[,i] <- fit.qmap$mhat.p
        }
    }
    m.c <- m.c.qmap
    m.p <- m.p.qmap
    # Pearson correlation to assess convergence
    if(cor.thresh > 0){
        cor.i <- cor(m.c)
        cor.i[is.na(cor.i)] <- 0
    }
    o.c.chol <- o.p.chol <- as.matrix(chol(nearPD(cov(o.c))$mat))
    # Iterative MBC/QDM
    for(i in seq(iter)){
        m.c.chol <- m.p.chol <- as.matrix(chol(nearPD(cov(m.c))$mat))
        fit.mbc <- MRS(o.c=o.c, m.c=m.c, m.p=m.p, o.c.chol=o.c.chol,
                       o.p.chol=o.p.chol, m.c.chol=m.c.chol, m.p.chol=m.p.chol)
        m.c <- fit.mbc$mhat.c
        m.p <- fit.mbc$mhat.p
        for(j in seq(ncol(o.c))){
            fit.qmap <- QDM(o.c=o.c[,j], m.c=m.c[,j], m.p=m.p[,j], ratio=FALSE,
                            n.tau=n.tau, pp.type=pp.type)
            m.c[,j] <- fit.qmap$mhat.c
            m.p[,j] <- fit.qmap$mhat.p
        }
        # Check on Pearson correlation convergence
        if(cor.thresh > 0){
            cor.j <- cor(m.c)
            cor.j[is.na(cor.j)] <- 0
            cor.diff <- mean(abs(cor.j-cor.i))
            cor.i <- cor.j
        } else{
            cor.diff <- 0
        }
        if(!silent) cat(i, cor.diff, '')
        if(cor.diff < cor.thresh) break
    }
    if(!silent) cat('\n')
    # Replace with shuffled QDM elements
    for(i in seq(ncol(o.c))){
        m.c[,i] <- sort(m.c.qmap[,i])[rank(m.c[,i], ties.method=ties)]
        m.p[,i] <- sort(m.p.qmap[,i])[rank(m.p[,i], ties.method=ties)]
    }
    list(mhat.c=m.c, mhat.p=m.p)
}

################################################################################
# Multivariate bias correction based on iterative application of random
# orthogonal rotation and quantile mapping (N-dimensional pdf transfer)
# Pitié, F., Kokaram, A.C., and Dahyot, R. 2005. N-dimensional probability
#  density function transfer and its application to color transfer.
#  In Tenth IEEE International Conference on Computer Vision, 2005. ICCV 2005.
#  (Vol. 2, pp. 1434-1439). IEEE.
# Pitié, F., Kokaram, A.C., and Dahyot, R. 2007. Automated colour grading
#  using colour distribution transfer. Computer Vision and Image Understanding,
#  107(1), 123-137.

rot.random <-
# Random orthogonal rotation
function(k) {
  rand <- matrix(rnorm(k * k), ncol=k)
  QRd <- qr(rand)
  Q <- qr.Q(QRd)
  R <- qr.R(QRd)
  diagR <- diag(R)
  rot <- Q %*% diag(diagR/abs(diagR))
  return(rot)
}

MBCn <- 
# Multivariate quantile mapping bias correction (N-dimensional pdf transfer)
# Cannon, A.J., 2017. Multivariate quantile mapping bias correction: An 
#  N-dimensional probability density function transform for climate model
#  simulations of multiple variables. Climate Dynamics.
#  doi:10.1007/s00382-017-3580-6
function(o.c, m.c, m.p, iter=30, ratio.seq=rep(FALSE, ncol(o.c)),
         trace=0.05, trace.calc=0.5*trace, jitter.factor=0, n.tau=NULL,
         ratio.max=2, ratio.max.trace=10*trace, ties='first',
         qmap.precalc=FALSE, rot.seq=NULL, silent=FALSE, n.escore=0,
         return.all=FALSE, subsample=NULL, pp.type=7){
    if(!is.null(rot.seq)){
        if(length(rot.seq)!=iter){
            stop('length(rot.seq) != iter')
        }
    }
    if(length(trace.calc)==1)
        trace.calc <- rep(trace.calc, ncol(o.c))
    if(length(trace)==1)
        trace <- rep(trace, ncol(o.c))
    if(length(jitter.factor)==1)
        jitter.factor <- rep(jitter.factor, ncol(o.c))
    if(length(ratio.max))
        ratio.max <- rep(ratio.max, ncol(o.c))
    if(length(ratio.max.trace)==1)
        ratio.max.trace <- rep(ratio.max.trace, ncol(o.c))
    # Energy score (rescaled)
    escore.iter <- rep(NA, iter+2)
    if(n.escore > 0){
        n.escore <- min(nrow(o.c), nrow(m.c), n.escore)
        escore.cases.o.c <- unique(suppressWarnings(matrix(seq(nrow(o.c)),
                                   ncol=n.escore)[1,]))
        escore.cases.m.c <- unique(suppressWarnings(matrix(seq(nrow(m.c)),
                                   ncol=n.escore)[1,]))
        escore.iter[1] <- escore(x=o.c[escore.cases.o.c,,drop=FALSE],
                                 y=m.c[escore.cases.m.c,,drop=FALSE],
                                 scale.x=TRUE)
        if(!silent) cat('RAW', escore.iter[1], ': ')        
    }
    m.c.qmap <- m.c
    m.p.qmap <- m.p
    if(!qmap.precalc){
        # Quantile delta mapping bias correction
        for(i in seq(ncol(o.c))){
            fit.qmap <- QDM(o.c=o.c[,i], m.c=m.c[,i], m.p=m.p[,i],
                            ratio=ratio.seq[i], trace.calc=trace.calc[i],
                            trace=trace[i], jitter.factor=jitter.factor[i],
                            n.tau=n.tau, ratio.max=ratio.max[i],
                            ratio.max.trace=ratio.max.trace[i],
                            subsample=subsample, pp.type=pp.type)
            m.c.qmap[,i] <- fit.qmap$mhat.c
            m.p.qmap[,i] <- fit.qmap$mhat.p
        }
    }
    m.c <- m.c.qmap
    m.p <- m.p.qmap
    # Energy score (QDM)
    if(n.escore > 0){
        escore.iter[2] <- escore(x=o.c[escore.cases.o.c,,drop=FALSE],
                                 y=m.c[escore.cases.m.c,,drop=FALSE],
                                 scale.x=TRUE)
        if(!silent) cat('QDM', escore.iter[2], ': ')        
    }
    # Standardize observations
    m.iter <- vector('list', iter)
    o.c.mean <- colMeans(o.c)
    o.c.sdev <- apply(o.c, 2, sd)
    o.c.sdev[o.c.sdev < .Machine$double.eps] <- 1    
    o.c <- scale(o.c, center=o.c.mean, scale=o.c.sdev)
    # Standardize model
    m.c.p <- rbind(m.c, m.p)
    m.c.p.mean <- colMeans(m.c.p)
    m.c.p.sdev <- apply(m.c.p, 2, sd)
    m.c.p.sdev[m.c.p.sdev < .Machine$double.eps] <- 1
    m.c.p <- scale(m.c.p, center=m.c.p.mean, scale=m.c.p.sdev)    
    Xt <- rbind(o.c, m.c.p)
    for(i in seq(iter)){
        if(!silent) cat(i, '')
        # Random orthogonal rotation
        if(is.null(rot.seq)){
            rot <- rot.random(ncol(o.c))
        } else{
            rot <- rot.seq[[i]]
        }
        Z <- Xt %*% rot
        Z.o.c <- Z[1:nrow(o.c),,drop=FALSE]
        Z.m.c <- Z[(nrow(o.c)+1):(nrow(o.c)+nrow(m.c)),,drop=FALSE]
        Z.m.p <- Z[(nrow(o.c)+nrow(m.c)+1):nrow(Z),,drop=FALSE]
        # Bias correct rotated variables using QDM
        for(j in seq(ncol(Z))){
            Z.qdm <- QDM(o.c=Z.o.c[,j], m.c=Z.m.c[,j], m.p=Z.m.p[,j],
                         ratio=FALSE, jitter.factor=jitter.factor[j],
                         n.tau=n.tau, pp.type=pp.type)
            Z.m.c[,j] <- Z.qdm$mhat.c
            Z.m.p[,j] <- Z.qdm$mhat.p
        }
        # Rotate back
        m.c <- Z.m.c %*% t(rot)
        m.p <- Z.m.p %*% t(rot)
        Xt <- rbind(o.c, m.c, m.p)
        # Energy score (MBCn)
        if(n.escore > 0){
            escore.iter[i+2] <- escore(x=o.c[escore.cases.o.c,,drop=FALSE],
                                       y=m.c[escore.cases.m.c,,drop=FALSE],
                                       scale.x=TRUE)
            if(!silent) cat(escore.iter[i+2], ': ')
        }
        if(return.all){
            m.c.i <- sweep(sweep(m.c, 2, attr(m.c.p, 'scaled:scale'), '*'), 2,
                           attr(m.c.p, 'scaled:center'), '+')
            m.p.i <- sweep(sweep(m.p, 2, attr(m.c.p, 'scaled:scale'), '*'), 2,
                           attr(m.c.p, 'scaled:center'), '+')
            m.iter[[i]] <- list(m.c=m.c.i, m.p=m.p.i)
        }
    }
    if(!silent) cat('\n')
    # Rescale back to original units
    m.c <- sweep(sweep(m.c, 2, m.c.p.sdev, '*'), 2, m.c.p.mean, '+')
    m.p <- sweep(sweep(m.p, 2, m.c.p.sdev, '*'), 2, m.c.p.mean, '+')
    # Replace npdft ordinal ranks with QDM outputs
    for(i in seq(ncol(o.c))){
        m.c[,i] <- sort(m.c.qmap[,i])[rank(m.c[,i], ties.method=ties)]
        m.p[,i] <- sort(m.p.qmap[,i])[rank(m.p[,i], ties.method=ties)]
    }
    names(escore.iter)[1:2] <- c('RAW', 'QM')
    names(escore.iter)[-c(1:2)] <- seq(iter)
    list(mhat.c=m.c, mhat.p=m.p, escore.iter=escore.iter, m.iter=m.iter)
}

################################################################################
