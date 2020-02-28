gapStatistic <- structure(function( #Gap Statistic for Estimating the Number of Clusters
###Given a genotype matrix, this function performs a constrained hierarchical clustering of the variables and infers the optimal number of clusters using the Gap statistic method.               
                                   X,
### genotype matrix coded in 1, 2 or 3. The variables to be clustered should be in columns.
                                   min.nc=1,
### integer; the minimum number of clusters to consider. Defaults to 1.
                                   max.nc=ncol(X)-1,
### integer; the maximum number of clusters to consider. Must be at least 2 and at most (p-1) with p the number of variables to be clustered.
                                   traceWB=NULL,
### An optional \code{(p-1) x B} matrix, giving within-cluster dispersion measures for B reference data sets.  Such measures can be obtained using  \code{\link{getTraceW}} for a given reference data set.
                                   B=500,
### integer; number of the bootstrap samples. Defaults to 500.
                                   flv=c("sample", "sampleX"),
### string; either "sample" or "sampleX". If flv="sample", each column of a reference data set is drawn independently from a uniform distribution over the values {1, 2, 3}. Else if flv="sampleX", each reference data set is resampled over the range of observed values (of X) for this column. Default "sample",
                                   smooth=FALSE,
### perform local smoothing of gap statistic (using \code{\link[stats]{smooth}})
                                   ...
### optionally further arguments for the function \code{\link{cWard}}.
                                   ){
  ## check arguments
  p <- ncol(X)
  N <- nrow(X)
  stopifnot(min.nc<max.nc)
  stopifnot(min.nc>=1)
  stopifnot(max.nc<p)
  ncs <- seq.int(from=min.nc, to=max.nc)
  flv <- match.arg(flv)
  
  best.k <- NULL
  result <- array(0, dim=c(length(ncs), 3), dimnames=list(nc=ncs,
                                                    gap=c("E.Wks", "SE.Wks", "W")))
  
  ## W.k
  res <- cWard(X, ...)
  W <- res$traceW[ncs, 2]  

  ## trace(W)
  if (!is.null(traceWB)) {
    stopifnot(nrow(traceWB)==(p-1))
    Wks <- traceWB[ncs, ]
  } else {
    Wks <- matrix(0, length(ncs), B)
    ## unstructured data
    if (flv=="sampleX") {
      for (bb in 1:B){
        X.sim <- apply(X, 2, function(v){
          sample(v)
        })
        Wks[, bb] <- getTraceW(N, p, X.sim, ...)[ncs]
      }
    } else if (flv=="sample") {
      for (bb in 1:B){
        Wks[, bb] <- getTraceW(N, p, ...)[ncs]
      }
    }
  }
  
  ## best K
  result[,1] <- rowMeans(Wks)
  result[,2] <- sqrt((1 + 1/B)*apply(Wks, 1, var))
  result[,3] <- W
  
  W <- result[,3]
  E.Wks <- result[,1]
  SE.Wks <- result[,2]
  gapK <-  E.Wks - W
  m <- diff(gapK, lag=1, differences=1) - SE.Wks[2:length(SE.Wks)]
  
  if (smooth) {
    m <- smooth(m)
  }
  best.k <- min(as.integer(names(which(m<0)))) - 1   ## names start at min.nc+1
  
  
  ##value<< A list with elements
  list(stats=result, ##<< a matrix with \code{max.nc-min.nc+1} rows and 3 columns named "E.Wks", "SE.Wks" and "W" corresponding respectively to the mean and the standard deviation of the reference data and the within-group dispersion measures at each step of clustering.
       best.k=best.k, ##<< estimated number of clusters. 
       tree=res, ##<<object of class \code{hclust} related to the clustering of X.
       stat=m  ## A numeric vector of length p-1, the statistic to be minimized
       )
  ##references<< R. Tibshirani, G. Walther, and T. Hastie. Estimating the number of clusters in a data set via the gap statistic. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 63(2):411-423, 2001.
  ##seealso<< \code{\link{hclust}}, \code{\link{cWard}}, \code{\link{grplassoCWard}}.
}, ex=function() {
  set.seed(5)

  N <- 100
  p <- 30

  ## 1- No group structure: k=1 should be optimal
  ## we use B=50 to keep it fast..
  X <- matrix(ceiling(runif(N*p)*3), N, p)
  gapStatistic(X, B=50)$best.k
  
  ## 2- Structured data: 6 clusters
  r2 <- 0.2
  corr <- 0.7
  blockSizes <- c(2:6, 5)
  p <- sum(blockSizes)
  sig.blocks <- 5
  nb.per.block <- 1
  betas <- simBeta(blockSizes, sig.blocks, nb.per.block)
  sim <- simulation(N, betas$betaMat[,"betaSNP"], betas$blockSizes, corr, r2, minMAF=0.45)
  Z <- sim$X
  ncs <- 1:10
  gapS <- gapStatistic(Z, min.nc=min(ncs), max.nc=max(ncs), B=50)
  print(gapS$best.k)

  ## plots like in Tibshirani et. al. 2001
  matplot(gapS$stats[, c(3, 1)], pch=c("O", "E"), lty=1, type="b",
          ylab="Observed vs Expected Wk", xlab="number of clusters k")

  gapK <- gapS$stats[,3]-gapS$stats[,1]
  plot(ncs, gapK, type="o", lty=1, pch=20,
       ylab="Gap", xlab="number of clusters k")

  ## externalization of the calculation of trace(W)
  B <- 50
  traceWB <- replicate(B, getTraceW(N, p))
  gapS <- gapStatistic(Z, min.nc=min(ncs), max.nc=max(ncs), traceWB=traceWB)
})

############################################################################
## HISTORY:
## 2014-07-19
## o Added argument 'smooth'.
## o Added return value 'stat' (mainly for testing).
## 2014-07-16
## o Externalization of the calculation of trace(W).
## 2013-...
## o Created.
############################################################################

