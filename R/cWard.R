####################
## WARD constrained
####################
## statistic: "R.squared", "D.prime", "D" ..
## returns a vector of p-1 similarities
.similarite <- function(X, statistic){
  l.d <- ld(X, depth=1, stats=statistic, symmetric = FALSE)
  sim <- slot(l.d, "x")
  return(sim)
}

simR2 <- structure(function( #Calculation of the \eqn{r^2} linkage disequilibrium measure
### This function returns the pairs of \eqn{r^2} LD measures between each column of \eqn{X[,i]} and each column of \eqn{X[,j]}. If i and j are not provided, this function returns a vector of the \eqn{p-1} \eqn{r^2} LD measures between the adjacent columns of X.
                            X,
### genotype matrix coded in {1, 2, 3} with the SNPs in columns.
                            i,
### (optional) vector of integers. Indices of the SNPs of the first group among the columns of X of length \eqn{l_1}. \eqn{l_1} must be between 1 and p.
                            j
### (optional) vector of integers. Indices of the SNPs of the second group among the columns of X of length \eqn{l_2}. \eqn{l_2} must be between 1 and p.
                            ){
  if (class(X)!="SnpMatrix"){
    X <- .snpM(X)
  }
  stats <- "R.squared"
  if (missing(i) && missing(j)) {
    res <- ld(X, depth=1, stats=stats)
    res <- slot(res, "x")
  } else if (length(i)>1 && length(j)>1) {
    res <- ld(X[,i], X[,j], depth=ncol(X)-1, stats=stats)
  } else {
    res <- ld(X[,i], X[,j], stats=stats)
  }
  ##value<< If \code{i} and \code{j} are provided, a \code{l_1 x l_2} matrix of the \eqn{r^2} LD measures is returned. Else, a (p-1) vector of the \eqn{r^2} LD measures between the adjacent columns of X is returned.
  res
  ##seealso<< \code{\link{simDp}}, \code{\link{cWard}}.
  ##references<< B. Devlin and N. Risch. A comparison of linkage disequilibrium measures for fine-scale mapping. Genomics, 29(2):311-322, 1995.
}, ex=function(){
  N <- 100
  p <- 20
  X <- matrix(ceiling(runif(N*p)*3), N, p)
  simR2(X, 1:5, 1:10)
  simR2(X)
})

.simR2Matrix <- function(X, i, j) {
  stats <- "R.squared"
  if (missing(i) && missing(j)) {
    res <- ld(.snpM(X), depth=1, stats=stats)
    res <- slot(res, "x")
  } else if (length(i)>1 && length(j)>1) {
    res <- ld(.snpM(X[, i, drop=FALSE]), .snpM(X[, j, drop=FALSE]), depth=ncol(X)-1, stats=stats)
  } else {
    res <- ld(.snpM(X[, i, drop=FALSE]), .snpM(X[, j, drop=FALSE]), stats=stats)
  }
  res
}

simDp <- structure(function( #Calculation of the \eqn{D'} linkage disequilibrium measure
### This function return the pairs of \eqn{D'} LD measures between each column of \eqn{X[,i]} and each column of \eqn{X[,j]}. If i and j are not provided, this function returns
###a vector of the \eqn{p-1} \eqn{D'} LD measures between the adjacent columns of X.
                            X,
### genotype matrix coded in {1, 2, 3} with the SNPs in columns.
                            i,
### (optional) vector of integers. Indices of the SNPs of the first group among the columns of X of length \eqn{l_1}. \eqn{l_1} must be between 1 and p.
                            j
### (optional) vector of integers. Indices of the SNPs of the second group among the columns of X of length \eqn{l_2}. \eqn{l_2} must be between 1 and p.
                            ){

  if (class(X)!="SnpMatrix"){
    X <- .snpM(X)
  }
  stats <- "D.prime"
  if (missing(i) && missing(j)) {
    res <- ld(X, depth=1, stats=stats)
    res <- slot(res, "x")
  } else if (length(i)>1 && length(j)>1) {
    res <- ld(X[,i], X[,j], depth=ncol(X)-1, stats=stats)
  } else {
    res <- ld(X[,i], X[,j], stats=stats)
  }
  ##value<< If \code{i} and \code{j} are provided, a \code{l_1 x l_2} matrix of the \eqn{D'} LD measures is returned. Else, a (p-1) vector of the \eqn{D'} LD measures between the adjacent columns of X is returned.
  res
  ##seealso<< \code{\link{simR2}}, \code{\link{cWard}}.
  ##references<< B. Devlin and N. Risch. A comparison of linkage disequilibrium measures for fine-scale mapping. Genomics, 29(2):311-322, 1995.
}, ex=function(){
  N <- 100
  p <- 20
  X <- matrix(ceiling(runif(N*p)*3), N, p)
  simDp(X, 1:5, 1:10)
  simDp(X)
})

.simDpMatrix <- function(X, i, j) {
  stats <- "D.prime"
  if (missing(i) && missing(j)) {
    res <- ld(.snpM(X), depth=1, stats=stats)
    res <- slot(res, "x")
  }
  if (length(i)>1 && length(j)>1) {
    res <- ld(.snpM(X[, i, drop=FALSE]), .snpM(X[, j, drop=FALSE]), depth=ncol(X)-1, stats=stats)
  } else {
    res <- ld(.snpM(X[, i, drop=FALSE]), .snpM(X[, j, drop=FALSE]), stats=stats)
  }
  res
}

.snpM <- function(X){
  ## if (!inherits(Xs, "SnpMatrix")) {
  colnames(X) <- 1:ncol(X)
  rownames(X) <- 1:nrow(X)
  mode(X) <- "raw"
  X <- new("SnpMatrix", X)
  ## }
  return(X)
}

## X must be coded in 1, 2, 3
## 0 and p+1 are the buffer snps
## i and j are 2 classes
## !!!! NA values in the calculus of LD ???? !!!!
## X is an object of class SnpMatrix
.somme <- function(i, j, X, sim){
  p <- ncol(X)
  if (min(i)==0 | max(i)==(p+1) )
    {
      stop("A buffer SNP has been merged!")
    } else {
      mesure <- sim(X, i, j)
      if (length(which(is.na(mesure)))!=0){
        stop("NA similarity values are not allowed")
      }
      mu <- sum(mesure)
    }
  return(mu)
}

## h controls the maximum lag between columns of X considered. Thus, LD measures are calculated between X[,i] and X[,j] only if i and j differ by no more than h
## the ld function generates NA on real data !!!!!
.sommeh <- function(i, j, X, h, sim, flavor="v1"){
  p <- ncol(X)
  i0 <- min(i)
  i1 <- max(i)
  j0 <- min(j)
  j1 <- max(j)
  if (i0==0 | i1==(p+1)) {
    stop("A buffer SNP has been merged!")
  }
  if (i1+h < j0)  {
    somme <- 0
  } else if ( i0+h >= j1) {
    s <- sim(X, i, j)
    ind <- which(is.na(s), arr.ind=TRUE)
    browser("NA generated", expr=length(ind)!=0)
    somme <- sum(s)
  } else if (length(i)==1) {
    s <- sim(X, i, j0:min(j1, i0+h))
    ind <- which(is.na(s), arr.ind=TRUE)
    browser("NA generated", expr=length(ind)!=0)
    somme <- sum(s)
  } else if (length(j)==1) {
    s <- sim(X, j, seq(i1, max(i0, j1-h)))
    ind <- which(is.na(s), arr.ind=TRUE)
    browser("NA generated", expr=length(ind)!=0)
    somme <- sum(s)
  } else {  ## (only happens if h is small ?)
    if (flavor=="v1") {
      d1 <- max(j0-h, i0)
      f2 <- min(j1, i1+h)

      if ( (d1 == (j0-h)) && (f2 == (i1 + h)) ) {
        ## 1 triangle
        somme1 <- 0
        for (kk in 0:(i1-d1)) {
          somme1 <- somme1 + sum(sim(X, d1+kk, j0:(j0+kk)))
        }
        somme <- somme1
      } else if ( (d1 == (j0-h)) && (f2 == (j1)) ) {
        ## 1 triangle + 1 rectangle
        somme1 <- 0
        for (kk in 0:(j1-h-d1)) {
          somme1 <- somme1 +  sum(sim(X, d1+kk, j0:(j0+kk)))
          ## PN [2013-11-21]: WHY NOT USING '.somme' HERE AND LATER ???
        }
        somme <- somme1 + sum(sim(X, seq(j1-h+1, i1, by=1), j0:f2))
      } else if ( (d1 == i0) && (f2 == (i1+h)) ) {
        ## 1 triangle + 1 rectangle
        somme1 <- 0
        for (kk in 0:(i1-i0)) {
          somme1 <- somme1 + sum(sim(X, i0+kk ,(i0+h):(i0+h+kk)))
        }
        somme <- somme1 + sum(sim(X, d1:i1, seq(j0, i0+h-1, by=1)))
      } else {
        ## ( (d1 == i0) && (f2 == j1) )
        ## 1 triangle + 3 rectangles
        somme1 <- 0
        for (kk in 0:(j1-h-i0)){
          somme1 <- somme1 + sum(sim(X, i0+kk, (i0+h):(i0+h+kk)))
        }
        s1 <- sum(sim(X, seq(d1, j1-h, by=1), seq(j0, i0+h-1, by=1)))
        s2 <- sum(sim(X, seq(j1-h+1, i1, by=1), seq(j0, i0+h-1, by=1)))
        s3 <- sum(sim(X, seq(j1-h+1, i1, by=1), seq(i0+h, f2, by=1)))
        somme <-  somme1 + s1 + s2 + s3
      }
    } else if (flavor=="v2") {
      m0 <- max(j0-h, i0)
      m1 <- min(j1-h, i1)

      ## triangle
      iTri <- seq(m0, m1)
      jTri <- iTri+h
      matTri <- sim(X, iTri, jTri)
      sumTri <- sum(matTri[!upper.tri(matTri)])

      sumR <- 0
      ## rectangle #1
      if (j0-h<i0) {
        jR <- seq(j0, m0+h-1)
        matR <- sim(X, iTri, jR)
        sumR <- sumR + sum(matR)
      }

      ## rectangle #2
      if (j1-h<i1) {
        iR <- seq(m1+1, i1)
        matR <- sim(X, iR, jTri)
        sumR <- sumR + sum(matR)
      }

      ## rectangle #3
      if ((j0-h<i0) && (j1-h<i1)) {
        matR <- sim(X, iR, jR)
        sumR <- sumR + sum(matR)
      }
      somme <- sumTri+sumR
    }
  }
  return(somme)
}

## classes est une liste d'indices
.majClasses <- function(j0, classes){
  classes[[j0]] <- c(classes[[j0]], classes[[j0+1]])
  classes[[j0+1]] <- NULL
  return(classes)
}

## D: distance between classes
## Sdiag and Sdiag1 must be updated !
## classes must be updated !!
## the index in Sdiag1 is the index of the line !!
.majD <- function(j0, D, Sdiag, Sdiag1, classes){
  nM1 <- length(classes[[j0-1]])
  n0 <- length(classes[[j0]])
  n1 <- length(classes[[j0+1]])
  D[j0-1] <- -((nM1*n0)/(nM1+n0))*( (2/(nM1*n0))*Sdiag1[j0-1] - (1/(nM1^2))*Sdiag[j0-1] - (1/(n0^2))*Sdiag[j0] )
  D[j0] <-  -((n1*n0)/(n1+n0))*( (2/(n0*n1))*Sdiag1[j0] - (1/(n0^2))*Sdiag[j0] - (1/(n1^2))*Sdiag[j0+1] )
  D <- D[-(j0+1)]
  return(D)
}

## j0 is the first class (of the left)
.majLabels <- function(j0, labels, ii){
  labels[j0] <- ii
  labels <- labels[-(j0+1)]
  return(labels)
}

## Sdiag: sum of the LD measures within the class
## update before the merging of j0 and j0+1
## Sdiag1 is not updated !!!
.majSdiag <- function(j0, Sdiag, Sdiag1){
  Sdiag[j0] <- Sdiag[j0] + Sdiag[j0+1] + 2*Sdiag1[j0]
  Sdiag <- Sdiag[-(j0+1)]
  return(Sdiag)
}

## Sdiag1 : sum of the subdiagonal elements
## Sdiag1 is not updated!!!
## update before the merging of j0 and j0+1
## Remark: the same dendrogram as chclust but not the same heights !
.majSdiag1 <- function(j0, Sdiag1, classes, X, h, sim, flavor){
  l <- length(Sdiag1)
  if (j0==l) {
    stop ("A buffer SNP has been merged !")
  }
  if (j0==2) {
    S2 <- .sommeh(classes[[j0]], classes[[j0+2]], X, h, sim, flavor)
    Sdiag1[j0] <- Sdiag1[j0+1] + S2
    Sdiag1[j0-1] <- -Inf  ## not required
  } else if (j0==(l-1)) {
    S1 <- .sommeh(classes[[j0-1]], classes[[j0+1]], X, h, sim, flavor)
    Sdiag1[j0-1] <- Sdiag1[j0-1] + S1
    Sdiag1[j0] <- -Inf
  } else {
    S1 <- .sommeh(classes[[j0-1]], classes[[j0+1]], X, h, sim, flavor)
    S2 <- .sommeh(classes[[j0]], classes[[j0+2]], X, h, sim, flavor)
    Sdiag1[j0-1] <- Sdiag1[j0-1] + S1
    Sdiag1[j0] <- Sdiag1[j0+1] + S2
  }
  Sdiag1 <- Sdiag1[-(j0+1)]
  return(Sdiag1)
}
.cWardOld <- function(X, p, blMin, h, sim, flavor, trace.time){

  if (trace.time) {
    t0 <- Sys.time()
    ts <- 0
  }
  ## initialization
  step <- 0
  gains <- rep(0, p-blMin)
  classes <- as.list(0:(p+1))
  l <- length(classes) ## p + 2 buffer snps
  merge <- matrix(0, nrow=p-blMin, ncol=2)  ## matrix of the merges
  labels <- seq(-1, -(p+1), by=-1)  ## vector of the labels (same length as G)
  trW <- numeric(p-blMin-1)  ## vector of trace(W)

  if (FALSE) {  ## too slow
    sd1 <- sapply(1:(p-1), function(ii) sim(X, ii, ii+1))
  } else {
    sd1 <- sim(X)
  }

  Sdiag1 <- c(-Inf, sd1, -Inf)

  Sdiag <- rep(1, p+2)  ## p + 2 buffer snps
  D <- 1-Sdiag1
  gains[step] <- min(D)

  steps <- 1:(p-2)
  for (step in steps) {

    if (trace.time) {
      td <- as.numeric(Sys.time()-t0, units="secs")  ## time difference in seconds
      ts <- c(ts, td)
    }
    j0 <- which.min(D)
    if (j0==1 | j0==(l-1)){
      stop("the min of the distances should not be verified by the buffer snp ! ")
    }
    ## browser("NA generes dans le calcul du LD", expr=(j0==1 | j0==(l-1)))

    Sdiag <- .majSdiag(j0, Sdiag, Sdiag1)
    Sdiag1 <- .majSdiag1(j0, Sdiag1, classes, X, h, sim, flavor)

    gains[step] <- min(D)
    classes <- .majClasses(j0, classes)
    D <- .majD(j0, D, Sdiag, Sdiag1, classes)

    idxs <- c(2:(length(classes)-1))  ## indices of non-buffer SNPs
    weights <-  1/sapply(classes, length)
    trW[step] <- sum(Sdiag[idxs]*weights[idxs])

    merge[step,] <- c(labels[j0-1], labels[j0])  ## the first buffer snp !!!
    labels <- .majLabels(j0-1, labels, step)
  }

  traceW <- cbind(rev(p-steps), rev(p-trW))
  traceW1 <- p - ((sum(Sdiag)-2+2*Sdiag1[2])/p)     ## step "p-1": all the SNPs in a single class; Take care of the buffer snps + the unique term of Sdiag1 !!!
  traceW <- rbind(c(1, traceW1), traceW)

  ## merging the remaining 2 classes
  j0 <- 2
  step <- step+1
  gains[step] <- D[2]
  classes <- .majClasses(2, classes)
  merge[step, ] <-  c(labels[j0-1], labels[j0])
  labels <- .majLabels(j0-1, labels, step)

  height <- cumsum(gains)
  tree <- list(traceW=traceW,
               gains=gains,
               merge = merge,
               height = height,
               seqdist = height,
               order = 1:p,
               labels = paste("",1:p),
               method = "cWard",
               call = match.call(),
               dist.method = attr(D, "method"))
  class(tree) <- "hclust"
  if (trace.time) {
    attr(tree, "time") <- ts
  }
  return(tree)
}

## blMin: the minimal number of blocks (no clusters if length(classes)<=blMin)
## Remark: cannot use (cutree, plot...) if blMin>1
cWard <- structure(function( #Constrained Hierarchical Clustering using the LD similarity
### This function performs a constrained hierarchical clustering based on the Ward's incremental sum of squares algorithm and using the LD as a measure of similarity.
                            ##details<< \code{cWard} performs a constrained hierarchical clustering of a genotype matrix, with clusters constrained by variables order. It returns an object of class \code{\link[stats]{hclust}} which can be plotted and interrogated.
                            X,
### genotype matrix coded in 1, 2 or 3. The variables to be clustered should be in columns.
                            blMin=1,
### depth of the clustering. Default 1. See details.
                            ##details<<The standard agglomerative hierarchical clustering starts with p groups of size 1 and successively merges the pair of groups leading to the minimal increase in within-group dispersion. This merging process is repeated until \code{blMin} groups remain.
                            h=p-1,
### integer; controls the maximum lag between columns of X considered. Thus, LD measures are calculated between X[,i] and X[,j] only if i and j differ by no more than \code{h}. h must be less or equal to (p-1)
                            sim=simR2,
### a funct name specifying the LD measure to be used for the clustering. This should contain one of these funct names: simR2 or simDp. Default: simR2
                            ##details<<The LD measures are calculated using the function \code{\link[snpStats]{ld}} which handles missing values of the genotype.
                            heaps=TRUE,
### logical; if TRUE the version with binary heaps will be used. This latter is faster than the version without heaps.
                            flavor="v2",
### flavor parameter; two code versions "v1" and "v2" are implemented. "v2" is faster when h is less than the number of variables.
                            trace.time=FALSE
### logical. If TRUE, the CPU time in seconds used at each step of the clustering. Is returned as an attribute ("time") of the output list.
                            ){
  p <- ncol(X)
  stopifnot(h <= p-1)

  if (class(X)!="SnpMatrix"){
    X <- .snpM(X)
  }

  if (heaps==TRUE){
    if (identical(sim, simR2)){
      LDsim <- "R.squared"
    }
    else if (identical(sim, simDp)){
      LDsim <- "D.prime"
    }
    tree <- .cWHeaps(X=X, p=p, h=h, blMin=blMin, LDsim=LDsim, trace.time=trace.time)
  }
  else {
    tree <- .cWardOld(X=X, p=p, blMin=blMin, h=h, sim=sim, flavor=flavor, trace.time=trace.time)
  }
  ##value<< An object of class \code{hclust} which describes the tree produced by the clustering process. The \code{(p-1) x 2} matrix of the within dispersion measures \eqn{W_k} at each step k of the clustering are also added to the returned list.
  return(tree)

  ##references<< D. Clayton. snpStats: SnpMatrix and XSnpMatrix classes and methods, 2013. R package version 1.12.0.
  ##references<< J. H. Ward Jr. Hierarchical grouping to optimize an objective function. Journal of the American statistical association, 58(301):236-244, 1963.
  ##seealso<< \code{\link{gapStatistic}}, \code{\link{select}}.
}, ex=function() {

  set.seed(5)
  ## unstructured data
  N <- 100
  p <- 50
  X <- matrix(ceiling(runif(N*p)*3), N, p)
  unstruct.R2 <- cWard(X, sim=simR2)
  plot(unstruct.R2, main="Dendrogram")

  ## structured data
  N <- 100
  r2 <- 0.5
  corr <- 0.6
  p <- 120
  blockSizes <- c(1:15)
  sig.blocks <- 5
  nb.per.block <- 1
  betas <- simBeta(blockSizes, sig.blocks, nb.per.block)
  sim <- simulation(N, betas$betaMat[,"betaSNP"], betas$blockSizes, corr, r2, minMAF=0.45)
  Z <- sim$X
  struct.R2 <- cWard(Z, sim=simR2)
  struct.Dp <- cWard(Z, sim=simDp)
  par(mfrow=c(1,2))
  plot(unstruct.R2, main="Unstructured data")
  ## the underlying structure can be observed in the dendrogram
  plot(struct.R2, main="Structured data")

  groups <- betas$betaMat[,"groups"]
  nBlocks <- max(groups)
  groupsR2 <- cutree(struct.R2, nBlocks)
  names(groupsR2) <- NULL
  groupsDp <- cutree(struct.Dp, nBlocks)
  names(groupsDp) <- NULL
  identical(groups, as.numeric(groupsDp))   ## not always TRUE if corr<0.5 !
  identical(groups, as.numeric(groupsR2))  ## not always TRUE if corr<0.5 !

})


############################################################################
## HISTORY:
## 2015-06-03
## o accelerated version with pencils and heaps added
## 2014-07-16
## o BUG FIX in calculation of trace(W): buffer SNPs were included in the
## calculation.
## 2013-...
## o Created.
############################################################################

