library(Matrix)
library(mvtnorm)

getSigmaBlock <- function(p, sizes=rmultinom(1,p,rep(p/K,K)), rho=rep(0.75,4)) {
  K <- length(rho)
  Sigma <- bdiag(lapply(1:K, function(k) return(matrix(rho[k],sizes[k],sizes[k]))))
  diag(Sigma) <- 1
  return(Sigma)
}

rPoisLN <- function(n,Sigma) {
  p <- ncol(Sigma)
  Z <- rmvnorm(n,sigma=as.matrix(Sigma))
  return(matrix(rpois(n*p, exp(Z)), n, p))

}

simu.paper <- function(n=50, p1=100, p2=80, K1=4, K2=5, s=1, sigma=1) {

  ## FIRST DATA
  ## covariance defined blockwise
  rhos <- runif(K1,.5,.95) # correlation within groups
  grp1 <- rmultinom(1,p1,rep(p1/K1,K1)) ## group sizes
  Sigma1 <- getSigmaBlock(p1,grp1,rhos)
  grp1 <- rep(1:length(grp1), grp1)
  X1 <- scale(rmvnorm(n, sigma=as.matrix(Sigma1)))

  ## SECOND DATA
  ## covariance defined blockwise
  rhos <- runif(K2,.5,.95) # correlation within groups
  grp2 <- rmultinom(1,p2,rep(p2/K2,K2)) ## group sizes
  Sigma2 <- getSigmaBlock(p2,grp2,rhos)
  grp2 <- rep(1:length(grp2), grp2)
  X2 <- scale(rmvnorm(n, sigma=as.matrix(Sigma2)))

  inter <- rlm.inter(X1, X2, grp1, grp2, s, sigma)

  return(list(y=inter$y, X1=X1, X2=X2, grp1=grp1, grp2=grp2,
              theta.inter.full = inter$theta.full, r2=inter$r2))
}

rlm.inter <- function(X1, X2, grp1, grp2, s, sigma) {

  n <- nrow(X1)
  simple.effects.1 <- computeCompressedDataFrame(X1, grp1, compression = "mean")
  simple.effects.2 <- computeCompressedDataFrame(X2, grp2, compression = "mean")

  interactions <- getPhiCompressed(simple.effects.1,simple.effects.2)

  ## vector of true parameter at the correct level of the hierarchies
  dim1 <- ncol(simple.effects.1)
  dim2 <- ncol(simple.effects.2)
  dim.inter <- dim1 * dim2

  ind.inter <- sample(1:dim.inter, s)
  theta.inter <- rep(0, dim.inter)
  theta.inter[ind.inter] <- runif(s, min=2, max=4) ## interactions : sparse
  theta.block.inter <- Matrix(theta.inter, dim1, dim2)

  ind.simple.effects <- which(theta.block.inter!=0, arr.ind=TRUE)
  theta1 <- rep(0, dim1)
  theta1[ind.simple.effects[,1]] <- runif(nrow(ind.simple.effects), min=1,max=2)

  theta2 <- rep(0, dim2)
  theta2[ind.simple.effects[,2]] <- runif(nrow(ind.simple.effects), min=1,max=2)

  theta.full <- theta.block.inter[rep(1:nrow(theta.block.inter), table(grp1)),
                                  rep(1:ncol(theta.block.inter), table(grp2))]

  epsilon <- rnorm(n) * sigma
  mu <- 3
  y <- mu + interactions %*% theta.inter + simple.effects.1 %*% theta1 + simple.effects.2 %*% theta2 + epsilon
  r2 <- 1-sum(epsilon^2) / sum((y-mean(y))^2)

  return(list(y=y, theta.inter.full = theta.full, theta1 = theta1, theta1 = theta2, r2=r2,
              simple.effects.1=simple.effects.1, simple.effects.2=simple.effects.2))
}
