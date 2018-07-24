## ---- message=FALSE, warning=FALSE---------------------------------------
set.seed(1234)
library(SIComORe)
library(Matrix)
library(mvtnorm)

## ------------------------------------------------------------------------
n <- 100

## ------------------------------------------------------------------------
## FIRST DATA
p1 <- 50 # number of variables 
K1 <- 5  # number of groups
rhos <- runif(K1,.5,.95) # correlation within groups
grp1.size <- rmultinom(1,p1,rep(p1/K1,K1)) ## group sizes
Sigma1 <- bdiag(lapply(1:K1, function(k) return(matrix(rhos[k],grp1.size[k],grp1.size[k]))))
diag(Sigma1) <- 1
grp1 <- rep(1:length(grp1.size), grp1.size)
X1 <- scale(rmvnorm(n, sigma=as.matrix(Sigma1)))

## ------------------------------------------------------------------------
p2 <- 40
K2 <- 10
rhos <- runif(K2,.5,.95) # correlation within groups
grp2.size <- rmultinom(1,p2,rep(p2/K2,K2)) ## group sizes
Sigma2 <- bdiag(lapply(1:K2, function(k) return(matrix(rhos[k],grp2.size[k],grp2.size[k]))))
diag(Sigma2) <- 1
grp2 <- rep(1:length(grp2.size), grp2.size)
X2 <- scale(rmvnorm(n, sigma=as.matrix(Sigma2)))

## ------------------------------------------------------------------------
simple.effects.1 <- t(rowsum(t(X1), grp1)/tabulate(grp1))
simple.effects.2 <- t(rowsum(t(X2), grp2)/tabulate(grp2))

## ------------------------------------------------------------------------
interactions <- do.call(cbind, lapply(1:ncol(simple.effects.2), 
      function(i) sweep(simple.effects.1, 1, simple.effects.2[, i], "*")))

## ------------------------------------------------------------------------
s <- 2 ## number of non-null interactions
# vector of true parameter at the correct level of the hierarchies
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

## ------------------------------------------------------------------------
sigma <- 3
epsilon <- rnorm(n) * sigma
mu <- 3
y <- mu + interactions %*% theta.inter + simple.effects.1 %*% theta1 + simple.effects.2 %*% theta2 + epsilon
r2 <- 1-sum(epsilon^2) / sum((y-mean(y))^2)
r2

## ------------------------------------------------------------------------
res <- sicomore(y, list(X1, X2))

## ------------------------------------------------------------------------
theta.full <- as.matrix(theta.block.inter[rep(1:nrow(theta.block.inter), table(grp1)),
                                  rep(1:ncol(theta.block.inter), table(grp2))])

par(mfrow=c(1,2))
res$plot(main="estimated")
image(theta.full, col=terrain.colors(3, alpha = 0.4),main="TRUE",axes=FALSE)
title(outer=TRUE, main="\nPerformance for interaction recovery", )

## ------------------------------------------------------------------------
getSelPerf <- function(theta, theta.star) {
  ones.true <- which(theta.star != 0)
  zero.true <- which(theta.star == 0)

  ones.hat <- which(theta != 0)
  zero.hat <- which(theta == 0)

  tp <- sum(ones.hat %in% ones.true)
  tn <- sum(zero.hat %in% zero.true)

  fp <- sum(ones.hat %in% zero.true)
  fn <- sum(zero.hat %in% ones.true)

  if ((tp+fn)==0) recall    <- 0 else recall    <- tp/(tp+fn) 
  if ((fp+tn)==0) fallout   <- 0 else fallout   <- fp/(fp+tn) 
  if ((tp+fp)==0) precision <- 0 else precision <- tp/(tp+fp) 

  return(c(recall=recall, fallout=fallout, precision=precision))
}

## ------------------------------------------------------------------------
fast <- sicomore(y, list(X1, X2), selection="fast")
mlgl <- sicomore(y, list(X1, X2), selection="mlgl")
hcar <- sicomore(y, list(X1, X2), selection="hcar")
  
theta.mlgl <- mlgl$getSignificance()
theta.hcar <- hcar$getSignificance()
theta.fast <- fast$getSignificance()

seq.thres <- seq(0,1,len=100)
res.fast <- data.frame(t(sapply(seq.thres, function(threshold) getSelPerf(theta.fast>=threshold, theta.full))))
res.mlgl <- data.frame(t(sapply(seq.thres, function(threshold) getSelPerf(theta.mlgl>=threshold, theta.full))))
res.hcar <- data.frame(t(sapply(seq.thres, function(threshold) getSelPerf(theta.hcar>=threshold, theta.full))))

par(mfrow=c(1,1))
plot(res.fast$fallout, res.fast$recall, type="l")
lines(res.mlgl$fallout, res.mlgl$recall,col=2)
lines(res.hcar$fallout, res.hcar$recall,col=3)
legend("bottomright", legend = c("fast", "MLGL", "HCAR"), col=1:3, lty=1)

