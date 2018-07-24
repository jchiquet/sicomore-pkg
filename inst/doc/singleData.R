## ---- message=FALSE, warning=FALSE---------------------------------------
library(SIComORe)
library(Matrix)
library(mvtnorm)

## ------------------------------------------------------------------------
p <- 200  # the number of features
K <- 113  # number of groups 
## 2 relevants groups with 30 variables
## 1 irrelevant group with 30 variables
## 110 irrelevant groups with a single variable
grp.size <- c(30, 30, 30,rep(1,110))

## ------------------------------------------------------------------------
## covariance defined blockwise
n <- 100
rhos <- runif(K,.5,.95) # correlation within groups
Sigma <- bdiag(lapply(1:K, function(k) return(matrix(rhos[k],grp.size[k],grp.size[k]))))
diag(Sigma) <- 1
grp <- rep(1:length(grp.size), grp.size)
X <- scale(rmvnorm(n, sigma=as.matrix(Sigma)))

## ------------------------------------------------------------------------
image(Matrix(cor(X)))

## ------------------------------------------------------------------------
X.comp <- t(rowsum(t(X), grp)/tabulate(grp))

## ------------------------------------------------------------------------
dim.theta <- ncol(X.comp)
theta <- rep(0, dim.theta)
theta[c(1,2)] <- runif(2, min=5, max=10) ## simple effects on the first two groups

## ------------------------------------------------------------------------
sigma <- 5
epsilon <- rnorm(n) * sigma
epsilon.test <- rnorm(n) * sigma
y <- X.comp %*% theta + epsilon
r2 <- 1-sum(epsilon^2) / sum((y-mean(y))^2)
r2

## ------------------------------------------------------------------------
hierarchy <- hclust(dist(t(scale(X))), method="ward.D2")
plot(hierarchy)

## ------------------------------------------------------------------------
out.hcar  <- getHierLevel(X, y, hierarchy, choice="lambda.1se", selection="hcar")
out.mlgl  <- getHierLevel(X, y, hierarchy, choice="lambda.1se", selection="mlgl")
out.fhcar <- getHierLevel(X, y, hierarchy, choice="lambda.1se", selection="fast")

## ------------------------------------------------------------------------
library(ggplot2)
all.cv <- rbind(cbind(out.hcar$cv.error , method="hcar"),
                cbind(out.fhcar$cv.error, method="fast"),
                cbind(out.mlgl$cv.error , method="MLGL"))
print(ggplot(all.cv, aes(x=lambda, y=mean, colour=method, group=method)) + geom_smooth(aes(ymin=mean-sd, ymax=mean+sd), stat="identity") + coord_trans(x="log"))

## ------------------------------------------------------------------------
out.mlgl$getGrp()
out.hcar$getGrp()
out.fhcar$getGrp()

