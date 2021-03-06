## x : the LD matrix
.toMatLeft <- function(x, p, h) {
  #x <- slot(t(x), "x") # we now send a symmetric matrix from LD, so no transposition is needed
  x <- slot(x, "x")
  stopifnot(length(x)==(p-1)*h-h*(h-1)/2)
  x1 <- head(x, (p-1)*h-h*(h-1))
  m1 <- matrix(x1, ncol=h, byrow=TRUE)
  x2 <- tail(x, h*(h-1)/2)
  m2t <- matrix(nrow = h, ncol = h-1)
  idxs <- seq(along=x2)
  for (jj in 1:(h-1)) {
    idxsJJ <- head(idxs, h-jj)
    idxs <- tail(idxs, -(h-jj))
    m2t[1:(h-jj), jj] <- x2[idxsJJ]
  }
  mat <- rbind(m1, t(m2t), NA)
  mat[is.na(mat)] <- 0
  mat
}

## x : the LD matrix
.toMatRight <- function(x, p, h){
  x <- slot(x, "x")
  stopifnot(length(x)==(p-1)*h-h*(h-1)/2)
  x1 <- tail(x, (p-1)*h-h*(h-1))
  m1 <- matrix(x1, ncol=h, byrow=TRUE)
  x2 <- head(x, h*(h-1)/2)
  m2t <- matrix(nrow=h-1, ncol=h)
  idxs <- seq(along=x2)
  for (jj in (h-1):1) {
    idxsJJ <- tail(idxs, jj)
    idxs <- head(idxs, -jj)
    m2t[jj, (h-jj+1):h] <- x2[idxsJJ]
  }
  mat <- rbind(NA, m2t, m1)
  mat[is.na(mat)] <- 0
  mat
}

.rotate <- function(mat) { t(mat[nrow(mat):1,,drop=FALSE]) }

## checked !! (with PROTECT)
.buildHeap <- function(heap, D, l){
  for (ii in floor(l/2) : 1){
    heap <- .Call("percDown", heap, D, as.integer(l), as.integer(ii), PACKAGE="sicomore")
  }
  heap
}
