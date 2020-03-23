##################################################
## new algorithm of cWard with h and binary heaps
##################################################
## res$merge identical to that of the old version !
###########################################
## OK for traceW : except for traceW[1,]
###########################################

.cWHeaps <- function(X, p, h, blMin, LDsim, trace.time){
  LD <- ld(X, depth=h, stats=LDsim, symmetric = TRUE)
  .cWLD(LD, p=p, h=h, blMin=blMin, trace.time=trace.time)
}

.cWLD <- function(simMat, p, h, blMin=1, trace.time=FALSE){
  ## sum of the "rectangles" beginning from the left
  matLD <- .toMatLeft(simMat, p, h)  ## a matrix p x h (with zeros at the bottom) of the LD values
  rCum <- rowCumsums(matLD)  ## p x h matrix
  rcCumLeft <- colCumsums(rCum)  ## p x h matrix

  ## sum of the "rectangles" beginning from the right
  matLDright <- .toMatRight(simMat, p, h)
  rotatedMat <- .rotate(.rotate(matLDright))
  rCumRight <- rowCumsums(rotatedMat)  ## p x h matrix
  rcCumRight <- colCumsums(rCumRight)  ## p x h matrix

  rm(simMat)

  if (trace.time) {
    t0 <- Sys.time()
    ts <- 0
  }

  ## initialization
  gains <- rep(0, p-blMin)
  merge <- matrix(0, nrow=p-blMin, ncol=2)  ## matrix of the merges
  traceW <- matrix(0, nrow=p-blMin, ncol=2)  ## matrix of traceW
  sd1 <- matLD[1:(p-1),1]

  ## initialization of the heap
  heap <- as.integer(rep(-1, 3*p))
  lHeap <- length(heap)
  heap[1:(p-1)] <- 1:(p-1)
  D <- rep(-1, 3*p)
  D[1:(p-1)] <- 1-sd1
  ## initialization of the length of the Heap
  lHeap <- p-1
  ## each element contains a vector: c(cl1, cl2, label1, label2, posL, posR, valid)
  chainedL <- matrix(-1, nrow=12, ncol=3*p)
  rownames(chainedL) <- c("minCl1", "maxCl1", "minCl2", "maxCl2", "lab1", "lab2", "posL", "posR", "Sii", "Sjj", "Sij", "valid")
  v <- 1:(p-1)
  w <- as.integer(v+1)
  chainedL[1,v] <- v
  chainedL[2,v] <- v
  chainedL[3,v] <- w
  chainedL[4,v] <- w
  chainedL[5,v] <- -v
  chainedL[6,v] <- -w
  chainedL[7,v] <- v-1
  chainedL[8,v] <- w
  chainedL[9,v] <- 1
  chainedL[10,v] <- 1
  chainedL[11,v] <- sd1
  chainedL[12,v] <- 1
  chainedL[7,1] <- -1
  chainedL[8,p-1] <- -1
  heap <- .buildHeap(heap, D, lHeap)

  res <- .Call("cWardHeaps", rcCumRight, rcCumLeft, as.integer(h), as.integer(p), chainedL, heap, D, as.integer(lHeap), merge, gains, traceW, as.integer(blMin), PACKAGE="sicomore")

  height <- cumsum(gains)
  tree <- list(traceW=traceW,
               gains=gains,
               merge = res,
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
