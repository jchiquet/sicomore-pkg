computeCompressedDataFrame <- function(X, group, compression="mean"){
  if (compression == "mean"){
    X.h <- t(rowsum(t(X), group)/tabulate(group))
  } else if (compression=="SNP.dist"){
    ## Calculate median haplotype for each column of Genotype
    haplo.median <- apply(X,2,median)
    #X.h <-   t(rowsum(abs(t(X)-haplo.median), group)/unlist(sapply(1:max(group), function(i) length(group[which(group==i)]))))
    X.h <-   t(rowsum(abs(t(X)-haplo.median), group)/tabulate(group))
    ### Is this meaningful to normalize with tabulate ????
  } else if (compression == "sum"){
    X.h <- t(rowsum(t(X), group)/unlist(sapply(1:max(group), function(i) length(group[which(group==i)]))))
  } else stop("No valid compression method")
}

getPhiCompressed <- function(X1k, X2l) {
  return(do.call(cbind, lapply(1:ncol(X2l), function(i) sweep(X1k, 1, X2l[, i], "*"))))
}

computeCompressedDataFrameFromVariables <- function(X, variables, compression="mean"){

  if (length(variables)==1) Xvar<-cbind(X[,variables])
  else Xvar<-X[,variables]

  if (compression == "mean"){
    X.h <- cbind(as.numeric(rowMeans(Xvar)))
  } else if (compression=="SNP.dist"){
    ## Calculate median haplotype for each column of Genotype
    haplo.median <- apply(Xvar,2,median)
    X.h <-  cbind(as.numeric(rowMeans(abs((Xvar)-haplo.median))))
    ### Is this meaningful to normalize with tabulate ????
  } else if (compression == "sum"){
    X.h <- cbind(as.numeric(rowSums(Xvar)/length(variables)))
  } else stop("No valid compression method")
}

computePval<-function(X,y,threshold=1){
  # Compute a pvalue for each column of X
  pval.raw <-apply(X,2,function(x){fstat<-summary(lm(y~x+0))$fstatistic;  return(1-pf(fstat[1],fstat[2],fstat[3]))})
  # Count the number of significative pvalues
  pval.adjusted <- -log10(p.adjust(pval.raw,"BH"))
  count <- sum( pval.adjusted > threshold)
  return(list(pval.raw=pval.raw, pval.adjusted = pval.adjusted,   count=count))
}
