#' @title sicomore
#' @description Selection of Interaction effects in COmpressed  Multiple Omics REpresentation
#'
#' From a set of input matrices and phenotype related to the same set of individual, sicomore is a two-step method which
#' 1. find and select groups of correlated variables in each input matrix which are good predictors for the common phenotype
#' 2. find the most predictive interaction effects between the set of data by testing for interaction between the selected groups of each input matrix
#'
#' @param y response variable (phenotype)
#' @param X.list a list of input matrices. Must have the same number of rows (number of observations). May have different number of columns (predictors)
#' @param compressions a vector of string for the compression methods for each data sets. Default used the mean to compressed predictors in the same group.
#' @param method.clus a vector of string specifying the method used for the hierarchical clustering, one for each input matrix in X.list.
#' By default, the hierarchy is obtain by a WARD clustering on the scaled input matrix.
#' To use an SNP-specific spatially contrained hierarchical clustering \insertCite{dehman}{SIComORe} from package adjclust, specify "snpClust".
#' It is also possible to specify no hierarchy with "noclust".
#' @param selection a string for the method used for variable selection for each data set.
#' Specify "sicomore" to use the method specifically developed for the package, "hcar" to use the method developped by \insertCite{park;textual}{SIComORe} or
#' "mlgl" for the method of \insertCite{grimonprez;textual}{SIComORe}.
#' @param cuts a list of numeric vector defining the cut levels to be considered for each hierarchy. By default a sequence of 100 levels is used.
#' @param choice a string (either "lambda.min" or "lambda.1se"). Indicates how the tuning parameter is chosen in the penalized regression approach.
#' @param depth.cut a vector of integers specifying the depth of the search space for the variable selection part of the algorithm.
#' This argument allows to increase the speed of the algorithm by restraining the search space without affecting too much the performance.
#' A value between 3 and 6 is recommended, the smaller the faster.
#' @param taxonomy a hierarchical tree object constructed using taxonomical unit.
#' This argument could be useful if one want to provid a specific hierarchical tree based on taxonomy rather than euclidean distance.
#' Recommended for those who want to detect genomic-metagenomic interactions.
#' @param mc.cores an integer for the number of cores to use in the parallelization of the cross-validation and some other functions.
#'
#' @return a RC object with class \code{sicomore} with methods \code{plot()}, \code{getSignificance()} and the following fields
#' \itemize{
#'  \item{pval:}{A matrix of p-values for each interactions effects between the compressed variables originating from 2 input matrices.}
#'  \item{pval.beta1:}{A matrix of p-values for the corresponding main effects of the compressed variables originating from the first input matrix}
#'  \item{pval.beta2:}{A matrix of p-values for the corresponding main effects of the compressed variables originating from thesecond input matrix}
#'  \item{models}{A class'sicomore-model' RC object obtain from \code{getHierLevel} function and with methods \code{nGrp()}, \code{nVar()}, \code{getGrp()}, \code{getVar()},
#'  \code{getCV()}, \code{getX.comp()}, \code{getCoef()}}
#'  \item{tuplets:}{A list of integer vector specifying the indexes of the compressed variables which are fitted together in a linear model with interaction.}
#'  \item{dim:} A array of integers specifying the dimension of the compressed matrices.
#' }
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{}
#'
sicomore <- function(y, X.list,
                     compressions = rep("mean", length(X.list)),
                     selection    =  c("sicomore", "hcar", "mlgl"),
                     cuts = lapply(X.list, function(X) unique(round(10^seq(log10(ncol(X)), log10(2), len=100)))),
                     choice=c("lambda.min", "lambda.min"), method.clus = c("ward.D2","ward.D2"),
                     depth.cut = c(3,3), mc.cores = NULL, taxonomy = NULL, verbose = TRUE)
{

  selection <- match.arg(selection, c("sicomore", "hcar", "mlgl"))
  ndata <- length(X.list)
  models <- vector("list", ndata)
  hierarchies  <- vector("list" , length(X.list))
  for (i in seq_along(X.list)) {
    if (verbose == TRUE) cat("\nConsidering hierarchy",i, " - ", selection, "selection method.")

    if (method.clus[i] == "ward.D2"){
      hierarchies[[i]] <- hclust(dist(t(scale(X.list[[i]]))), method="ward.D2")
      models[[i]] <- getHierLevel(X.list[[i]], y, hierarchies[[i]], cut.levels = cuts[[i]], compression=compressions[i],
                                  selection=selection, choice=choice[i], depth.cut = depth.cut[i], mc.cores=mc.cores)
    }
    if (method.clus[i] == "snpClust") {
      if (ncol(X.list[[i]]) > 600) h <- 600
      else h <- ncol(X.list[[i]])-1
      hierarchies[[i]] <- snpClust(X.list[[i]], h=h)
      models[[i]] <- getHierLevel(X.list[[i]], y, hierarchies[[i]], cut.levels = cuts[[i]], compression=compressions[i],
                                  selection=selection, choice=choice[i], depth.cut = depth.cut[i], mc.cores=mc.cores)
    }
    if (method.clus[i] == "noclust"){
      if (!is.null(mc.cores)) doMC::registerDoMC(cores=mc.cores)
      cv <- cv.glmnet(X.list[[i]], y, parallel=!is.null(mc.cores))
      groups <- predict(cv, s=choice[i], type="nonzero")[[1]]
      coefficients <- predict(cv, s=choice[i], type="coefficients")[-1] # without intercept
      coefficients <- coefficients[groups]
      models[[i]] <- new("sicomore-model", coefficients = coefficients, groups = as.list(groups),
                         X.comp = X.list[[i]][,groups], cv.error = data.frame(mean = cv$cvm, sd = cv$cvsd, lambda  = cv$lambda))
    }
    if (method.clus[i] == "taxonomy"){
      if (is.null(taxonomy)) stop("Please provide taxonomy for Metagenomic data", call. = F)

      taxon <- colnames(taxonomy)[-1]
      Xcomp <- list()
      group.name <- list()
      for (t in seq(taxon)){
        tax <- taxonomy[order(taxonomy[,t+1]),t+1]
        group <- unlist(sapply(1:nlevels(as.factor(tax)), function(l) rep(l, length(which(tax== unique(tax)[l])))))
        group <- c(group, rep(max(group)+1, length(which(is.na(tax)))))
        name <- taxonomy[order(taxonomy[,t+1]),1]
        Xcomp[[t]] <- computeCompressedDataFrame(X.list[[i]][,name], group, compressions[i])
        colnames(Xcomp[[t]]) <- unique(tax)
        group.name[[t]] <- split(taxonomy[order(taxonomy[,t+1]),1], group)
        names(group.name[[t]]) <- unique(tax)
      }

      Xcomp <- cbind(do.call(cbind, Xcomp), X.list[[i]])
      group.name$species <- split(colnames(X.list[[i]]), 1:ncol(X.list[[i]]))
      names(group.name) <- c(colnames(taxonomy)[-1], "species")
      group.name <- do.call(c,group.name)

      if (!is.null(mc.cores)) doMC::registerDoMC(cores=mc.cores)
      sel <- NULL
      while (length(sel) == 0){
        cv <- cv.glmnet(Xcomp, y, parallel=!is.null(mc.cores))
        sel <- predict(cv, s=choice[i], type="nonzero")[[1]]
        group.select <- group.name[sel]
        sel <- sel[-which(duplicated(group.select))]  ## Remove duplicated species
        group.select[which(duplicated(group.select))] <- NULL
      }
      ## Get coefficients
      coefficients <- predict(cv, s=choice[i], type="coefficients")[-1] # without intercept
      coefficients <- coefficients[sel]

      ## If no main effets are selected, keep the single variables with no compression
      if (length(sel) == 0) {
        warning("All coeffs are null, single variables are kept as main effects", call.=F, immediate.=TRUE)
        X.comp <- X.list[[i]]
        group.select <- as.list(1:ncol(X.list[[i]]))
      }
      else X.comp <-  Xcomp[,sel]

      models[[i]] <- new("sicomore-model", coefficients = coefficients, groups = group.select, X.comp = X.comp,
                         cv.error = data.frame(mean = cv$cvm, sd = cv$cvsd, lambda  = cv$lambda))
    }
  }

  sequences <- lapply(models, function(hier) 1:hier$nGrp())
  tuplets   <- as.list(data.frame(t(expand.grid(sequences, KEEP.OUT.ATTRS = FALSE))))  # get all the combinaisons

  pval <- sapply(tuplets, function(tuplet) {
    tuplet <- unname(unlist(tuplet))
    data.tmp <- as.data.frame(cbind(sapply(1:ndata, function(m) models[[m]]$getX.comp()[,tuplet[m]]), y))
    colnames(data.tmp) <- c("X1", "X2", "y")
    pval <- summary(lm(y ~ X1*X2, data=data.tmp))$coefficients[-1,4]
    if (!("X1:X2" %in% names(pval))) pval<- c(pval,1) ; names(pval) <- c("X1", "X2", "X1:X2")
    return(pval)
  })

  ## TODO : If no simple effects are reported as significant,
  ## infer the cut level on the predictor only by hierarchical clustering

  ## convert interaction pval to an array
  dims <- sapply(models, function(model) model$nGrp())
  pval.Array <- array(NA, dim = dims)
  rownames(pval.Array) <- paste0("X1comp_",1:nrow(pval.Array))
  colnames(pval.Array) <- paste0("X2comp_",1:ncol(pval.Array))

  pval.beta1 <- pval.Array ; pval.beta2 <- pval.Array ; pval.inter <- pval.Array
  pval.beta1[do.call(rbind, tuplets)] <- pval[1,]
  pval.beta2[do.call(rbind, tuplets)] <- pval[2,]
  pval.inter[do.call(rbind, tuplets)] <- pval[3,]

  return(new("sicomore-fit", pval = pval.inter, pval.beta1 = pval.beta1, pval.beta2 = t(pval.beta2), tuplets = tuplets, models=models, dim=sapply(X.list, ncol)))
}

