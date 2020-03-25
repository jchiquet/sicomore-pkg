#' sicomore
#'
#' Selection of Interaction effects in COmpressed  Multiple Omics REpresentation
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
#' To use an SNP-specific spatially contrained hierarchical clustering \insertCite{dehman}{sicomore} from package adjclust, specify "snpClust".
#' It is also possible to specify no hierarchy with "noclust".
#' @param selection method used to perform variable selection. Either 'sicomore', 'mlgl' or 'rho-sicomore' (see details). Default is 'sicomore'.
#' @param cuts a list of numeric vector defining the cut levels to be considered for each hierarchy. By default a sequence of 100 levels is used.
#' @param choice a string (either "lambda.min" or "lambda.1se"). Indicates how the tuning parameter is chosen in the penalized regression approach.
#' @param depth.cut a vector of integers specifying the depth of the search space for the variable selection part of the algorithm.
#' This argument allows to increase the speed of the algorithm by restraining the search space without affecting too much the performance.
#' A value between 3 and 6 is recommended, the smaller the faster.
#' @param taxonomy a hierarchical tree object constructed using taxonomical unit.
#' This argument could be useful if one want to provid a specific hierarchical tree based on taxonomy rather than euclidean distance.
#' Recommended for those who want to detect genomic-metagenomic interactions.
#' @param mc.cores an integer for the number of cores to use in the parallelization of the cross-validation and some other functions. Default is 1.
#' @param verbose not yet documented
#' @param stab A boolean indicating if the algorithm perform a lasso stability selection using stabsel function from stabs package.
#' @param stab.param A list of parameter for the stabsel function if stab = TRUE.
#' The parameters to choose are the FWER (1 by default), cut-off (0.75 by default) and bootstrap number (200 by default).
#' @param grp.min Minimum number of groups to consider for the highest level in the hierarchy.
#' Correspond to the highest allowed cut in the hierarchy. If NULL, no restriction is given.
#' @details The methods for variable selection are variants of Lasso or group-Lasso designed to perform selection of interaction between multiple hierarchies:
#' 'sicomore' and 'rho-sicomore' (see \insertCite{sicomore;textual}{sicomore}) use a LASSO penalty on compressed groups of variables along the hierarchies to select interactions.
#' rho-sicomore is a variant where a more sound weighting scheme is used dependending on the level of the hierarchy considered. The method 'mlgl' of \insertCite{grimonprez_PhD;textual}{sicomore}
#' uses a group-Lasso penalty which does not require compression but requires heavier computational resources.

#' @return an RC object with class \code{sicomore} with methods \code{plot()}, \code{getSignificance()} and the following fields
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
sicomore <- function(y,
                     X.list,
                     compressions = rep("mean", length(X.list)),
                     selection =  c("rho-sicomore", "sicomore", "mlgl"),
                     choice=c("lambda.min", "lambda.min"),
                     method.clus = c("ward.D2","ward.D2"),
                     depth.cut = c(3,3),
                     grp.min = NULL,
                     mc.cores = 1,
                     taxonomy = NULL,
                     verbose = TRUE,
                     stab = TRUE,
                     stab.param  = list(B = c(100,100), cutoff = c(.75,.75), PFER = c(1,1)))
{

  selection <- match.arg(selection, c("rho-sicomore", "sicomore", "mlgl"))
  ndata <- length(X.list)
  models      <- vector("list", ndata)
  hierarchies <- vector("list" , ndata)
  for (i in seq_along(X.list)) {
    if (verbose == TRUE) cat("\nConsidering hierarchy",i, " - ", selection, "selection method.")

    if (method.clus[i] == "ward.D2"){
      hierarchies[[i]] <- hclust(dist(t(scale(X.list[[i]]))), method="ward.D2")
      models[[i]] <- getHierLevel(X = X.list[[i]], y,
                                  hc.object = hierarchies[[i]],
                                  compression=compressions[i],
                                  selection=selection,
                                  choice=choice[i],
                                  depth.cut = depth.cut[i],
                                  mc.cores=mc.cores,
                                  stab, stab.param = lapply(stab.param, function(x) x[[i]]))
    }
    if (method.clus[i] == "snpClust") {
      if (ncol(X.list[[i]]) > 600) h <- 600
      else h <- ncol(X.list[[i]]) - 1
      #hierarchies[[i]] <- adjclust::snpClust(X.list[[i]], h=h)
      hierarchies[[i]] <- cWard(X.list[[i]], h)

      models[[i]] <- getHierLevel(X = X.list[[i]], y = y,
                                  hc.object = hierarchies[[i]],
                                  compression=compressions[i],
                                  selection=selection, choice=choice[i],
                                  depth.cut = depth.cut[i], mc.cores=mc.cores,
                                  stab = stab, stab.param = lapply(stab.param, function(x) x[[i]]))

      sapply(models[[i]]$getGrp(), length)
    }
    if (method.clus[i] == "noclust"){
      if (stab == TRUE){
        stab.lasso <- stabs::stabsel(x = X.list[[i]], y = y, B = 50,
                              fitfun = stabs::glmnet.lasso, cutoff = 0.75,
                              PFER = 1, mc.cores = mc.cores)
        groups <- as.numeric(stab.lasso$selected)
        models[[i]] <- new("sicomore-model",
                           groups = as.list(groups),
                           X.comp = X.list[[i]][,groups],
                           probs = stab.lasso$max)
      } else{
        cv <- cv.glmnet(X.list[[i]], y, parallel=!is.null(mc.cores))
        groups <- predict(cv, s=choice[i], type="nonzero")[[1]]
        coefficients <- predict(cv, s=choice[i], type="coefficients")[-1] # without intercept
        coefficients <- coefficients[groups]
        models[[i]] <- new("sicomore-model",
                           coefficients = coefficients,
                           groups = as.list(groups),
                           X.comp = X.list[[i]][,groups],
                           cv.error = data.frame(mean = cv$cvm, sd = cv$cvsd, lambda  = cv$lambda))
      }
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

      sel <- NULL
      while (length(sel) == 0){
        if(stab == TRUE){
          stab.lasso <- stabs::stabsel(x = Xcomp, y = y, B = 50,
                                fitfun = glmnet.lasso, cutoff = 0.75,
                                PFER = 1, mc.cores = mc.cores)
          sel <- as.numeric(stab.lasso$selected)
          group.select <- group.name[sel]
          sel <- sel[-which(duplicated(group.select))]  ## Remove duplicated species
          group.select[which(duplicated(group.select))] <- NULL
          probs <- stab.lasso$max
        } else {
          cv <- cv.glmnet(Xcomp, y, parallel=!is.null(mc.cores))
          sel <- predict(cv, s=choice[i], type="nonzero")[[1]]
          group.select <- group.name[sel]
          sel <- sel[-which(duplicated(group.select))]  ## Remove duplicated species
          group.select[which(duplicated(group.select))] <- NULL
          ## Get coefficients
          coefficients <- predict(cv, s=choice[i], type="coefficients")[-1] # without intercept
          coefficients <- coefficients[sel]

          models[[i]] <- new("sicomore-model", coefficients = coefficients, groups = group.select, X.comp = X.comp,
                             cv.error = data.frame(mean = cv$cvm, sd = cv$cvsd, lambda  = cv$lambda))
        }
      }

      ## If no main effets are selected, keep the single variables with no compression
      if (length(sel) == 0) {
        warning("All coeffs are null, single variables are kept as main effects", call.=F, immediate.=TRUE)
        X.comp <- X.list[[i]]
        group.select <- as.list(1:ncol(X.list[[i]]))
        probs <- NULL
        models[[i]] <- new("sicomore-model",
                           groups = group.select,
                           probs = probs)
      }
      else X.comp <-  Xcomp[,sel]

      models[[i]]$X.comp <- X.comp
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
    pval
  })

  ## TODO : If no simple effects are reported as significant,
  ## infer the cut level on the predictor only by hierarchical clustering

  ## convert interaction pval to an array
  dims <- sapply(models, function(model) model$nGrp())
  pval.Array <- array(NA, dim = dims)
  rownames(pval.Array) <- paste0("X1comp_", 1:nrow(pval.Array))
  colnames(pval.Array) <- paste0("X2comp_", 1:ncol(pval.Array))

  pval.beta1 <- pval.Array ; pval.beta2 <- pval.Array ; pval.inter <- pval.Array
  pval.beta1[do.call(rbind, tuplets)] <- pval[1,]
  pval.beta2[do.call(rbind, tuplets)] <- pval[2,]
  pval.inter[do.call(rbind, tuplets)] <- pval[3,]

  new("sicomore-fit", pval = pval.inter, pval.beta1 = pval.beta1, pval.beta2 = t(pval.beta2),
             tuplets = tuplets, models=models, dim=sapply(X.list, ncol))
}

