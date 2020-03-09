#' @title getHierLevel
#' @description A function to select groups of variables which are good for predicting a given phenotype.
#' The groups considered corresponds various cut levels in a user defined hierarchy.
#' The selection is performed by various penalty-based regression methods.
#'
#' @param X input matrix
#' @param y response variable
#' @param hierarchy the results of a hierarchical clustering algorithm, typically \code{hclust}. Must be an "hclust" object
#' @param selection a string for the method used for variable selection for each data set.
#' Specify "sicomore" to use the method specifically developed for the package, "hcar" to use the method developped by \insertCite{park;textual}{SIComORe} or
#' "mlgl" for the method of \insertCite{grimonprez;textual}{SIComORe}.
#' @param compression a string (either "mean" or "SNP.dist"). Indicates how the groups of variables are compressed before variable selection at each level of the hierarchy.
#' @param cut.levels a numeric vector, the level consider in the hierarchy. By default a sequence of 100 levels is used.
#' @param choice a string (either "lambda.min" or "lambda.1se"). Indicates how the tuning parameter is chosen in the penalized regression approach
#' @param depth.cut an integer specifying the depth of the search space for the variable selection part of the algorithm.
#' This argument allows to increase the speed of the algorithm by restraining the search space without affecting too much the performance.
#' A value between 3 and 6 is recommended, the smaller the faster.
#' @param mc.cores an integer for the number of cores to use in the parallelization of the cross-validation and some other functions.
#'
#' @return a RC object with class 'sicomore-model', with methods \code{nGrp()}, \code{nVar()}, \code{getGrp()}, \code{getVar()}, \code{getCV()}, \code{getX.comp()}, \code{getCoef()} and with the following fields:
#' \itemize{
#'  \item{groups:}{a list with the selected groups of predictors}
#'  \item{coefficients:}{a vector with the estimated coefficients (one per selected group)}
#'  \item{X.comp:}{The compressed version of the original input matrix (as many columns as number of selected groups)}
#'  \item{cv.error:}{for the best grouping, a data frame showing the cross-validation error used in the variable selection procedure}
#'  \item{selection:}{the selection method used}
#'  \item{compression:}{the compression method used}
#' }
#' @include utils.R
#' @import glmnet MLGL
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{}
getHierLevel <- function(X,
                         y,
                         hierarchy,
                         selection = c("rho-sicomore", "sicomore", "mlgl"),
                         compression = "mean",
                         depth.cut = 3,
                         cut.levels  = unique(round(10^seq(log10(ncol(X)), log10(2), len=100))),
                         choice=c("lambda.min", "lambda.1se"),
                         mc.cores=NULL,
                         stab = FALSE,
                         stab.param = NULL) {

  ## _________________________________________________________________
  ##
  ## PRETREATMENT
  choice <- match.arg(choice)

  ## forcing the matrix type for glmnet
  if (!is.matrix(X)) X <- as.matrix(X)

  if (class(hierarchy) == "chac") class(hierarchy) <- "hclust"
  stopifnot(class(hierarchy) == "hclust")

  if (is.element(1,cut.levels)) {
    message("A cut level of 1 is not allowed: removing it from the list.")
    cut.levels <- cut.levels[!(cut.levels %in% 1)]
  }

  ## _________________________________________________________________
  ##
  ## GET THE BEST COMPRESSION LEVEL IN THE HIERARCHY
  .getHier <- switch(selection,
                     "rho-sicomore" = getHierLevel.rhosicomore,
                     "sicomore" = getHierLevel.sicomore,
                     "mlgl" = getHierLevel.MLGL)

  out <- .getHier(X, y, hierarchy, compression, cut.levels, choice, depth.cut, mc.cores, stab, stab.param)

  if (length(out$groups) > 0) {
    grouping <- rep(1:length(out$groups), sapply(out$groups, length))
    X.comp   <- cbind(computeCompressedDataFrame(X[, unlist(out$groups)], grouping, compression))
    if(!stab) out$coefficients <- setNames(out$coefficients, paste("group", 1:length(out$groups)))
    out$groups       <- setNames(out$groups, NULL)
  } else {
    grouping <- cutree(hierarchy, k=1+which.max(rev(diff(hierarchy$height))))
    X.comp <- cbind(computeCompressedDataFrame(X, grouping, compression))
    out$groups <- split(1:ncol(X), grouping)
  }
  if (selection=="pval") {
    criterion <- out$pval.count
  }

  res <- new("sicomore-model",
             groups       = out$groups      ,
             X.comp       = X.comp          ,
             compression  = compression     ,
             selection    = selection       )
  if(!stab){
    res$coefficients = out$coefficients
    res$cv.error = out$cv
  }

  return(res)
}

## __________________________________________________________________
##
## CROSS-VALIDATED COMPRESSED LASSO ALONG THE HIERARCHY
##
## explore the hierarchy to find the best level of compression (cut.levels are the best)
## The bottom of the hiearchy is discarded.... (too many variables :( )
getHierLevel.rhosicomore <- function(X, y, hc.object, compression, cut.levels, choice, depth.cut, mc.cores, stab, stab.param) {

  weights <- sqrt(1/abs(diff(hc.object$height)))
  cut.levels <- order(rev(c(max(weights),weights)))+1
  cut.levels <- cut.levels[cumsum(cut.levels) <= depth.cut*ncol(X)]
  weights <- c(0,rev(weights),0) # Reorder the weight to have a correspondance with cut.levels

  ## Build a data frame with all compressed variables from interesting cut levels
  hierarchy <- lapply(apply(cutree(hc.object, k = cut.levels), 2, list), unlist, recursive=FALSE)

  # Xcomp[,j] is a compressed variables and Xcomp.variables[[j]] is the corresponding vector of variables
  Xcomp.variables <- unlist(lapply(hierarchy,function(group) {lapply(1:max(group),function(k){which(group==k)})}), recursive=FALSE)
  uniqueIndex <- !duplicated(Xcomp.variables)
  Xcomp.variables <- Xcomp.variables[uniqueIndex] # Eliminate the group which are present at different level of the hierarchy
  Xcomp <- do.call(cbind,lapply(Xcomp.variables,function(variables) {computeCompressedDataFrameFromVariables(X, variables, compression)}))

  ## Adjusting the model
  if (!is.null(mc.cores)) doMC::registerDoMC(cores=mc.cores)
  else mc.cores = 1
  penalty.factor <- rep(weights[cut.levels],cut.levels)[uniqueIndex]

  if (stab == TRUE){
    stab.lasso <- stabs::stabsel(x = Xcomp, y = y, B = stab.param$B,
                          cutoff = stab.param$cutoff, PFER = stab.param$PFER,
                          fitfun = glmnet.lasso,
                          args.fitfun = list(penalty.factor = penalty.factor),
                          mc.cores = mc.cores, sampling.type = "MB")
    selected.groups <- as.numeric(stab.lasso$selected)
    groups <- Xcomp.variables[selected.groups]   # elements of each group of the best model
    res <- list(groups = groups)
  } else {
    cv <- glmnet::cv.glmnet(Xcomp, y, penalty.factor = penalty.factor, parallel=!is.null(mc.cores))
    selected.groups <- predict(cv, s=choice, type="nonzero")[[1]]
    groups <- Xcomp.variables[selected.groups]   # elements of each group of the best model
    coefficients <- predict(cv, s=choice, type="coefficients")[-1] # without intercept
    coefficients <- coefficients[selected.groups]
    res <- list(cv = data.frame(mean = cv$cvm, sd = cv$cvsd, lambda  = cv$lambda),
                groups = groups, coefficients = coefficients)

  }

  return(res)
}

## __________________________________________________________________
##
## CROSS-VALIDATED COMPRESSED LASSO ALONG THE HIERARCHY
##
## explore the hierarchy to find the best level of compression
getHierLevel.sicomore <- function(X, y, hc.object, compression, cut.levels, choice, depth.cut, mc.cores, stab, stab.param) {

  ## recovering the hierarchies at each cuting level.
  hierarchy <- lapply(apply(cutree(hc.object, k = cut.levels), 2, list), unlist, recursive=FALSE)

  if (!is.null(mc.cores)) doMC::registerDoMC(cores=mc.cores)
  else mc.cores = 1

  Xcomp.variables <- unlist(lapply(hierarchy,function(group) {lapply(1:max(group),function(k){which(group==k)})}), recursive=FALSE)
  uniqueIndex <- !duplicated(Xcomp.variables)
  Xcomp.variables <- Xcomp.variables[uniqueIndex] # Eliminate the group which are present at different level of the hierarchy
  Xcomp <- do.call(cbind,lapply(Xcomp.variables, function(variables) {computeCompressedDataFrameFromVariables(X, variables, compression)}))

  if (stab == TRUE){
       stab.lasso <- stabs::stabsel(x = Xcomp, y = y, B = stab.param$B,
                          cutoff = stab.param$cutoff, PFER = stab.param$PFER,
                          fitfun = glmnet.lasso, mc.cores = mc.cores, sampling.type = "MB")
    selected.groups <- as.numeric(stab.lasso$selected)
    groups <- Xcomp.variables[selected.groups]   # elements of each group of the best model
    res <- list(groups = groups)

  } else {
    cv <- glmnet::cv.glmnet(Xcomp, y, parallel=!is.null(mc.cores))
    selected.groups <- predict(cv, s=choice, type="nonzero")[[1]]
    groups <- Xcomp.variables[selected.groups]   # elements of each group of the best model
    coefficients <- predict(cv, s=choice, type="coefficients")[-1] # without intercept
    coefficients <- coefficients[selected.groups]
    res <- list(cv = data.frame(mean = cv$cvm, sd = cv$cvsd, lambda  = cv$lambda),
                groups = groups, coefficients = coefficients)

    # all.cv <- lapply(hierarchy, function(group) {
    #   return(glmnet::cv.glmnet(computeCompressedDataFrame(X, group, compression), y, parallel=!is.null(mc.cores)))
    # })
    # ## piles up the cv.error of each level of the hierarchy for all lambda
    # cv.error <- do.call(rbind,
    #                     lapply(all.cv, function(cv) {
    #                       data.frame(mean = cv$cvm, sd = cv$cvsd, lambda  = cv$lambda, ngroup = cv$glmnet.fit$dim[1])
    #                     }))
    #
    # ## Extract the best model
    # ibest <- match(cv.error$ngroup[which.min(cv.error$mean)], cut.levels)
    # selected.group <- predict(all.cv[[ibest]], s=choice, type="nonzero")[[1]]
    # groups    <- lapply(selected.group, function(x) which(hierarchy[[ibest]] == x))
    # coefficients <- predict(all.cv[[ibest]], s=choice, type="coefficients")[-1,] # without intercept
    # coefficients <- coefficients[selected.group]
    #
    # res <- list(cv = data.frame(mean = all.cv[[ibest]]$cvm, sd = all.cv[[ibest]]$cvsd, lambda  = all.cv[[ibest]]$lambda),
    #             groups = groups, coefficients = coefficients)
  }

  return(res)
}

## __________________________________________________________________
##
## CROSS-VALIDATED MULTILAYER GROUP-LASSO ALONG THE HIERARCHY
##
getHierLevel.MLGL <- function(X, y, hc.object, compression, cut.levels, choice, depth.cut, mc.cores, stab, stab.param) {

  ## Grimponprez's weights
  #weights <- c(0, sqrt(1/abs(diff(hc.object$height))), 0)
  #weights[-(ncol(X)+1  - cut.levels)] <- 0

  ## Adjusting the model + cross-validation
  #fit <- MLGL::MLGL(X, y, hc=hc.object, weightLevel = weights)
  fit <- MLGL::MLGL(X, y)
  #cv.error <- MLGL::cv.MLGL(X, y, hc=hc.object, weightLevel = weights)
  cv.error <- MLGL::cv.MLGL(X, y)

  ## Extracting the best model
  ibest <- switch(choice,
                  "lambda.min" = match(cv.error$lambda.min, fit$lambda),
                  "lambda.1se" = match(cv.error$lambda.1se, fit$lambda))
  variables <- unname(unlist(fit$var[ibest]))
  double <- which(duplicated(variables))
  if (length(double) != 0){
    variables <- variables[-double]
    fit.beta <- fit$beta[[ibest]][-double]
    fit.group <- fit$group[[ibest]][-double]
  } else {
    fit.beta <- fit$beta[[ibest]]
    fit.group <- fit$group[[ibest]]
  }
  o <- order(variables)
  groups       <- split(variables, fit.group)
  coefficients <- as.numeric(sparseVector(fit.beta,unique(variables[o]),ncol(X)))

  return(list(cv = data.frame(mean=cv.error$cvm, sd=cv.error$cvsd, lambda=cv.error$lambda),
              groups = groups, coefficients = coefficients))
}



# getHierLevel.PVAL <- function(X, y, hc.object, compression, cut.levels, choice, threshold=1, depth.cut, mc.cores) {
#
#   ## recovering the hierarchies at each cuting level.
#   hierarchy <- lapply(apply(cutree(hc.object, k = cut.levels), 2, list), unlist, recursive=FALSE)
#
#   ## recover the grid of penalties used for comparison
#   all.pval <- lapply(hierarchy, function(group) {
#     return(computePval(computeCompressedDataFrame(X, group, compression), y, threshold=threshold)) # return pvalue for each supervariable
#                                                                                # and number of significative pvalue
#                                                                                # after multiple pvalue correction
#
#   })
#   browser()
#   ## piles up the pvalue number  of each level of the hierarchy for all cut levels
#   pval.count  <- do.call(rbind,
#                       lapply(all.pval, function(pval) {
#                         data.frame( count = pval$count , nbgroup = length(pval$pval.raw))
#                       }))
#
#   ## Extract the best model
#   ibest <- match(pval.count$nbgroup[which.max(pval.count$count)], cut.levels)
#   selected.group <- which(all.pval[[ibest]]$pval.adjusted > threshold)    #
#   groups    <- lapply(selected.group, function(x) which(hierarchy[[ibest]] == x))
#
#   return(list(pval.count = pval.count,  groups = groups, coefficients = all.pval[[ibest]]$pval.raw))
# }
