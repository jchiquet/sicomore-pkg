#' getHierLevel
#'
#' A function to select groups of variables which are good at predicting a given phenotype.
#' The groups considered corresponds to various cut levels in a user defined hierarchy.
#' The selection is performed by various penalty-based regression methods (weighted-LASSO or group-LASSO).
#'
#' @param X input matrix
#' @param y response variable
#' @param hc.object output of a hierarchical clustering algorithm in the \code{hclust} format (must be an "hclust" object)
#' @param selection method used to perform variable selection. Either 'sicomore', 'rho-sicomore'  or 'mlgl' (see details). Default is 'rho-sicomore'.
#' @param compression a string (either "mean" or "SNP.dist"). Indicates how groups of variables are compressed before variable selection is performed at each level of the hierarchy. Only relevant for 'sicomore' or 'rho-sicomore'.
#' @param choice a string (either "lambda.min" or "lambda.1se"). Indicates how the tuning parameter is chosen in the penalized regression approach
#' @param depth.cut an integer specifying the depth of the search space for the variable selection part of the algorithm.
#' This argument allows to increase the speed of the algorithm by restraining the search space without affecting too much the performance.
#' A value between 3 and 6 is recommended, the smaller the faster.
#' @param stab A boolean indicating if the algorithm perform a lasso stability selection using stabsel function from stabs package.
#' @param stab.param A list of parameter for the stabsel function if stab = TRUE.
#' The parameters to choose are the FWER (1 by default), cut-off (0.75 by default) and bootstrap number (200 by default).
#' @param grp.min Minimum number of groups to consider for the highest level in the hierarchy.
#' Correspond to the highest allowed cut in the hierarchy. If NA, no restriction is given.
#' @param mc.cores an integer for the number of cores to use in the parallelization of the cross-validation and some other functions. Default is 1.
#'
#' @details The methods for variable selection are variants of the LASSO or the group-LASSO designed to perform selection of interaction between multiple hierarchies:
#' 'sicomore' and 'rho-sicomore' (see \insertCite{sicomore;textual}{sicomore}, \insertCite{park;textual}{sicomore}) use a LASSO penalty on compressed groups of variables along the hierarchies to select the
#' interactions. The rho-sicomore variant is a weighted version of sicomore, which weights depend on the levels in the hierarchies. The method 'mlgl' of \insertCite{grimonprez_PhD;textual}{sicomore}
#' uses a weigthed group-Lasso penalty which does not require compression but is more computationally demanding.
#'
#' @return an RC object with class 'sicomore-model', with methods \code{nGrp()}, \code{nVar()}, \code{getGrp()}, \code{getVar()}, \code{getCV()}, \code{getX.comp()}, \code{getCoef()} and with the following fields:
#' \itemize{
#'  \item{groups:}{a list with the selected groups of predictors}
#'  \item{coefficients:}{a vector with the estimated coefficients (one per selected group) if stab=FALSE}
#'  \item{X.comp:}{The compressed version of the original input matrix (as many columns as number of selected groups)}
#'  \item{cv.error:}{for the best grouping, a data frame showing the cross-validation error used in the variable selection procedure if stab=FALSE}
#'  \item{selection:}{the selection method used}
#'  \item{compression:}{the compression method used}
#'  \item{group_inference:}{the group selection infered by "lasso" or "hclust" if no selection by lasso.}
#' }
#' @include utils.R
#' @import glmnet MLGL
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{}
getHierLevel <- function(X,
                         y,
                         hc.object,
                         selection = c("rho-sicomore", "sicomore", "mlgl"),
                         compression = "mean",
                         depth.cut = 3,
                         grp.min = NA,
                         choice=c("lambda.min", "lambda.1se"),
                         mc.cores=NULL,
                         stab = FALSE,
                         stab.param = list(B = c(100,100), cutoff = c(.75,.75), PFER = c(1,1))) {

  ## _________________________________________________________________
  ##
  ## PRETREATMENT
  ##
  choice <- match.arg(choice)

  ## forcing the matrix type for glmnet
  if (!is.matrix(X)) X <- as.matrix(X)

  if (class(hc.object) == "chac") class(hc.object) <- "hclust"
  stopifnot(class(hc.object) == "hclust")

  ## _________________________________________________________________
  ##
  ## GET THE BEST COMPRESSION LEVEL IN THE HIERARCHY
  ##
  weights <- switch(selection,
                     "rho-sicomore" = sqrt(1/abs(diff(hc.object$height))),
                     "sicomore"     = rep(1, length(hc.object$height)),
                     "mlgl"         = NULL)
  if (selection == 'mlgl')
    out <- .MLGL(X, y, hc.object, choice)
  else
    out <- .sicomore(X, y, weights, hc.object, compression, choice, depth.cut,
                     mc.cores, stab, stab.param, grp.min)

  ## _________________________________________________________________
  ##
  ## POST-TREATMENTS
  ##

  res <- new("sicomore-model",
             compression  = compression     ,
             selection    = selection       )

  if (length(out$groups) > 0) {
    grouping   <- rep(1:length(out$groups), sapply(out$groups, length))
    X.comp     <- cbind(computeCompressedDataFrame(X[, unlist(out$groups)], grouping, compression))
    out$groups <- setNames(out$groups, NULL)
    if(!stab) out$coefficients <- setNames(out$coefficients, paste("group", 1:length(out$groups)))
    res$group_inference <- "lasso"
  } else {
    grouping   <- cutree(hc.object, k=1+which.max(rev(diff(hc.object$height))))
    X.comp     <- cbind(computeCompressedDataFrame(X, grouping, compression))
    out$groups <- split(1:ncol(X), grouping)
    res$group_inference <- "hclust"
  }
  res$groups <- out$groups
  res$X.comp <- X.comp
  if(!stab){
    res$coefficients <- out$coefficients
    res$cv.error     <- out$cv
  }

  res
}

## __________________________________________________________________
##
## CROSS-VALIDATED COMPRESSED LASSO ALONG THE HIERARCHY
##
## explore the hierarchy to find the best level of compression (cut.levels are the best)
## The bottom of the hiearchy is discarded.... (too many variables :( )
.sicomore <- function(X, y, weights, hc.object, compression, choice, depth.cut,
                      mc.cores, stab, stab.param, grp.min) {

  cut.levels <- order(rev(c(max(weights),weights))) + 1
  cut.levels <- cut.levels[cumsum(cut.levels) <= depth.cut*ncol(X)]
  n.grp.levels <- sapply(cut.levels, function(k) length(unique(cutree(hc.object, k = k))))
  if (!is.na(grp.min) & any(n.grp.levels <= grp.min)) cut.levels[-which(n.grp.levels <= grp.min)] ## Remove out high levels in the hierarchy
  weights <- c(0,rev(weights),0) # Reorder the weight to have a correspondance with cut.levels

  ## Build a data frame with all compressed variables from interesting cut levels
  hierarchy <- lapply(apply(cutree(hc.object, k = cut.levels), 2, list), unlist, recursive=FALSE)

  # Xcomp[,j] is a compressed variables and Xcomp.variables[[j]] is the corresponding vector of variables
  Xcomp.variables <- unlist(lapply(hierarchy,function(group) {
    lapply(1:max(group), function(k){which(group==k)})
    }), recursive=FALSE)
  uniqueIndex <- !duplicated(Xcomp.variables)
  Xcomp.variables <- Xcomp.variables[uniqueIndex] # Eliminate the group which are present at different level of the hierarchy
  Xcomp <- do.call(cbind,lapply(Xcomp.variables, function(variables) {
    computeCompressedDataFrameFromVariables(X, variables, compression)
    }))

  ## Adjusting the model
  penalty.factor <- rep(weights[cut.levels],cut.levels)[uniqueIndex]

  if (stab == TRUE){
    stab.lasso <- stabs::stabsel(x = Xcomp, y = y, B = stab.param$B,
                          cutoff = stab.param$cutoff, PFER = stab.param$PFER,
                          fitfun = stabs::glmnet.lasso,
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

  res
}

## __________________________________________________________________
##
## CROSS-VALIDATED MULTILAYER GROUP-LASSO ALONG THE HIERARCHY
##

.MLGL <- function(X, y, hc.object, choice, mc.cores) {

  ## Adjusting the model + cross-validation
  fit      <- MLGL::MLGL(X, y)
  cv.error <- MLGL::cv.MLGL(X, y)

  ## Extracting the best model
  ibest <- switch(choice,
                  "lambda.min" = match(cv.error$lambda.min, fit$lambda),
                  "lambda.1se" = match(cv.error$lambda.1se, fit$lambda))
  variables <- unname(unlist(fit$var[ibest]))
  double    <- which(duplicated(variables))
  if (length(double) != 0){
    variables <- variables[-double]
    fit.beta  <- fit$beta[[ibest]][-double]
    fit.group <- fit$group[[ibest]][-double]
  } else {
    fit.beta  <- fit$beta[[ibest]]
    fit.group <- fit$group[[ibest]]
  }
  o <- order(variables)
  groups       <- split(variables, fit.group)
  coefficients <- as.numeric(Matrix::sparseVector(fit.beta,unique(variables[o]),ncol(X)))

  list(cv = data.frame(mean=cv.error$cvm, sd=cv.error$cvsd, lambda=cv.error$lambda),
              groups = groups, coefficients = coefficients)
}

