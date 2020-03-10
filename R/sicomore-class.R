#' Sicomore-model RC object
#'
#' @docType class
#' @import ggplot2
#' @export
sicomore.model <-
  setRefClass("sicomore-model",
    fields  = list(
          coefficients = "numeric",
          groups       = "list",
          cv.error     = "data.frame",
          X.comp       = "matrix",
          compression  = "character",
          selection    = "character"
        ),
              methods = list(
                nGrp      = function() return(length(groups)),
                nVar      = function() return(length(unique(unlist(groups)))),
                getGrp    = function() return(groups),
                getVar    = function() return(unique(unlist(groups))),
                getX.comp = function() return(X.comp),
                getCV     = function() return(cv.error),
                getCoef   = function() return(coefficients)
              )
  )

sicomore.model$methods(plotCV = function(plot=TRUE) {
  p <- ggplot(cv.error, aes(x=lambda, y=mean)) +
    geom_smooth(aes(ymin=mean-sd, ymax=mean+sd), stat="identity") + coord_trans(x="log")
  if (plot)
    print(p)
  p
})

#' Sicomore RC object
#'
#' @docType class
#' @export
sicomore.fit <-
  setRefClass("sicomore-fit",
              fields  = list(
                pval         = "array",
                pval.beta1  = "array",
                pval.beta2  = "array",
                models       = "list",
                tuplets      = "list",
                dim          = "numeric"
                ),
              methods = list(
                plot = function(threshold=0.05, main=paste0("Significant Interaction effects (thres.=",threshold,")")) {
                   if (length(models) > 2) {
                     stop("not implemented yet")
                   } else {
                     theta <- (res$getSignificance(effect="theta") > 1-threshold) + 0
                     ggplot(data =  reshape2::melt(theta), aes(Var1, Var2, fill=value)) +
                       geom_tile(show.legend=FALSE, colour = "grey80") +
                       scale_fill_gradient2(low = "grey60", high = "grey40") +
                       theme_minimal()+ # minimal theme
                       ggtitle(main) +
                       labs(x = 'Variables in X_1', y = 'Variables in X_2') +
                       coord_fixed()

                   }
                },
                getSignificance = function(effect=c("theta","beta1","beta2")) {
                  significance <- array(0, dim = dim)
                  for (grp.ind in tuplets) {
                    i.expand <- lapply(1:length(models), function(i) {
                      return(models[[i]]$getGrp()[[grp.ind[i]]])
                    } )
                    if (length(dim) == 2) {
                      if (effect == "theta") significance[i.expand[[1]], i.expand[[2]]] <- 1 - pval[rbind(grp.ind)]
                      if (effect == "beta1") significance[i.expand[[1]], i.expand[[2]]] <- 1 - pval.beta1[rbind(grp.ind)]
                      if (effect == "beta2") significance[i.expand[[1]], i.expand[[2]]] <- 1 - pval.beta2[rbind(grp.ind)]
                    } else {
                      stop("Need to think about this a bit more")
                    }
                  }
                  return(significance)
                }
              )
  )

# pval <- matrix(NA,ncol(res.hier[[1]]$X.comp),ncol(res.hier[[2]]$X.comp))
# significance <- Matrix(0,ncol(X.list[[1]]),ncol(X.list[[2]]))
#
