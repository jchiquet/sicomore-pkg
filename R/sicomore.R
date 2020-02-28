sicomore <- function(y,
                     X.list,
                     compressions = rep("mean", length(X.list)),
                     selection =  c("rho-sicomore", "sicomore", "mlgl"),
                     cuts = lapply(X.list, function(X) unique(round(10^seq(log10(ncol(X)), log10(2), len=100)))),
                     choice=c("lambda.min", "lambda.min"),
                     method.clus = c("ward.D2","ward.D2"),
                     depth.cut = c(3,3),
                     mc.cores = NULL,
                     taxonomy = NULL,
                     verbose = TRUE,
                     stab = TRUE,
                     stab.param  = list(B = c(100,100), cutoff = c(.75,.75), PFER = c(1,1)))
{

  selection <- match.arg(selection, c("rho-sicomore", "sicomore", "mlgl"))
  ndata <- length(X.list)
  models <- vector("list", ndata)
  hierarchies  <- vector("list" , length(X.list))
  for (i in seq_along(X.list)) {
    if (verbose == TRUE) cat("\nConsidering hierarchy",i, " - ", selection, "selection method.")

    if (method.clus[i] == "ward.D2"){
      hierarchies[[i]] <- hclust(dist(t(scale(X.list[[i]]))), method="ward.D2")
      models[[i]] <- getHierLevel(X = X.list[[i]], y, hierarchies[[i]], cut.levels = cuts[[i]], compression=compressions[i],
                                  selection=selection, choice=choice[i], depth.cut = depth.cut[i], mc.cores=mc.cores,
                                  stab, stab.param = lapply(stab.param, function(x) x[[i]]))
    }
    if (method.clus[i] == "snpClust") {
      if (ncol(X.list[[i]]) > 600) h <- 600
      else h <- ncol(X.list[[i]]) - 1
      #hierarchies[[i]] <- adjclust::snpClust(X.list[[i]], h=h)
      hierarchies[[i]] <- BALD::cWard(X.list[[i]], h=h, heaps = TRUE)
      models[[i]] <- getHierLevel(X.list[[i]], y, hierarchies[[i]], cut.levels = cuts[[i]], compression=compressions[i],
                                  selection=selection, choice=choice[i], depth.cut = depth.cut[i], mc.cores=mc.cores,
                                  stab, stab.param = lapply(stab.param, function(x) x[[i]]))
    }
    if (method.clus[i] == "noclust"){
      if (!is.null(mc.cores)) doMC::registerDoMC(cores=mc.cores)
      else mc.cores = 1

      if (stab == TRUE){
        stab.lasso <- stabs::stabsel(x = X.list[[i]], y = y, B = 50,
                              fitfun = glmnet.lasso, cutoff = 0.75,
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

      if (!is.null(mc.cores)) doMC::registerDoMC(cores=mc.cores)
      else mc.cores = 1
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
    return(pval)
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

  return(new("sicomore-fit", pval = pval.inter, pval.beta1 = pval.beta1, pval.beta2 = t(pval.beta2),
             tuplets = tuplets, models=models, dim=sapply(X.list, ncol)))
}

