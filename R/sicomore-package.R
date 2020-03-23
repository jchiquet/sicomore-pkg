#' sicomore
#'
#' Selection of Interaction effects in COmpressed  Multiple Omics REpresentation
#'
#' @author Julien Chiquet \email{julien.chiquet@@inrae.fr}
#' @author Florent Guinot \email{fly.guinot@@gmail.com}
#' @author Marie Szafranski  \email{marie.szafranski@@math.cnrs.fr}
#' @author Christophe Ambroise \email{christophe.ambroise@@genopole.cnrs.fr}
#'
#' @docType package
#' @useDynLib sicomore, .registration=TRUE
#' @importFrom stats cutree dist hclust lm median p.adjust pf predict setNames
#' @importFrom methods new slot
#' @importFrom Matrix sparseVector
#' @importFrom stabs glmnet.lasso stabsel
#' @importFrom utils head tail
#' @importFrom snpStats ld
#' @importFrom matrixStats rowCumsums colCumsums
#' @name sicomore
NULL
