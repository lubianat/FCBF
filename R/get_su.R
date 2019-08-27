#'  Symmetrical Uncertainty diagnostic
#'
#' This functions runs symmetrical uncertainty for a feature table and a class, returning
#' the scores of symmetrical uncertainty for all features
#'
#' @param x A table of features (observations in rows, variables in columns)
#' @param y A target vector, factor containing classes of the observations. Note: the
#' observations must be in the same order as the parameter x.
#' @param samples_in_rows A flag for the case in which samples are in rows and variables/genes in columns. Defaults to FALSE.
#' @param bar_of_progress A flag to show progress. Defaults to FALSE.
#' @return A dataframe containing the SU values for each feature
#' @import pbapply
#' @export
#' @examples
#' data(scDengue)
#' exprs <- SummarizedExperiment::assay(scDengue, 'logcounts')
#' discrete_expression <- as.data.frame(discretize_exprs(exprs))
#' infection <- SummarizedExperiment::colData(scDengue)
#' target <- infection$infection
#' su_values <- get_su(discrete_expression[,],target[])
#' su_values[1:10,]
get_su <- function(x, y, samples_in_rows = FALSE, bar_of_progress = FALSE) {
  if (!samples_in_rows){
    x <- t(x)
  }
  if (!is.data.frame(x)){
    x <- data.frame(x)
  }
  if (bar_of_progress){
    su_ic <- pbapply(x, 2, function(xx, yy) {
      SU(xx, yy)
    }, y)

  } else{
  su_ic <- apply(x, 2, function(xx, yy) {
    SU(xx, yy)
  }, y)
  }
  su_ic <- as.data.frame(sort(su_ic,decreasing = TRUE))
  su_ic$gene <- rownames(su_ic)
  su_ic
  }

