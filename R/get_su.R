#'  Symmetrical Uncertainty diagnostic
#'
#' This functions runs symmetrical uncertainty for a feature table and a class, returning
#' the scores of symmetrical uncertatinty for all features
#' @param x A table of features (observations in rows, variables in columns)
#' @param y A target vector, factor containing classes of the observations. Note: the
#' observations must be in the same order as the parameter x.
#' @return A dataframe containing the SU values for each feature
#' @export
#' @examples
#'
#' @examples
#' data(scDengue)
#' exprs <- SummarizedExperiment::assay(scDengue, 'logcounts')
#' discrete_expression <- as.data.frame(discretize_exprs(exprs))
#' infection <- SummarizedExperiment::colData(scDengue)
#' target <- infection$infection
#' su_values <- get_su(discrete_expression[,],target[])
#' su_values[1:10,]



get_su <- function(x, y) {
  x <- t(x)
  x <- data.frame(x)
  su_ic <- apply(x, 2, function(xx, yy) {
    SU(xx, yy)
  }, y)
  as.data.frame(sort(su_ic,decreasing = T))
}
