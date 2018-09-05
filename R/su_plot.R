

#'  Symmetrical Uncertainty diagnostic
#'
#' This functions runs symmetrical uncertainty for a feature table and a class, returning
#' an histogram of the scores
#' @param x A table of features (observations in rows, variables in columns)
#' @param y A target vector, factor containing classes of the observations. Note: the
#' observations must be in the same order as the parameter x.
#' @return Plots an histogram of symmetrical uncertainty values regarding the class.
#' @export
#' @importFrom graphics hist
#' @examples
#' data(single_cell_dengue_exprs)
#' discrete_expression <- as.data.frame(discretize_exprs(single_cell_dengue_exprs))
#' head(discrete_expression[,1:4])
#' data("single_cell_dengue_annot")
#' su_plot(discrete_expression,target)


su_plot <- function(x, y) {
  x <- data.frame(x)
  su_ic <- apply(x, 2, function(xx, yy) {
    SU(xx, yy)
  }, y)
  hist(su_ic)
}
