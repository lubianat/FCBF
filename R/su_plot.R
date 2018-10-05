

#'  Symmetrical Uncertainty diagnostic
#'
#' This functions runs symmetrical uncertainty for a feature table and a class, returning
#' an histogram of the scores
#' @param x A table of features (observations in rows, variables in columns)
#' @param y A target vector, factor containing classes of the observations. Note: the
#' observations must be in the same order as the parameter x.
#' @return Plots an histogram of symmetrical uncertainty values regarding the class.
#' @export
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @examples
#' data(scDengue)
#' exprs <- assay(scDengue, 'logcounts')
#' discrete_expression <- as.data.frame(discretize_exprs(exprs))
#' su_plot(discrete_expression,target)


su_plot <- function(x, y) {
  require(ggplot2)
  require(gridExtra)
  x <- t(x)
  x <- data.frame(x)
  su_ic <- apply(x, 2, function(xx, yy) {
    SU(xx, yy)
  }, y)
  su_ic <- as.data.frame(su_ic)
  p1 <- ggplot(su_ic, (aes(x=su_ic)))+
    geom_histogram(binwidth= 0.008) +
    xlab('Correlation by Symmetrical Uncertainty to target') +
    ylab('Variable counts') +
    ggtitle('Histogram of SU values of each variable to the target classes')

  p2 <- ggplot(su_ic, (aes(y=su_ic)))+
    geom_boxplot(outlier.alpha = 0.5) +
    ylab('Correlation by Symmetrical Uncertainty to target') +
    ggtitle('Boxplot of SU values of each variable to the target classes') +
    theme(axis.text.x=element_blank())

  grid.arrange(p1, p2, nrow = 1)
}
