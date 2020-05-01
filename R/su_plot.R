#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importMethodsFrom SummarizedExperiment assay
NULL

#'  Symmetrical Uncertainty diagnostic
#'
#' This functions runs symmetrical uncertainty for a feature table and a class, returning
#' an histogram of the scores
#' @param feature_table A table of features (observations in rows, variables in columns)
#' @param target_vector A target vector, factor containing classes of the observations. Note: the
#' observations must be in the same order as the parameter x.
#' @return Plots an histogram of symmetrical uncertainty values regarding the class.
#' @export
#' @examples
#' data(scDengue)
#' exprs <- SummarizedExperiment::assay(scDengue, 'logcounts')
#' discrete_expression <- as.data.frame(discretize_exprs(exprs))
#' infection <- SummarizedExperiment::colData(scDengue)
#' target <- infection$infection
#' su_plot(discrete_expression,target)



su_plot <- function(feature_table, target_vector) {
  feature_table <- t(feature_table)
  feature_table <- data.frame(feature_table)
  su_values_for_features_with_regards_to_class <- apply(feature_table, 2, function(feature_vector, target_vector) {
    get_SU_for_vector_pair(feature_vector, target_vector)
  }, target_vector)
  su_values_for_features_with_regards_to_class <- as.data.frame(su_values_for_features_with_regards_to_class)
  p1 <- ggplot(su_values_for_features_with_regards_to_class, (aes(x = su_values_for_features_with_regards_to_class))) +
    geom_histogram(binwidth = 0.008) +
    xlab('Correlation by Symmetrical Uncertainty to target') +
    ylab('Variable counts') +
    ggtitle('Histogram of SU values of each variable to the target classes')

  p2 <- ggplot(su_values_for_features_with_regards_to_class, (aes(y = su_values_for_features_with_regards_to_class))) +
    geom_boxplot(outlier.alpha = 0.5) +
    ylab('Correlation by Symmetrical Uncertainty to target') +
    ggtitle('Boxplot of SU values of each variable to the target classes') +
    theme(axis.text.x = element_blank())

  grid.arrange(p1, p2, nrow = 1)
}
