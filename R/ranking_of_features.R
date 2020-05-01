#'  Symmetrical Uncertainty diagnostic
#'
#' This functions runs symmetrical uncertainty for a feature table and a class, returning
#' the scores of symmetrical uncertainty for all features
#'
#' @param feature_table A table of features (observations in rows, variables in columns)
#' @param target_vector A target vector, factor containing classes of the observations. Note: the
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
get_su_for_feature_table_and_vector <- function(feature_table, target_vector, samples_in_rows = FALSE, bar_of_progress = FALSE) {
  if (!samples_in_rows){
   feature_table <- t(feature_table)
  }
  if (!is.data.frame(feature_table)){
   feature_table <- data.frame(feature_table)
  }
  if (bar_of_progress){
    su_values_for_features_with_regards_to_class <- pbapply(feature_table, 2, function(feature_vector, target_vector) {
      get_SU_for_vector_pair(feature_vector, target_vector)
    }, target_vector)

  } else{
  su_values_for_features_with_regards_to_class <- apply(feature_table, 2, function(feature_vector, target_vector) {
    get_SU_for_vector_pair(feature_vector, target_vector)
  }, target_vector)
  
  }
  su_values_for_features_with_regards_to_class <- as.data.frame(sort(su_values_for_features_with_regards_to_class,decreasing = TRUE))
  su_values_for_features_with_regards_to_class$gene <- rownames(su_values_for_features_with_regards_to_class)
  su_values_for_features_with_regards_to_class
}

#' @importMethodsFrom SummarizedExperiment assay colData
NULL
#'  Get information gain
#'
#' This functions runs information gain for a feature table and a class, returning
#' the scores of information gain for all features
#' @param feature_table A table of features (observations in rows, variables in columns)
#' @param target_vector A target vector, factor containing classes of the observations. Note: the
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
#' ig_values <- get_ig(discrete_expression[,],target[])
#' ig_values[1:10,]



get_ig_for_feature_table_and_vector <- function(feature_table, target_vector) {
 feature_table <- t(feature_table)
 feature_table <- data.frame(feature_table)
  ig_values_for_features_with_regards_to_class <- apply(feature_table, 2, function(feature_vector, target_vector) {
    get_IG_for_vector_pair(feature_vector, target_vector)
  }, target_vector)
  as.data.frame(sort(ig_values_for_features_with_regards_to_class, decreasing = TRUE))
}


####### Deprecated functions ########

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
#' @name get_su-deprecated
#' @seealso \code{\link{FCBF-deprecated}}
#' @keywords internal
NULL

#' @rdname FCBF-deprecated
#' @section \code{get_su}:
#' For \code{get_su}, use \code{\link{get_su_for_feature_table_and_vector}}.
#'
#' @export
get_su <- function(x, y) {
  .Deprecated("get_su_for_feature_table_and_vector")
  "get_su function is no longer available"
}


#'  Get information gain
#'
#' This functions runs information gain for a feature table and a class, returning
#' the scores of information gain for all features
#' @param x A table of features (observations in rows, variables in columns)
#' @param y A target vector, factor containing classes of the observations. Note: the
#' observations must be in the same order as the parameter x.
#' @return A dataframe containing the SU values for each feature
#' @rdname FCBF-deprecated
#' @section \code{get_su}:
#' For \code{get_su}, use \code{\link{get_ig_for_feature_table_and_vector}}.
#'
#' @export
get_ig <- function(x, y) {
  .Deprecated("get_ig_for_feature_table_and_vector")
  "get_ig function is no longer available"
}



