#' @importMethodsFrom SummarizedExperiment assay
NULL


#' discretize_exprs
#' Simple discretizing of gene expression
#'
#' This function takes the range of values for each gene in a
#' previously normalized expression table (genes/variables in rows, samples/
#' observations in columns) and uses it for a width-based discretization.
#' Each feature is divide into "n" bins of equal width. The first bin is
#' attributed the class 'low' and the next bins are assigned to "high".
#' It transposes the original expression table.
#'
#' @param expression_table  A previously normalized expression table
#' Note: this might drastically change the number of selected features.
#' @param method Method applied to all genes for discretization. Methods available: "varying_width"
#'  (Varying width binarization, default, described in function description. Modulated by the number_of_bins param),
#' "mean" (Split in ON/OFF by each gene mean expression),
#' "median" (Split in ON/OFF by each gene median expression),
#' "mean_sd"(Split in low/medium/high by each assigning "medium" to the interval between mean +- standard_deviation.
#' Modulated by the alpha param, which enlarges (>1) or shrinks (<1) the "medium" interval. ),
#' ),
#' "kmeans"(Split in different groups by the kmeans algorithm. As many groups as specified by the centers param) and
#' "min_max_%" (Similat to the "varying width", a binarization threshold in a % of the min-max range is set. (minmax% param))
#' @param number_of_bins Number of equal-width bins for discretization.
#' Note: it is a binary discretization, with the
#' first bin becoming one class ('low') and the other bins, another class ('high').
#' Defaults to 3.
#' @param alpha Modulator for the "mean_sd" method.Enlarges (>1) or shrinks (<1) the "medium" interval. Defaults to 1.
#' @param centers Modulator for the "kmeans" method. Defaults to 3.
#' @param min_max_cutoff <- Modulator for the "min_max_%" method. Defaults to 0.25.
#' @return A data frame with the discretized features in the same order as previously
#' @export
#' @examples
#' data(scDengue)
#' exprs <- SummarizedExperiment::assay(scDengue, 'logcounts')
#' discrete_expression <- as.data.frame(discretize_exprs(exprs))
#' head(discrete_expression[,1:4])


discretize_exprs <- function(expression_table, number_of_bins = 3, method = "varying_width", alpha = 1, centers = 3, min_max_cutoff = 0.25) {

  if (method == "varying_width"){
    .split_fun <-.split_vector_in_two_varying_width
  } else   if (method == "mean"){
    .split_fun <-.split_vector_in_two_by_mean
  }  else   if (method == "median"){
    .split_fun <-.split_vector_in_two_by_median
  } else   if (method == "mean_sd"){
    .split_fun <-.split_vector_in_three_by_mean_sd
    formals(.split_fun) <- alpha
  } else   if (method == "kmeans"){
    .split_fun <-.split_vector_by_kmeans
    formals(.split_fun) <- centers
  }  else   if (method == "min_max_%"){
    .split_fun <-.split_vector_in_two_by_min_max_thresh
    formals(.split_fun) <- centers
  }


  discrete_expression <-
    apply(expression_table, 1, .split_fun, n_of_bins = number_of_bins)
  discrete_expression <- as.data.frame(t(discrete_expression))
  rownames(discrete_expression) <- rownames(expression_table)
  return(discrete_expression)
}
