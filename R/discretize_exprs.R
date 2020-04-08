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
#' "min_max_percent" (Similat to the "varying width", a binarization threshold in a percent of the min-max range is set. (minmaxpercent param)),
#' "GMM" (A Gaussian Mixture Model as implemented by the package mclust, trying to fit 2:5 Gaussians)
#' @param number_of_bins Number of equal-width bins for discretization.
#' Note: it is a binary discretization, with the
#' first bin becoming one class ('low') and the other bins, another class ('high').
#' Defaults to 3.
#' @param alpha Modulator for the "mean_sd" method.Enlarges (>1) or shrinks (<1) the "medium" interval. Defaults to 1.
#' @param centers Modulator for the "kmeans" method. Defaults to 3.
#' @param min_max_cutoff <- Modulator for the "min_max_percent" method. Defaults to 0.25.
#' @param progress_bar Enables a progress bar for the discretization. Defaults to TRUE.
#' @return A data frame with the discretized features in the same order as previously
#' @import mclust
#' @import pbapply
#' @export
#' @examples
#' data(scDengue)
#' exprs <- SummarizedExperiment::assay(scDengue, 'logcounts')
#' discrete_expression <- as.data.frame(discretize_exprs(exprs))
#' head(discrete_expression[, 1:4])

discretize_exprs <-
  function(expression_table,
           number_of_bins = 3,
           method = "varying_width",
           alpha = 1,
           centers = 3,
           min_max_cutoff = 0.25,
           progress_bar = TRUE) {
    if (method == "varying_width") {
      .split_fun <- .split_vector_in_two_varying_width
      formals(.split_fun)[2] <- number_of_bins
    } else   if (method == "mean") {
      .split_fun <- .split_vector_in_two_by_mean
    }  else   if (method == "median") {
      .split_fun <- .split_vector_in_two_by_median
    } else   if (method == "mean_sd") {
      .split_fun <- .split_vector_in_three_by_mean_sd
      formals(.split_fun)[2] <- alpha
    } else   if (method == "kmeans") {
      .split_fun <- .split_vector_by_kmeans
      formals(.split_fun)[2] <- centers
    }  else   if (method == "min_max_percent") {
      .split_fun <- .split_vector_in_two_by_min_max_thresh
      formals(.split_fun)[2] <- min_max_cutoff
    } else   if (method == "GMM") {
      .split_fun <- function(X) {
        sink('aux')
        fit <- Mclust(X, G = 2:5, warn = FALSE)
        sink(NULL)
        return(fit$classification)
      }
    }

    if (progress_bar){
      discrete_expression <-
        pbapply(expression_table, 1, .split_fun)
    } else {
      discrete_expression <-
        apply(expression_table, 1, .split_fun)
    }

    discrete_expression <- as.data.frame(t(discrete_expression), stringsAsFactors=TRUE)
    rownames(discrete_expression) <- rownames(expression_table)
    return(discrete_expression)
  }
