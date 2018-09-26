#' Simple discretizing of gene expression
#'
#' This function takes the range of values for each gene in a
#' previously normalized expression table (genes/variables in rows,
#' samples/observations in columns) and uses it for a width-based discretization.
#' Each feature is divide into "n" bins of equal width. The first bin is
#' attributed the class 'low' and the next bins are assigned to "high".
#' It transposes the original expression table.
#'
#' @param expression_table  A previously normalized expression table (genes/variables in rows,
#' samples/observations in columns)
#' Note: this might drastically change the number of selected features.
#' @param number_of_bins Number of equal-width bins for discretization.
#' Note: it is a binary discretization, with the
#' first bin becoming one class ('low') and the other bins, another class ('high').
#' Defaults to 3.
#' @export
#' @return A data frame with the discretized features (observations in rows, variables in columns).
#' Notice: it is transposed regarding the data frame inserted as input.
#' @examples
#' data(single_cell_dengue_exprs)
#' discrete_expression <- as.data.frame(discretize_exprs(single_cell_dengue_exprs))
#' head(discrete_expression[,1:4])

discretize_exprs <- function(expression_table, number_of_bins = 3) {
  apply(expression_table, 1, split_vector_in_two, n_of_bins = number_of_bins)
}


split_vector_in_two <-
  function(gene_expression_across_samples,
           n_of_bins = 3) {
    gene_expression_across_samples <-
      as.numeric(gene_expression_across_samples)
    max_expression = max(gene_expression_across_samples)
    min_expression = min(gene_expression_across_samples)
    break_size = (max_expression - min_expression) / n_of_bins
    return(ifelse(gene_expression_across_samples < (min_expression+break_size),  'low', 'high'))
  }

# exprs <- data.frame(fread('../RawCounts_Infection.txt'))
# exprs <- column_to_rownames(exprs, 'Probes')
# exprs <- exprs[1:100, ]
#
# discretize_exprs(exprs)
