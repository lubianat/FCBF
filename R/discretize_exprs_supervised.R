#' @importMethodsFrom SummarizedExperiment assay colData
#' @import pbapply
NULL
#' supervised_disc_df
#'
#' Uses several discretizations and selects the one that is best for a given variable (gene)
#' in comparison to a target class by equivocation
#'
#' @param expression_table  A previously normalized expression table
#' @param target A series of labels matching each of the values in the gene vector
#' (genes in rows, cells/samples in columns)
#' @param parallel Set calculations in parallel. May be worth it if the number of rows and columns is really large. Do watchout for memory overload.
#' @export
#' @return A data frame with the discretized features in the same order as previously
#' @examples
#' data(scDengue)
#' exprs <- as.data.frame(SummarizedExperiment::assay(scDengue, 'logcounts'))
#' exprs <- exprs [1:200, 1:120]
#' infection <- SummarizedExperiment::colData(scDengue)
#' target <- infection$infection
#' discrete_expression <- as.data.frame(discretize_exprs_supervised(exprs,target))
#' fcbf(discrete_expression,target, thresh = 0.05, verbose = TRUE)

discretize_exprs_supervised <-
  function(expression_table, target, parallel = FALSE) {
    if (parallel) {
      ncores <- parallel::detectCores() - 2
      cl <- parallel::makeCluster(ncores)
      discrete_expression <-
        parallel::parApply(cl,
                           expression_table,
                           1,
                           discretize_gene_supervised,
                           target)
      parallel::stopCluster(cl)
    }
    else{
      discrete_expression <-
        pbapply(expression_table, 1, discretize_gene_supervised, target)
    }
    discrete_expression <- as.data.frame(t(discrete_expression))
    rownames(discrete_expression) <- rownames(expression_table)
    return(discrete_expression)
  }
