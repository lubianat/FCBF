#' @importMethodsFrom SummarizedExperiment assay colData
#' @importFrom stats kmeans median sd
NULL

source('R/discretization_methods.R')

#' discretize_gene_supervised
#'
#' Uses several discretizations and selects the one that is best for a given variable (gene)
#' in comparison to a target class by equivocation
#' Note that set.seed() should be used for reproducing the results. The inner
#' kmeans #' function would, otherwise, provide different results each time.
#'
#' Note that a seed for random values has to bew set for reproducibility.
#' Otherwise, the "kmeans" value might vary from iteration to iteration.
#'
#' @param gene  A previously normalized gene expression vector
#' @param target A series of labels matching each of the values in the gene vector
#' @param output If it is equal to 'discretized_vector', the output is the vector. I it is 'su', returns a dataframe. Defaults to 'discretized_vector'
#' @param vw_params cuttof parameters for the varying width function. Defaults to 0.25, 0.5 and 0.75
#' @param kmeans_centers Numeric vector with the number of centers to use for kmeans. Defaults to 2, 3 and 4
#' @param discs Defaults to   c(
#' ".split_vector_in_two_by_median",
#' split_vector_in_two_by_mean",
#' ".split_vector_by_kmeans",
#' ".split_vector_in_three_by_mean_sd",
#' ".split_vector_in_two_by_vw")
#' @param sd_alpha Parameter for adusting the 'medium' level of the mean +- sd discretization. Defaults to sd_alpha = c(0.75, 1, 1.25))
#' @export
#' @return A data frame with the discretized features in the same order as previously
#' @examples
#'  data(scDengue)
#'  exprs <- as.data.frame(SummarizedExperiment::assay(scDengue, 'logcounts'))
#'  gene <- exprs['ENSG00000166825',]
#'  infection <- SummarizedExperiment::colData(scDengue)
#'  target <- infection$infection
#'  set.seed(3)
#'  discrete_expression <- as.data.frame(discretize_gene_supervised(gene, target))
#'  table(discrete_expression)
#'
discretize_gene_supervised <-
  function(gene,
           target,
           output = 'discretized_vector',
           discs =   c(
             ".split_vector_in_two_by_median",
             ".split_vector_in_two_by_mean",
             ".split_vector_by_kmeans",
             ".split_vector_in_three_by_mean_sd",
             ".split_vector_in_two_by_min_max_thresh"
           )
           ,
           vw_params = c(0.25, 0.5, 0.75),
           kmeans_centers = c(2, 3, 4),
           sd_alpha = c(0.75, 1, 1.25)) {
    #print(discs)

    su_top <- 0
    disc_top <- 'none'
    vec_top <- gene
    for (disc in discs) {
      if (length(grep('min_max_thresh', disc)) != 0) {
        f <- get(disc)
        for (param in vw_params) {
          discretized_gene <- f(gene, param)
          su_now <- get_SU_for_vector_pair(discretized_gene, target)
          #print(paste(su_now,disc))
          if (su_now > su_top) {
            su_top <- su_now
            disc_top <- paste0(disc, '_', param)
            vec_top <- discretized_gene
          }

        }
      } else if (length(grep('kmeans', disc)) != 0) {
        f <- get(disc)
        for (param in kmeans_centers) {
          discretized_gene <- f(gene, param)
          su_now <- get_SU_for_vector_pair(discretized_gene, target)
          #print(paste(su_now,disc))
          if (su_now > su_top) {
            su_top <- su_now
            disc_top <- paste0(disc, '_', param)
            vec_top <- discretized_gene
          }

        }
      } else if (length(grep('sd', disc)) != 0) {
        f <- get(disc)
        for (param in sd_alpha) {
          discretized_gene <- f(gene, param)
          su_now <- get_SU_for_vector_pair(discretized_gene, target)
          #print(paste(su_now,disc))
          if (su_now > su_top) {
            su_top <- su_now
            disc_top <- paste0(disc, '_', param)
            vec_top <- discretized_gene
          }
        }
      } else {
        f <- get(disc)
        discretized_gene <- f(gene)
        su_now <- get_SU_for_vector_pair(discretized_gene, target)
        #print(paste(su_now,disc))
        if (su_now > su_top) {
          su_top <- su_now
          disc_top <- disc
          vec_top <- discretized_gene
        }

      }

    }
    if (output == 'discretized_vector') {
      return(vec_top)
    } else if (output == 'su') {
      return(data.frame(su_ic = su_top, disc = disc_top))
    }
  }
