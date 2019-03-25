#' discretize_gene_supervised
#'
#' Uses several discretizations and selects the one that is best for a given variable (gene)
#' in comparison to a target class by equivocation
#'
#' @param gene  A previously normalized gene expression vector
#' @param target A series of labels matching each of the values in the gene vector
#' @vw_params cuttof parameters for the varying width function. Defaults to 0.25, 0.5 and 0.75
#' @return A data frame with the discretized features in the same order as previously
#' @examples
# data(scDengue)
# exprs <- as.data.frame(SummarizedExperiment::assay(scDengue, 'logcounts'))
# gene <- exprs['ENSG00000166825',]
# infection <- SummarizedExperiment::colData(scDengue)
# target <- infection$infection
# discrete_expression <- as.data.frame(discretize_gene_supervised(gene, target))
# table(discrete_expression)

source('R/entropy.R')

rm(list = grep('split', names(.GlobalEnv), value = TRUE))

split_vector_in_two_by_mean <-
  function(gene_expression_across_samples) {
    gene_expression_across_samples <-
      as.numeric(gene_expression_across_samples)
    return(ifelse(
      gene_expression_across_samples < mean(gene_expression_across_samples),
      'low',
      'high'
    ))
  }

split_vector_in_two_by_median <-
  function(gene_expression_across_samples) {
    gene_expression_across_samples <-
      as.numeric(gene_expression_across_samples)
    return(ifelse(
      gene_expression_across_samples < median(gene_expression_across_samples),
      'low',
      'high'
    ))
  }


split_vector_by_kmeans <-
  function(gene_expression_across_samples,
           centers) {
    set.seed(3)
    gene_expression_across_samples <-
      as.numeric(gene_expression_across_samples)
    tryCatch(
      return(as.factor(
        kmeans(gene_expression_across_samples, centers = centers)$cluster
      )),

      error = function(e) {
        return(as.factor(gene_expression_across_samples * 0))
      }
    )
  }

split_vector_in_three_by_mean_sd <-
  function(gene_expression_across_samples, alpha = 1) {
    gene_expression_across_samples <-
      as.numeric(gene_expression_across_samples)
    avg <- mean(gene_expression_across_samples)
    modified_sd <- sd(gene_expression_across_samples) * alpha
    x <-
      ifelse(
        gene_expression_across_samples < (avg - modified_sd),
        'low',
        ifelse(
          gene_expression_across_samples < (avg + modified_sd),
          'medium',
          'high'
        )
      )
  }

split_vector_in_two_by_vw <-
  function(gene_expression_across_samples,
           cutoff) {
    gene_expression_across_samples <-
      as.numeric(gene_expression_across_samples)
    max_expression = max(gene_expression_across_samples)
    min_expression = min(gene_expression_across_samples)
    break_size = (max_expression - min_expression) * cutoff
    return(ifelse(
      gene_expression_across_samples < (min_expression + break_size),
      'low',
      'high'
    ))
  }

list_of_discs <- grep('split', names(.GlobalEnv), value = TRUE)

discretize_gene_supervised <-
  function(gene,
           target,
           output = 'discretized_vector',
           discs = list_of_discs,
           vw_params = c(0.25, 0.5, 0.75),
           kmeans_centers = c(2, 3, 4),
           sd_alpha = c(0.75, 1, 1.25)) {
    #print(discs)

    su_top <- 0
    disc_top <- 'none'
    vec_top <- gene
    for (disc in discs) {
      if (length(grep('vw', disc)) != 0) {
        f <- get(disc)
        for (param in vw_params) {
          discretized_gene <- f(gene, param)
          su_now <- SU(discretized_gene, target)
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
          su_now <- SU(discretized_gene, target)
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
          su_now <- SU(discretized_gene, target)
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
        su_now <- SU(discretized_gene, target)
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
