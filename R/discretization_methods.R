# Discretization methods from varied sources.



.split_vector_in_two_varying_width <-
  function(gene_expression_across_samples,
           n_of_bins = 3) {
    gene_expression_across_samples <-
      as.numeric(gene_expression_across_samples)
    max_expression = max(gene_expression_across_samples)
    min_expression = min(gene_expression_across_samples)
    break_size = (max_expression - min_expression) / n_of_bins
    return(ifelse(
      gene_expression_across_samples < (min_expression + break_size),
      'low',
      'high'
    ))
  }

.split_vector_in_two_by_mean <-
  function(gene_expression_across_samples) {
    gene_expression_across_samples <-
      as.numeric(gene_expression_across_samples)
    return(ifelse(
      gene_expression_across_samples < mean(gene_expression_across_samples),
      'low',
      'high'
    ))
  }

.split_vector_in_two_by_median <-
  function(gene_expression_across_samples) {
    gene_expression_across_samples <-
      as.numeric(gene_expression_across_samples)
    return(ifelse(
      gene_expression_across_samples < median(gene_expression_across_samples),
      'low',
      'high'
    ))
  }


.split_vector_by_kmeans <-
  function(gene_expression_across_samples,
           centers) {
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

.split_vector_in_three_by_mean_sd <-
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

.split_vector_in_two_by_min_max_thresh <-
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
