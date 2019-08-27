#' @importMethodsFrom SummarizedExperiment assay colData
NULL
# Adapted by Tiago Lubiana (tiago.lubiana.alves@usp.br)
# From Rajarshi Guha <rajarshi@presidency.com>'s implementation of FCBF
# Implementation of the Fast Correlation Based Filter
# described in
#
# Yu, L. and Liu, H.; Feature Selection
# for High-Dimensional Data A Fast Correlation Based Filter Solution,
# Proc. 20th Intl. Conf. Mach. Learn. (ICML-2003), Washington DC, 2003
#
# Require the functions in entropy.R which can be obtained from
# http://blue.chem.psu.edu/~rajarshi/code/R/#entropy

source('R/entropy.R')

.get.next.elem <- function(s, first_prime) {
  index <- which(s == first_prime)
  if (index == length(s)) {
    NA
  } else {
    s[index + 1]
  }
}

#' Fast Correlation Based Filter function.
#'
#' This functions allows selection of variables from a feature table
#' of discrete/categorial variables and a target class.
#' The function is based on the algorithm described in
#' Yu, L. and Liu, H.; Feature Selection
#' for High-Dimensional Data A Fast Correlation Based Filter Solution,
#' Proc. 20th Intl. Conf. Mach. Learn. (ICML-2003), Washington DC, 2003
#'
#' Obs: For gene expression, you will need to run discretize_exprs first
#'
#' @param x A table of features (samples in rows, variables in columns, and each observation in each cell)
#' @param y A target vector, factor containing classes of the observations. Note: the
#' observations must be in the same order as the parameter x
#' @param thresh A threshold for the minimum correlation (as determined by symettrical uncertainty)
#' between each variable and the class. Defaults to 0.25.
#' Note: this might drastically change the number of selected features.
#' @param samples_in_rows A flag for the case in which samples are in rows and variables/genes in columns. Defaults to FALSE.
#' @param verbose Adds verbosity. Defaults to FALSE.
#' @param n_genes Sets the number of genes to be selected in the first part of the algorithm.
#' If left unchanged, it defaults to NULL and the thresh parameter is used.
#' Caution: it overrides the thresh parameter altogether.
#' @param balance_classes Balances number of instances in the target vector y by sampling the number of instances in the
#' minor class from all others. The number of samplings is controlled by resampling_number. Defaults to FALSE.

#' @return Returns a data frame with the selected features index (first row) and their symmetrical uncertainty values regarding the class (second row). Variable names are present in rownames
#' @export
#' @examples
#'  data(scDengue)
#'  exprs <- SummarizedExperiment::assay(scDengue, 'logcounts')
#'  discrete_expression <- as.data.frame(discretize_exprs(exprs))
#'  head(discrete_expression[,1:4])
#'  infection <- SummarizedExperiment::colData(scDengue)
#'  target <- infection$infection
#'  fcbf(discrete_expression,target, thresh = 0.05, verbose = TRUE)
#'  fcbf(discrete_expression,target, n_genes = 100)

fcbf <-
  function(x,
           y,
           thresh = 0.25,
           n_genes = NULL,
           verbose = FALSE,
           samples_in_rows = FALSE,
           balance_classes = FALSE) {
    if (!samples_in_rows) {
      x <- t(x)
    }
    if (!balance_classes) {
      x <- data.frame(x)
      nvar <- ncol(x)
      if (verbose) {
        message("Calculating symmetrical uncertainties")
      }
      su_ic <- apply(x, 2, function(xx, yy) {
        SU(xx, yy)
      }, y)

      if (length(n_genes)) {
        thresh <- sort(su_ic, decreasing = TRUE)[n_genes - 1]
      }

      s_prime <-
        data.frame(f = (seq_len(nvar))[which(su_ic >= thresh)], su = su_ic[which(su_ic >= thresh)])

      s_prime <- s_prime[sort.list(s_prime$su, decreasing = TRUE), ]

      # s_prime is the list of selected features ranked by su_ic
      s_prime <- s_prime[, 1]

      if (length(s_prime) == 1) {
        s_prime
      } else if (length(s_prime) == 0) {
        stop("No prospective features for this threshold level. Threshold: ",
             thresh)
      }

      print(paste('Number of prospective features = ', length(s_prime)))


      first_prime  <- s_prime[1]
      cnt <- 1
      while (TRUE) {
        if (verbose) {
          cat("Round ")
          cat(cnt, "\n")
          cnt <- cnt + 1
          print(paste(
            'first_prime  round ( |s_prime| =',
            length(s_prime),
            ')' ,
            sep =
              ' '
          ))
        }

        next_prime <- .get.next.elem(s_prime, first_prime)
        if (!is.na(next_prime)) {
          while (TRUE) {
            prime_to_be_compared <- next_prime
            su1 = SU(x[, first_prime], x[, next_prime])
            su2 = SU(x[, next_prime], y)
            if (su1 > su2) {
              next_prime <- .get.next.elem(s_prime, next_prime)
              s_prime <-
                s_prime[-which(s_prime == prime_to_be_compared)]
              if (verbose) {
                cat("  ",
                    su1,
                    " ",
                    su2,
                    " ",
                    "Removed feature ",
                    prime_to_be_compared,
                    "\n")
              }
            }
            else {
              next_prime <- .get.next.elem(s_prime, next_prime)
            }
            if (is.na(next_prime)) {
              break
            }
          }
        }
        first_prime  <- .get.next.elem(s_prime, first_prime)

        if (is.na(first_prime)) {
          break
        }
      }
      if (length(s_prime) > 1) {
        suvalues <- apply(x[, s_prime], 2, function(xx, yy) {
          SU(xx, yy)
        }, y)
        data.frame(index = s_prime, SU = suvalues)
      } else {
        data.frame(index = s_prime, SU = SU(x[, s_prime], y))
      }
    }
    else{
      instances_in_minor_class <- min(table(y))
      n_x <- as.data.frame(cbind(x, y))
      final_x <- data.frame(n_x[1, ])
      final_x <- final_x[-1, ]
      for (i in levels(as.factor(n_x[, ncol(n_x)]))) {
        n_x_i <- n_x[n_x[, ncol(n_x)] == i, ]
        n_x_i <-
          n_x_i[sample(seq_len(length.out = nrow(n_x_i)), instances_in_minor_class), ]
        final_x <- rbind(n_x_i, final_x)
      }
      final_x$y <- NULL
      fcbf(t(final_x),
           y,
           thresh,
           verbose,
           samples_in_rows,
           balance_classes = FALSE)

    }
  }
