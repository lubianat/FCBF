# Based on implementation  from Rajarshi Guha <rajarshi@presidency.com>
# 13/05/2005
#
# Modified by Tiago Lubiana (22/08/2018)
# Functions to calculate forms of entropy for categorical variables ("factors")

# H(X) - entropy

####### get_SU_for_vector_pair ####### 

#'  get_SU_for_vector_pair

#'  Formula for symetrical uncertainty as described in  Yu, L. and Liu, H. , 2003.
#' This functions runs symmetrical uncertainty for two features,
#' returning the score
#'
#' @param x A vector containing a categorical feature
#' @param y A vector containing other categorical feature
#' @param base The base used for the logaritmic function. The default is exp(1) (~2.718)
#' @return A numerical value for the Symetrical Uncertainty score
#' @export
#' @examples
#'  data(scDengue)
#'  exprs <- SummarizedExperiment::assay(scDengue, 'logcounts')
#'  discrete_expression <- as.data.frame(discretize_exprs(exprs))
#'  discrete_expression_gene_1 <- discrete_expression$V1
#'  discrete_expression_gene_2 <- discrete_expression$V2
#'  get_SU_for_vector_pair(discrete_expression_gene_1,discrete_expression_gene_2)

get_SU_for_vector_pair <- function(x, y, base = 2) {
  if (is.character(x)) {
    x <- as.factor(x)
  }
  y <- as.factor(y) 
  if (!is.factor(x) || !is.factor(y)) {
    stop(
      "For calculating the symmetrical uncertainty, the vectors x & y must be factors.
      Using a continuous(numeric) feature set leads to this error."
    )
  }
  Ht <- get_IG_for_vector_pair(x, y, base)
  Hx <- get_entropy_for_vector(x, base)
  Hy <- get_entropy_for_vector(y, base)
  #cat(Ht,' ',Hx,' ',Hy,'\n')
  
  # Returns the symmetrical uncertainty value for the vector pair
  2 * (Hy + Hx - Ht) / (Hx + Hy)
  
}

####### get_IG_for_vector_pair ####### 

#'  Information Gain
#' This functions runs Information Gain for two features,
#' returning the score
#'
#' @param x A vector containing a categorical feature
#' @param y A vector containing other categorical feature
#' @param base The base used for the logaritmic function. The default is exp(1) (~2.718)
#' @return A numerical value for the Information Gain score
#' @export
#' @examples
#'   data(scDengue)
#'   exprs <- SummarizedExperiment::assay(scDengue, 'logcounts')
#'   discrete_expression <- as.data.frame(discretize_exprs(exprs))
#'   discrete_expression_gene_1 <- discrete_expression$V1
#'   discrete_expression_gene_2 <- discrete_expression$V2
#'   get_IG_for_vector_pair(discrete_expression_gene_1,discrete_expression_gene_2)

# Formula for Information Gain

get_IG_for_vector_pair <- function(x, y, base = 2) {
  if (is.character(x)) {
    x <- as.factor(x)
  }
  if (!is.factor(x) || !is.factor(y)) {
    stop(
      "For calculating the information gain, the vectors x & y must be factors.
      Using a continuous(numeric) feature set leads to this error."
    )
  }
  Ht <- get_joint_entropy_for_vectors(x, y, base)
  Hx <- get_entropy_for_vector(x, base)
  Hy <- get_entropy_for_vector(y, base)
  # Returns the information gain for the pair
  IG <- (Hy + Hx - Ht)
  IG
}


####### get_entropy_for_vector ####### 

get_entropy_for_vector <- function(x, base = 2) {
  if (!is.factor(x)) {
    stop("For calculating the entropy, the vector must be a factor")
  }
  t <- table(x)
  probabily_of_t <- t / sum(t)
  if (any(t == 0)) {
    probabily_of_t <- probabily_of_t[-which(t == 0)]
  }
  ent <- -1 * sum(probabily_of_t * log(probabily_of_t) / log(base))
  if (is.na(ent)) {
    ent <- 0
  }
  ent
}

####### get_joint_entropy_for_vectors ####### 
# H(X,Y) - joint entropy
get_joint_entropy_for_vectors <- function(x, y, base = 2) {
  if (!is.factor(x) || !is.factor(y)) {
    stop("For calculating the joint entropy, the vector x & y must be factors")
  }
  t <- table(x, y)
  probabily_of_t <- as.numeric(t / sum(t))
  if (any(probabily_of_t == 0)) {
    probabily_of_t <- probabily_of_t[-which(probabily_of_t == 0)]
  }
  ent <- -1 * sum(probabily_of_t * log(probabily_of_t) / log(base))
  if (is.na(ent)) {
    ent <- 0
  }
  ent
}

####### get_conditional_entropy_for_vectors ####### 
# H(X|Y) = H(X,Y) - H(Y) - conditional entropy
get_conditional_entropy_for_vectors <- function(x, y, base = 2) {
  if (!is.factor(x) || !is.factor(y)) {
    stop("For calculating the conditional entropy, the vectors x & y must be factors")
  }
  ent <- get_joint_entropy_for_vectors(x, y, base) - get_entropy_for_vector(y, base)
  if (is.na(ent)) {
    ent <- 0
  }
  ent
}






####### Deprecated Functions #####

## SU
# Formula for symetrical uncertainty as described in  Yu, L. and Liu, H. , 2003.
#' This functions runs symmetrical uncertainty for two features,
#' returning the score
#'
#' @param x A vector containing a categorical feature
#' @param y A vector containing other categorical feature
#' @param base The base used for the logaritmic function. The default is exp(1) (~2.718)
#' @return A numerical value for the Symetrical Uncertainty score
#' @name SU-deprecated
#' @seealso \code{\link{FCBF-deprecated}}
#' @keywords internal
NULL

#' @rdname FCBF-deprecated
#' @section \code{SU}:
#' For \code{SU}, use \code{\link{get_SU_for_vector_pair}}.
#'
#' @export
SU <- function(x, y, base = 2) {
  .Deprecated("SU")
  "SU function is no longer available"
}


#'  Information Gain
#' This functions runs Information Gain for two features,
#' returning the score
#'
#' @param x A vector containing a categorical feature
#' @param y A vector containing other categorical feature
#' @param base The base used for the logaritmic function. The default is exp(1) (~2.718)
#' @return A numerical value for the Information Gain score
#' @name IG-deprecated
#' @seealso \code{\link{FCBF-deprecated}}
#' @keywords internal
NULL

#' @rdname FCBF-deprecated
#' @section \code{IG}:
#' For \code{IG}, use \code{\link{get_IG_for_vector_pair}}.
#'
IG <- function(x, y, base = 2) {
  .Deprecated("IG")
  "IG function is no longer available"
}

