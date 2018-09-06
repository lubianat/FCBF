# Based on implementation  from Rajarshi Guha <rajarshi@presidency.com>
# 13/05/2005
#
# Modified by Tiago Lubiana (22/08/2018)
# Functions to calculate forms of entropy for categorical variables ("factors")

# H(X) - entropy
entropy <- function(x, base = exp(1)) {
  if (!is.factor(x)) {
    stop("For calculating the entropy, the vector must be a factor")
  }
  x <- factor(x)
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

# H(X,Y) - joint entropy
entropy.joint <- function(x, y, base = exp(1)) {
  if (!is.factor(x) || !is.factor(y)) {
    stop("For calculating the joint entropy, the vector x & y must be factors")
  }
  x <- factor(x)
  y <- factor(y)
  t <- table(x, y)
  probabily_of_t <- as.numeric(t / sum(t))
  if (any(probabily_of_t == 0)) {
    probabily_of_t <- probabily_of_t[-which(probabily_of_t== 0)]
  }
  ent <- -1 * sum(probabily_of_t * log(probabily_of_t) / log(base))
  if (is.na(ent)) {
    ent <- 0
  }
  ent
}

# H(X|Y) = H(X,Y) - H(Y) - conditional entropy
entropy.cond <- function(x, y, base = exp(1)) {
  if (!is.factor(x) || !is.factor(y)) {
    stop("For calculating the conditional entropy, the vectors x & y must be factors")
  }
  ent <- entropy.joint(x, y, base) - entropy(y, base)
  if (is.na(ent)) {
    ent <- 0
  }
  ent
}

# Formula for symetrical uncertainty as described in  Yu, L. and Liu, H. , 2003.
SU <- function(x, y, base = exp(1)) {
  if (is.character(x)) {
    x <- as.factor(x)
  }
  if (!is.factor(x) || !is.factor(y)) {
    stop(
      "For calculating the symmetrical uncertainty, the vectors x & y must be factors.
      Using a continuous(numeric) feature set leads to this error."
    )
  }

  Ht <- entropy.joint(x, y, base)
  Hx <- entropy(x, base)
  Hy <- entropy(y, base)
  #cat(Ht,' ',Hx,' ',Hy,'\n')

  # Returns the symmetrical uncertainty value for the vector pair
  2 * (Hy + Hx - Ht) / (Hx + Hy)

  }
