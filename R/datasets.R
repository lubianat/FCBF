#' Dengue infected macrophages; gene expression data from GEO study GSE110496
#'
#' Expression data from single cells, from  adengue virus infection study by
#' Zanini et al, #' 2018. The expression was filtered to get cells 12 hours
#' after infection with #' a multiplicity of infection (moi) of 10 (dengue)
#' or uninfected(ctrl). Gene counts were normalized via Bioconductor
#' package "SCNorm".
#'
#' Gene expression has to be discretized for use in FCBF.
#'

#' @name single_cell_dengue_exprs
#' @docType data
#' @usage data(single_cell_dengue_exprs)
#' @format An object of class \code{data.frame}
#' @keywords datasets
#' @references Zanini, F., Pu, S. Y., Bekerman, E., Einav, S., & Quake, S. R.
#' (2018). Single-cell transcriptional dynamics of flavivirus infection.
#' Elife, 7, e32942.
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/29451494}{PubMed}
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110496}{GEO}
#' @examples
#' data(single_cell_dengue_exprs)
#' # Discretize gene expression
#' \dontrun{discrete_expression <- discretize_exprs(single_cell_dengue_exprs)}
#' \dontrun{fcbf_features <- fcbf(discrete_expression,
#'                                single_cell_dengue_annot,
#'                                thresh = 0.05,
#'                                verbose = T)}
"single_cell_dengue_exprs"

#' Dengue infected Macrophages Annotation data
#'
#' Annotation from a dengue virus infection study by Zanini et al,
#' 2018. The expression was filtered to get cells 12 hours
#' after infection with #' a multiplicity of infection (moi) of 10 (dengue)
#' or uninfected(ctrl).This  target annotation is in the same order as
#' single_cell_dengue_exprs and can be used as target for function  fcbf
#'
#' @name single_cell_dengue_annot
#' @docType data
#' @usage data(single_cell_dengue_annot)
#' @format An object of class \code{Factor}
#' @keywords datasets
#' @references Zanini, F., Pu, S. Y., Bekerman, E., Einav, S., & Quake, S. R.
#' (2018). Single-cell transcriptional dynamics of flavivirus infection.
#' Elife, 7, e32942.
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/29451494}{PubMed}
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110496}{GEO}
#' @examples
#' data(single_cell_dengue_exprs)
#' # Discretize gene expression
#' \dontrun{discrete_expression <- discretize_exprs(single_cell_dengue_exprs)}
#' \dontrun{fcbf_features <- fcbf(discrete_expression,
#'                                single_cell_dengue_annot,
#'                                thresh = 0.05,
#'                                verbose = T)}
"single_cell_dengue_annot"
