#' Dengue infected macrophages; gene expression data from GEO study GSE110496
#'
#' Expression data from single cells, from  adengue virus infection study by
#' Zanini et al, #' 2018. The expression was filtered to get cells 12 hours
#' after infection with #' a multiplicity of infection (moi) of 1 (dengue)
#' or uninfected(ctrl). Gene counts were normalized via Bioconductor
#' package "SCNorm".
#'
#' Gene expression has to be discretized for use in FCBF.
#'

#' @name scDengue
#' @docType data
#' @usage data(scDengue)
#' @format An object of class \code{SingleCellExperiment}
#' @keywords datasets, dengue, single-cell
#' @references Zanini, F., Pu, S. Y., Bekerman, E., Einav, S., & Quake, S. R.
#' (2018). Single-cell transcriptional dynamics of flavivirus infection.
#' Elife, 7, e32942.
#' \href{https://www.ncbi.nlm.nih.gov/pubmed/29451494}{PubMed}
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110496}{GEO}
@examples
data(scDengue)
infection <- scDengue@colData
target <- infection@listData$infection
exprs <- scDengue@assays$data$logcounts
# Discretize gene expression
discrete_expression <- as.data.frame(discretize_exprs(exprs))
fcbf_features <- fcbf(discrete_expression,
                              target,
                              thresh = 0.05,
                              verbose = TRUE)
"scDengue"
