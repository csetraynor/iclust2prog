#' Gene Names
#'
#' @details These data are the Hugo Symbol and Entrez Id for the Microarray data
#' and CNA experiment. It matches Hugo Symbol and Entrez Id.
#'
#' @name gene_names
#' @aliases gene_names
#' @docType data
#' @return \item{gene_names}{a data frame}
#'
#' @keywords datasets
#' @examples
#' data(gene_names)
#' str(gene_names)
#' @references
#' Pereira, B. et al. Nat. Commun. 7, 2016.
NULL

#' iClust2 GLMNET
#'
#' @details These data are the cross-validation elastic net fit object for
#' only using intClust2 patients.
#'
#' @name iclust2_glmnet
#' @aliases iclust2_glmnet
#' @docType data
#' @return \item{iclust2_glmnet}{a glmnet object}
#'
#' @keywords datasets
#' @examples
#' data(iclust2_glmnet)
#' print(iclust2_glmnet)
#' plot(iclust2_glmnet)
#' @references
#' Noah Simon, Jerome Friedman, Trevor Hastie, Rob Tibshirani (2011).
#'  Regularization Paths for Cox's Proportional Hazards Model via Coordinate
#'  Descent. Journal of Statistical Software, 39(5), 1-13.
#'  URL http://www.jstatsoft.org/v39/i05/.
NULL

#' pooled GLMNET
#'
#' @details These data are the cross-validation elastic net fit object for
#' the pooled dataset.
#'
#' @name pooled_glmnet
#' @aliases pooled_glmnet
#' @docType data
#' @return \item{pooled_glmnet}{a glmnet object}
#'
#' @keywords datasets
#' @examples
#' data(pooled_glmnet)
#' print(pooled_glmnet)
#' plot(pooled_glmnet)
#' @references
#' Noah Simon, Jerome Friedman, Trevor Hastie, Rob Tibshirani (2011).
#'  Regularization Paths for Cox's Proportional Hazards Model via Coordinate
#'  Descent. Journal of Statistical Software, 39(5), 1-13.
#'  URL http://www.jstatsoft.org/v39/i05/.
NULL

#' intClust2 data
#'
#' @details These data are the response variables time/status and only includs
#'  selected variables by glmnet, either using intclust2 only or the pooled dataset.
#'
#' @name intclust2dat
#' @aliases intclust2dat
#' @docType data
#' @return \item{intclust2dat}{a data frame}
#'
#' @keywords datasets
#' @examples
#' data(intclust2dat)
#' str(intclust2dat)
#' @references
#' Pereira, B. et al. Nat. Commun. 7, 2016.
NULL
