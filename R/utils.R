#' Extract features
#'
#' Extract relevant features from a glmnet cox fit object.
#'
#' @param mod Coxph model object fitted with coxph (survival).
#' @return Features
#' @seealso [glmnet]
#' @keywords glmnet
#' @author Carlos S Traynor
#' @references
#'
#'  Noah Simon, Jerome Friedman, Trevor Hastie, Rob Tibshirani (2011).
#'   Regularization Paths for Cox's Proportional Hazards Model via
#'  Coordinate Descent. Journal of Statistical Software, 39(5), 1-13.
#'  URL  http://www.jstatsoft.org/v39/i05/.
#' @export extract_features

extract_features <- function(mod){
  # find optimised lambda
  optimal.coef <- as.matrix(glmnet::coef(mod, s = "lambda.min"))
  optimal.coef <- as.data.frame(optimal.coef)
  colnames(optimal.coef) <- "coef"
  optimal.coef <- tibble::rownames_to_column(optimal.coef, var = "feature")
  optimal.coef <-  optimal.coef[optimal.coef$coef != 0,]
  return(optimal.coef)
}
