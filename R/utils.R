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
  optimal.coef <- as.matrix(glmnet::coef.cv.glmnet(mod, s = "lambda.min"))
  optimal.coef <- as.data.frame(optimal.coef)
  colnames(optimal.coef) <- "coef"
  optimal.coef <- tibble::rownames_to_column(optimal.coef, var = "feature")
  optimal.coef <-  optimal.coef[optimal.coef$coef != 0,]
  return(optimal.coef)
}

#' Function par table
#'
#' Extract features of a classical fit model
#' @param
#' obj : survival coxph fit. \cr
#' @return a coxph fit object
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @import prodlim
#' @import survival
#' @author Sahota Tarj - caret

par.table <- function(fit){ ## needs glm object
  d1 <- summary(fit)$coefficients
  d1 <- as.data.frame(d1)
  dc <- as.data.frame(matrix(confint(fit),ncol=2))
  names(dc) <- c("lower","upper")
  d1 <- cbind(data.frame(Parameter=row.names(d1)),d1,dc)
  rownames(d1) <- NULL
  d1$description <- NA

  d1 <- d1 %>%
    dplyr::rename(Estimate = coef,
                  "se_Estimate" = "se(coef)",
                  HR = "exp(coef)"
    ) %>%
    dplyr::mutate(se = exp(se_Estimate),
                  lower = exp(lower),
                  upper = exp(upper)) %>%
    dplyr::select(Parameter, Estimate, HR, se_Estimate, se, dplyr::everything()) %>%
    dplyr::select(-Estimate,-se_Estimate)

  d1$description <-"Hazard ratio (relative SE)"

  return(d1)
}



#' Map gene_entrez_id
#'
#' Mapping the genes
#' @param Hugo_Symbol \cr
#' @return gene_entrez_id
#' @export mg

mg <- function(x){
  gene_names$Entrez_Gene_Id[match(x, gene_names$Hugo_Symbol)]
}

#' Round data frame
#'
#' Round data frame
#' @param decimals \cr
#' @return rounded
#' @export my_round
my_round <- function(x){
  if(is.character(x)){
    x
  }else{
    round(x,3)
  }
}
