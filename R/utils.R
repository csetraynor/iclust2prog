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

par.table <- function(fit, log_scale = FALSE){ ## needs glm object
  d1 <- summary(fit)$coefficients
  d1 <- as.data.frame(d1)
  dc <- as.data.frame(matrix(confint(fit),ncol=2))
  names(dc) <- c("lower","upper")
  d1 <- cbind(data.frame(Parameter=row.names(d1)),d1,dc)
  rownames(d1) <- NULL
  d1$description <- NA

  if(log_scale){
    d1 <- d1 %>%
      dplyr::rename(Estimate = coef,
                    "se_Estimate" = "se(coef)",
                    HR = "exp(coef)"
      ) %>%
      dplyr::mutate(se = se_Estimate,
                    lower = lower,
                    upper = upper,
                    HR = Estimate) %>%
      dplyr::select(Parameter, Estimate, HR, se_Estimate, se, dplyr::everything()) %>%
      dplyr::select(-Estimate,-se_Estimate)

    d1$description <-"log hazard ratio (relative SE)"
  }else{
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
  }
  d1
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

#' Change - string
#'
#' Change - string
#' @param string \cr
#' @return replacement
#' @export my_replace
my_replace <- function(x){

  x <- gsub("1-Sep", "Sep_1", x)
  x <- gsub("\\-", "_", x)
  gsub("`", "", x)
}


#' Function posterior table
#'
#' Extract features of a Bayesian model
#' @param
#' obj : posterior distribution and contrasts \cr
#' @return a table for publication
#' @export post_tab

post_tab <- function(diff_tab, ibrier_tab){
  ibrier_tab$model <- gsub("ibrier_", "",  ibrier_tab$model)
  ibrier_tab <- ibrier_tab[c(2,1,3),]
  diff_tab <- rbind(rep(NA_real_,3) , diff_tab)
  diff_tab$contrast <- NULL
  cbind(ibrier_tab, diff_tab)
}

get_geneTable <- function(tab){
  genedata <- data.frame(Hugo_Symbol = tab$Parameter,
                         coef = tab$HR)
  if("age_std" %in% genedata$Hugo_Symbol){
    geneTable <- genedata[(-match(c("age_std", "npi"),genedata$Hugo_Symbol)),]
  }else{
    geneTable <- genedata
  }
  geneTable
}

get_coeff_Tab <- function(samp){
  coefficient <- lapply(samp$GeneTab_iclust2, function(x){
    out <- x$coef
  }
  )
  coefficient_tab <- as.data.frame(do.call(rbind, coefficient))
  colnames(coefficient_tab) <- samp$GeneTab_iclust2 $`1`$Hugo_Symbol

  int_coeff <- cbind(samp, coefficient_tab)
  int_coeff <- int_coeff %>%
    dplyr::select(-dplyr::matches("^mod"), -dplyr::matches("^Cox"),  -dplyr::matches("^GeneTab"))
  tibble::as.tibble(int_coeff)
}

get_mat <- function(hlist){
  hmat <- matrix(0L,
                 nrow = length(names(hlist)), ncol = length(names(hlist) ))
  for(ic in seq_along(names(hlist)) ){
    for(iv in seq_along(names(hlist)) ){
      if( any(hlist[[ic]] %in% hlist[[iv]])) hmat[ic,iv] = 1
    }
  }
  colnames(hmat) <- tolower(names(hlist))
  colnames(hmat) <- gsub("\\..*","",colnames(hmat))
  colnames(hmat) <- gsub("_up|_dn","",colnames(hmat))
  for(i in seq_along(hlist)){
    if(sum(hmat[,i]) <= 1)  hmat[i,1] = -99
  }
  hmat = hmat[hmat[,1] != -99,hmat[,1] != -99]
  hmat
}
