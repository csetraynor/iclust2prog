#' Model fitting
#'
#'
#' Fits coxph model.
#'
#' @param x data
#' @param mod Coxph model object fitted with coxph (survival).
#' @return mod
#' @seealso [coxph]
#' @keywords coxph
#'
#' @author Carlos S Traynor
#' @references
#'
#'  Terry M. Therneau and Patricia M. Grambsch (2000).
#'   Modeling Survival Data: Extending the Cox Model.
#'   Springer, New York. ISBN 0-387-98784-3.
#'@export mod_fit

mod_fit <- function(x, form, iter = 0, inits = NA_character_,...) {
  x <- rsample::assessment(x)

  X <- x[,form$feature]

  # if(is.character(inits)){
  #   inits = rep(0, length(features))
  # }

  suppressWarnings(coxph(Surv(x$time, x$status)~ . , data = X, init = form$coef, control = coxph.control(iter.max = iter) ) )
}
