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
#'   _Modeling Survival Data: Extending the Cox Model_.
#'   Springer, New York. ISBN 0-387-98784-3.

mod_fit <- function(x, form, inits, iter = 0, ...) {
  suppressWarnings( coxph(form, data = x, init = inits, control = coxph.control(iter.max = iter), ...) )
}
