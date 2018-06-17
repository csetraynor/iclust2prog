#' Fit MAP estimates for the Poisson Model PEM.   Important aspects of the modelling #'approach: Dependent variable is the binary event indicator. Adding an intercept is #'recommendable to provide a decent and interpretable reference category for factor #'variables. Use the logarithm of the interval length as an offset.
#'
#'
#' @param
#' d a dataset \cr
#' formula  \cr
#' a log baseline hazard
#' t_dur duration time for individtual ith in the hth interval
#' lambda hazard rate
#' os_event censoring indicator (actual event) for the ith subject in the hth interval
#' @return a MAP fit
#' @export
#' @importFrom rethinking map
#' @importFrom rlang !!
map_pem <- function(dat , formula = tobacco_smoking_history_indicator, offset= t_dur){

  m1 <- rethinking::map(
    alist(
      status ~ dpois(lambda),
      log(lambda) <- a[t_id] + log_tdur,
      a[t_id] ~ dnorm(0, 1)
    ), data=mylist )
  return(m1)
}
