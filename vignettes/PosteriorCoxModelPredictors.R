devtools::document()
# library(iclust2prog)
library(glmnet)
library(purrr)
data("intclustdat")
intclustdat <- intclustdat %>%
  dplyr::rename(time = os_months,
                status = os_deceased) %>%
  dplyr::mutate(status = status == 1)
data("iclust2_glmnet")

set.seed(9666)
mc_samp <- bootstraps(intclustdat, strata = "status", times = 100)

iclust2_features <- extract_features(iclust2_glmnet)
iclust2_features$feature <-my_replace(iclust2_features$feature)
colnames(intclustdat) <- my_replace(colnames(intclustdat))

mc_samp$mod_iclust2 <- pmap(list(mc_samp$splits),
                            function(data){
                              mod_fit(x = data,
                                      form = iclust2_features,
                                      iter = 1)
                            })

mc_samp$CoxTab_iclust2 <- pmap(list(mc_samp$mod_iclust2),
                               function(mod){
                                 par.table(mod, log_scale = TRUE)
                               })

mc_samp$GeneTab_iclust2 <- pmap(list(mc_samp$CoxTab_iclust2),
                              function(tab){
                                get_geneTable(tab)}
          )

coeff_Tab <- get_coeff_Tab(mc_samp)

int_coeff <- tidyposterior::perf_mod(coeff_Tab, seed = 6507, iter = 5000)

int_coeff_tab <- tidy(int_coeff) %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(HR = exp(mean(posterior)),
            lower = exp(quantile(posterior, 0.05)),
            upper = exp(quantile(posterior, 0.95)) )
as.data.frame(int_coeff_tab) %>% dplyr::mutate_all(my_round)

devtools::use_data(int_coeff_tab)
