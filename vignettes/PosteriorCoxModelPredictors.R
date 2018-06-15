devtools::document()
# library(iclust2prog)
library(glmnet)
library(purrr)
data("ic2dat")

set.seed(9666)
mc_samp <- bootstraps(ic2dat, strata = "status", times = 100)

iclust2_features <- extract_features(iclust2_glmnet)
iclust2_features$feature <-my_replace(iclust2_features$feature)
colnames(ic2dat) <- my_replace(colnames(ic2dat))

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

gene_mod <- coeff_Tab %>%
  dplyr::mutate(kras = rowSums(.[c("NGF", "GPR4")])) %>%
  dplyr::rename(lef1 = MAP1B) %>%
  dplyr::select(splits, id, kras, lef1)


int_gene <- tidyposterior::perf_mod(gene_mod, seed = 6507, iter = 5000)

int_gene_tab <- tidy(int_gene) %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(HR = exp(mean(posterior)),
            lower = exp(quantile(posterior, 0.05)),
            upper = exp(quantile(posterior, 0.95)) )
as.data.frame(int_gene_tab) %>% dplyr::mutate_all(my_round)

devtools::use_data(int_coeff_tab, overwrite = T)

