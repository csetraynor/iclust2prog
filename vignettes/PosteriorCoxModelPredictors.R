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
  dplyr::mutate(kras = rowSums(.[kras])) %>%
  dplyr::select(splits, id, kras, lef1) %>%
  dplyr::rename(lef1 = MAP1B)

int_gene <- tidyposterior::perf_mod(gene_mod, seed = 6507, iter = 5000)

int_gene_tab <- tidy(int_gene) %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(HR = exp(mean(posterior)),
            lower = exp(quantile(posterior, 0.05)),
            upper = exp(quantile(posterior, 0.95)) )
as.data.frame(int_coeff_tab) %>% dplyr::mutate_all(my_round)

devtools::use_data(int_coeff_tab, overwrite = T)




coeff_tab <- data.frame('Oncogenic Signature' = c('K-ras', 'LEF1', 'PI3K/Akt'),
           HR = c(sum(geneTable[geneTable$Hugo_Symbol == 'NGF',]$HR,
                      geneTable[geneTable$Hugo_Symbol == 'GPR4',]$HR)/2,
                  geneTable[geneTable$Hugo_Symbol == 'MAP1B',]$HR,
                  geneTable[geneTable$Hugo_Symbol == 'NGF',]$HR),
           'lower 0.05' = c(sum(geneTable[geneTable$Hugo_Symbol == 'NGF',]$lower,
                                geneTable[geneTable$Hugo_Symbol == 'GPR4',]$lower)/2,
                            geneTable[geneTable$Hugo_Symbol == 'MAP1B',]$lower,
                            geneTable[geneTable$Hugo_Symbol == 'NGF',]$lower),
           'upper 0.95' = c(sum(geneTable[geneTable$Hugo_Symbol == 'NGF',]$upper,
                                geneTable[geneTable$Hugo_Symbol == 'GPR4',]$upper)/2,
                            geneTable[geneTable$Hugo_Symbol == 'MAP1B',]$upper,
                            geneTable[geneTable$Hugo_Symbol == 'NGF',]$upper) )
coeff_tab <- as.data.frame(coeff_tab) %>%
  mutate_at(vars('Oncogenic.Signature'), as.character) %>%
  mutate_all(my_round)
write.csv(coeff_tab, 'coeff_tab.csv', row.names = FALSE)
