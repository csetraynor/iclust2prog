library(iclust2prog)
library(dplyr)
library(tidyr)
library(tidyposterior)
theme_set(theme_bw())

###Load data
data("pooled_glmnet")
data("iclust2_glmnet")

#Get plots from fit glmnet
plot(pooled_glmnet)
plot(iclust2_glmnet)

#Extract features
pooled_features <- extract_features(pooled_glmnet)
iclust2_features <- extract_features(iclust2_glmnet)

############ Survival analysis vignette
data("intclustdat")

intclustdat <- intclustdat %>%
  rename(time = os_months,
         status = os_deceased,
         "Sep_1" = "1-Sep") %>%
  mutate(status = status == 1)

library(rsample)
set.seed(9666)
mc_samp <- mc_cv(intclustdat, strata = "status", times = 100)

library(purrr)
cens_rate <- function(x) mean(analysis(x)$status == 1)
summary(map_dbl(mc_samp$splits, cens_rate))


################ Create formulas
paste0(iclust2_features$feature, collapse = "+")
paste0(pooled_features$feature, collapse = "+")

iclust2 <- as.formula(Surv(time, status) ~ npi+MTCP1NB+C9orf95+DB451841+MAP1B+CSAG1+NDUFA4L2+NGF+BOLA1+CAMKK1+AK055323+GPR4+age_std)
pooled <- as.formula(Surv(time, status) ~ npi+SEMA3B+COX7B2+STIP1+BF511988+ADCK2+STAT5A+SRFBP1+CEACAM5+BQ007025+LIPH+MAP1B+SYTL1+POLR3A+ANGPT2+KIAA1191+RBP7+LARP1+GPR87+RMND5B+EFEMP2+C2orf74+PPIL3+Sep_1+PDLIM7+DDX41+OMD+USP30+PAOX+FGD3+BG218808+WDR67+CPNE1+TNIP2+ENO1+GSK3B+KIAA1161+RALGAPB+ORAI3+C10orf27+BX102609+CATSPERB+ENC1+AI783972+VANGL1+age_std+TRPM3_cna2)

############### Create models
iclust2_inits <- iclust2_features$coef
pooled_inits <- pooled_features$coef

mc_samp$mod_iclust2 <- pmap(list(mc_samp$splits),
                            function(data){
                              mod_fit(x = data,
                                      form = iclust2,
                                      inits = iclust2_inits,
                                      iter = 0)
                            })
mc_samp$mod_pooled <- pmap(list(mc_samp$splits),
                            function(data){
                              mod_fit(x = data,
                                      form = pooled,
                                      inits = pooled_inits,
                                      iter = 0)
                            })
mc_samp$mod_pooled_relaxed <- pmap(list(mc_samp$splits),
                           function(data){
                             mod_fit(x = data,
                                     form = pooled,
                                     inits = pooled_inits,
                                     iter = 14)
                           })
mc_samp$mod_iclust2_relaxed <- pmap(list(mc_samp$splits),
                            function(data){
                              mod_fit(x = data,
                                          form = iclust2,
                                          inits = iclust2_inits,
                                          iter = 14)
                            })

############### Get Brier

mc_samp$brier_pooled <- pmap(list(mc_samp$splits, mc_samp$mod_pooled),
                              function(data, model){
                                get_tdbrier(data = data,
                                            mod = model,
                                            inits = pooled_inits,
                                            iters = 0
                                            )
                              })
mc_samp$brier_iclust2 <- pmap(list(mc_samp$splits, mc_samp$mod_iclust2),
                              function(data, model){
                                get_tdbrier(data = data,
                                            mod = model,
                                            inits = iclust2_inits,
                                            iters = 0)
                              })
mc_samp$brier_pooled_relaxed <- pmap(list(mc_samp$splits, mc_samp$mod_pooled_relaxed),
                                      function(data, model){
                                        get_tdbrier(data = data,
                                                    mod = model,
                                                    inits = pooled_inits,
                                                    iters = 14
                                        )
                                      })

mc_samp$brier_iclust2_relaxed <- pmap(list(mc_samp$splits, mc_samp$mod_iclust2_relaxed),
                             function(data, model){
                               get_tdbrier(data = data,
                                           mod = model,
                                           inits = iclust2_inits,
                                           iters = 14
                               )
                             })

###integrate Brier
mc_samp$ibrier_iclust2 <- map_dbl(mc_samp$brier_iclust2, integrate_tdbrier)
mc_samp$ibrier_pooled <- map_dbl(mc_samp$brier_pooled, integrate_tdbrier)
mc_samp$ibrier_pooled_relaxed <- map_dbl(mc_samp$brier_pooled_relaxed, integrate_tdbrier)
mc_samp$ibrier_iclust2_relaxed <- map_dbl(mc_samp$brier_iclust2_relaxed, integrate_tdbrier)



int_brier <- mc_samp %>%
  select(-matches("^mod"), -starts_with("brier"),  -starts_with("cindex"), -starts_with("roc"), -starts_with("iroc"))

int_brier %>%
  select(-splits) %>%
  gather() %>%
  ggplot(aes(x = statistic, col = model)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(legend.position = "top")

int_brier <- perf_mod(int_brier, seed = 6507, iter = 5000, transform = logit_trans)

ggplot(tidy(int_brier)) +
  theme_bw()

ibrier_tab <- summary(tidy(int_brier))

my_round <- function(x){
  if(is.character(x)){
    x
  }else{
    round(x,3)
  }
}

as.data.frame(ibrier_tab) %>% mutate_all(my_round)

require(stargazer)
stargazer(ibrier_tab, type = "latex", summary = FALSE, digits.extra = 3,
          digits = 3, digit.separator = ".",
          title = "Bayesian analysis of resampling iBrier")


comparisons <- contrast_models(
  int_brier,
  list_1 = rep("ibrier_iclust2_relaxed", 3),
  list_2 = c("ibrier_iclust2", "ibrier_pooled_relaxed", "ibrier_pooled"),
  seed = 4654
)

ggplot(comparisons, size = 0.01) +
  theme_bw()

summary(comparisons, size = 0.01) %>%
  select(contrast, starts_with("pract"))


############### Get ROC

mc_samp$roc_iclust2 <- pmap(list(mc_samp$splits, mc_samp$mod_iclust2),
                              function(data, model){
                                get_tdroc(data = data,
                                            mod = model)
                              })
mc_samp$roc_pooled <- pmap(list(mc_samp$splits, mc_samp$mod_pooled),
                             function(data, model){
                               get_tdroc(data = data,
                                           mod = model
                               )
                             })
mc_samp$roc_iclust2_relaxed <- pmap(list(mc_samp$splits, mc_samp$mod_iclust2_relaxed),
                           function(data, model){
                             get_tdroc(data = data,
                                       mod = model
                             )
                           })

###integrate ROC
mc_samp$iroc_iclust2 <- map_dbl(mc_samp$roc_iclust2, integrate.tdroc)
mc_samp$iroc_pooled <- map_dbl(mc_samp$roc_pooled, integrate.tdroc)
mc_samp$iroc_iclust2_relaxed <- map_dbl(mc_samp$roc_iclust2_relaxed, integrate.tdroc)


library(dplyr)
int_roc <- mc_samp %>%
  select(-matches("^mod"), -starts_with("brier"),  -starts_with("ibrier"), -starts_with("roc"), -starts_with("cindex"))


int_roc %>%
  select(-splits) %>%
  gather() %>%
  ggplot(aes(x = statistic, col = model)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(legend.position = "top")


library(tidyposterior)
int_roc <- perf_mod(int_roc, seed = 6507, iter = 5000)

ggplot(tidy(int_roc)) +
  theme_bw()

roc_tab <- summary(tidy(int_roc))
as.data.frame(roc_tab) %>% mutate_all(my_round)

comparisons <- contrast_models(
  int_roc,
  list_1 = rep("iroc_pooled", 2),
  list_2 = c("iroc_iclust2", "iroc_iclust2_relaxed"),
  seed = 4654
)

ggplot(comparisons, size = 0.05) +
  theme_bw()

summary(comparisons, size = 0.05) %>%
  select(contrast, starts_with("pract"))


######Concordance Index
mc_samp$cindex_iclust2 <- pmap_dbl(list(mc_samp$splits, mc_samp$mod_iclust2),
                            function(data, model){
                              get_cindex(data = data,
                                        mod = model)
                            })
mc_samp$cindex_pooled <- pmap_dbl(list(mc_samp$splits, mc_samp$mod_pooled),
                           function(data, model){
                             get_cindex(data = data,
                                       mod = model
                             )
                           })
mc_samp$cindex_iclust2_relaxed <- pmap_dbl(list(mc_samp$splits, mc_samp$mod_iclust2_relaxed),
                                   function(data, model){
                                     get_cindex(data = data,
                                                mod = model)
                                   })

cindex_est <- mc_samp %>%
  select(-matches("^mod"), -starts_with("brier"),  -starts_with("ibrier"), -starts_with("roc"), -starts_with("iroc"))

cindex_est %>%
  select(-splits) %>%
  gather() %>%
  ggplot(aes(x = statistic, col = model)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(legend.position = "top")

mc_cindex_est <- perf_mod(cindex_est, seed = 6507, iter = 5000)

ggplot(tidy(mc_cindex_est)) +
  theme_bw()

ci_tab <- summary(tidy(mc_cindex_est))
as.data.frame(ci_tab) %>% mutate_all(my_round)

comparisons <- contrast_models(
  mc_cindex_est,
  list_1 = rep("cindex_pooled",2),
  list_2 = c("cindex_iclust2", "cindex_iclust2_relaxed"),
  seed = 4654
)

ggplot(comparisons, size = 0.05) +
  theme_bw()

summary(comparisons, size = 0.05) %>%
  select(contrast, starts_with("pract"))



