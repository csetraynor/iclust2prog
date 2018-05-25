library(iclust2prog)
library(dplyr)
library(tidyr)
library(tidyposterior)
theme_set(theme_bw())

###Load data
data("erpos_glmnet")
data("iclust2_glmnet")
data("brca_glmnet")

#Get plots from fit glmnet
plot(erpos_glmnet)
plot(iclust2_glmnet)
plot(brca_glmnet)

#Extract features
erpos_features <- extract_features(erpos_glmnet)
erpos_features$feature[match("1-Sep", erpos_features$feature)] <- "Sep_1"

iclust2_features <- extract_features(iclust2_glmnet)
iclust2_features$feature[match("`GCSAML-AS1_cna`1", iclust2_features$feature)] <- "GCSAML_AS1_cna1"

brca_features <- extract_features(brca_glmnet)
brca_features$feature[match("1-Sep", brca_features$feature)] <- "Sep_1"
brca_features$feature[match("RABEPK_cna-1", brca_features$feature)] <- "RABEPK_cna1"
brca_features$feature[match("ZBTB34_cna-1", brca_features$feature)] <- "ZBTB34_cna1"
brca_features$feature[match("HSPA5_cna-1", brca_features$feature)] <- "HSPA5_cna1"


############ Survival analysis vignette
brca <- readRDS("//mokey.ads.warwick.ac.uk/User41/u/u1795546/Documents/Abstract_WIN_Symposium/Rdata_brca/brca_data.RDS")
intclustdat <- brca[brca$intclust == 2, ]

colnames(intclustdat)[match("`GCSAML-AS1_cna`1", colnames(intclustdat))] <- "GCSAML-AS1_cna1"
colnames(intclustdat)[match("1-Sep", colnames(intclustdat))] <- "Sep_1"
colnames(intclustdat)[match("RABEPK_cna-1", colnames(intclustdat))] <- "RABEPK_cna1"
colnames(intclustdat)[match("ZBTB34_cna-1", colnames(intclustdat))] <- "ZBTB34_cna1"
colnames(intclustdat)[match("HSPA5_cna-1", colnames(intclustdat))] <- "HSPA5_cna1"

intclustdat <-  intclustdat[, unique(c(iclust2_features$feature, erpos_features$feature,  brca_features$feature, "os_months", "os_deceased"))]
rm(brca)
devtools::use_data(intclustdat, overwrite = T)

data("intclustdat")

library(rsample)
set.seed(9666)
mc_samp <- mc_cv(intclustdat, strata = "status", times = 100)

library(purrr)
cens_rate <- function(x) mean(analysis(x)$status == 1)
summary(map_dbl(mc_samp$splits, cens_rate))


############### Create models
colnames(intclustdat)
intclustdat[,brca_features$feature[30:54]]

mc_samp$mod_iclust2 <- pmap(list(mc_samp$splits),
                            function(data){
                              mod_fit(x = data,
                                      form = iclust2_features)
                            })
mc_samp$mod_erpos <- pmap(list(mc_samp$splits),
                          function(data){
                            mod_fit(x = data,
                                    form = erpos_features)
                          })
mc_samp$mod_pooled <- pmap(list(mc_samp$splits),
                          function(data){
                            mod_fit(x = data,
                                    form = brca_features)
                          })

############### Get Brier

mc_samp$brier_erpos <- pmap(list(mc_samp$splits, mc_samp$mod_erpos),
                            function(data, model){
                              get_tdbrier(data = data,
                                          mod = model,
                                          form = erpos_features
                              )
                            })
mc_samp$brier_iclust2 <- pmap(list(mc_samp$splits, mc_samp$mod_iclust2),
                              function(data, model){
                                get_tdbrier(data = data,
                                            mod = model,
                                            form = iclust2_features)
                              })
mc_samp$brier_pooled <- pmap(list(mc_samp$splits, mc_samp$mod_pooled),
                              function(data, model){
                                get_tdbrier(data = data,
                                            mod = model,
                                            form = brca_features)
                              })

###integrate Brier
mc_samp$ibrier_iclust2 <- map_dbl(mc_samp$brier_iclust2, integrate_tdbrier)
mc_samp$ibrier_erpos <- map_dbl(mc_samp$brier_erpos, integrate_tdbrier)
mc_samp$ibrier_pooled <- map_dbl(mc_samp$brier_pooled, integrate_tdbrier)



int_brier <- mc_samp %>%
  select(-matches("^mod"), -starts_with("brier"),  -starts_with("cindex"), -starts_with("roc"), -starts_with("iroc"))

int_brier %>%
  select(-splits) %>%
  gather() %>%
  ggplot(aes(x = statistic, col = model, fill = model)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(legend.position = "top")

int_brier <- perf_mod(int_brier, seed = 6507, iter = 5000, transform = logit_trans)

ggplot(tidy(int_brier)) +
  theme_bw()

ibrier_tab <- summary(tidy(int_brier))
as.data.frame(ibrier_tab) %>% mutate_all(my_round)

require(stargazer)
stargazer(ibrier_tab, type = "latex", summary = FALSE, digits.extra = 3,
          digits = 3, digit.separator = ".",
          title = "Bayesian analysis of resampling iBrier")


comparisons <- contrast_models(
  int_brier,
  list_1 = rep("ibrier_iclust2", 2),
  list_2 = c("ibrier_erpos", "ibrier_pooled"),
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
mc_samp$roc_erpos <- pmap(list(mc_samp$splits, mc_samp$mod_erpos),
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
mc_samp$roc_erpos_relaxed <- pmap(list(mc_samp$splits, mc_samp$mod_erpos_relaxed),
                                  function(data, model){
                                    get_tdroc(data = data,
                                              mod = model
                                    )
                                  })

###integrate ROC
mc_samp$iroc_iclust2 <- map_dbl(mc_samp$roc_iclust2, integrate_tdroc)
mc_samp$iroc_erpos <- map_dbl(mc_samp$roc_erpos, integrate_tdroc)
mc_samp$iroc_iclust2_relaxed <- map_dbl(mc_samp$roc_iclust2_relaxed, integrate_tdroc)
mc_samp$iroc_erpos_relaxed <- map_dbl(mc_samp$roc_erpos_relaxed, integrate_tdroc)

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
  list_1 = rep("iroc_iclust2_relaxed", 3),
  list_2 = c("iroc_iclust2", "iroc_erpos", "iroc_erpos_relaxed"),
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
mc_samp$cindex_erpos <- pmap_dbl(list(mc_samp$splits, mc_samp$mod_erpos),
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
mc_samp$cindex_erpos_relaxed <- pmap_dbl(list(mc_samp$splits, mc_samp$mod_erpos_relaxed),
                                         function(data, model){
                                           get_cindex(data = data,
                                                      mod = model
                                           )
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
  list_1 = rep("cindex_iclust2_relaxed",3),
  list_2 = c("cindex_iclust2", "cindex_erpos_relaxed","cindex_erpos"),
  seed = 4654
)

ggplot(comparisons, size = 0.05) +
  theme_bw()

summary(comparisons, size = 0.05) %>%
  select(contrast, starts_with("pract"))



