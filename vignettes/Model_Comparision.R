if(!require("devtools")) install.packages("devtools")
if(!require("iclust2prog")) devtools::install_github("iclust2prog")

library(iclust2prog)
library(glmnet)
library(purrr)
library(dplyr)
library(tidyr)
library(rsample)
library(tidyposterior)
theme_set(theme_bw())

###Load data
data("erpos_glmnet")
data("iclust2_glmnet")

#Plot glmnet fits
glmnet::plot.cv.glmnet(erpos_glmnet)
glmnet::plot.cv.glmnet(iclust2_glmnet_72)

#Extract features, see my_replace in utils
erpos_features <- extract_features(erpos_glmnet)
erpos_features$feature <-my_replace(erpos_features$feature)

iclust2_features <- extract_features(iclust2_glmnet)
iclust2_features$feature <-my_replace(iclust2_features$feature)

############ Survival analysis vignette
# brca <- readRDS("/home/mtr/rfactory/brca_data.RDS")
# cna <- readRDS("/home/mtr/rfactory/cna_expression.RDS")
# brca <- cbind(brca, cna)
# intclustdat <- brca[brca$intclust == 2, ]
# #
#  colnames(intclustdat) <- my_replace(colnames(intclustdat))
# #
# intclustdat <-  intclustdat[, unique(c(iclust2_features$feature, iclust2_features_72$feature, erpos_features$feature,  "os_months", "os_deceased"))]
# rm(brca)
# rm(cna)
# devtools::use_data(intclustdat, overwrite = T)

data("intclustdat")
intclustdat <- intclustdat %>%
  dplyr::rename(time = os_months,
         status = os_deceased) %>%
  dplyr::mutate(status = status == 1)


set.seed(9666)
mc_samp <- mc_cv(intclustdat, strata = "status", times = 100)


cens_rate <- function(x) mean(analysis(x)$status == 1)
summary(map_dbl(mc_samp$splits, cens_rate))


############### Create models
mc_samp$mod_iclust2 <- pmap(list(mc_samp$splits),
                            function(data){
                            mod_fit(x = data,
                                      form = iclust2_features,
                                    iter = 5)
                            })
mc_samp$mod_erpos <- pmap(list(mc_samp$splits),
                          function(data){
                            mod_fit(x = data,
                                    form = erpos_features,
                                    iter = 5)
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

###integrate Brier
mc_samp$iClust2 <- map_dbl(mc_samp$brier_iclust2, integrate_tdbrier)
mc_samp$'ER+/HER2-' <- map_dbl(mc_samp$brier_erpos, integrate_tdbrier)
mc_samp$Reference <- map_dbl(mc_samp$brier_erpos, integrate_tdbrier_reference)



int_brier <- mc_samp %>%
  dplyr::select(-matches("^mod"), -starts_with("brier"),  -starts_with("cindex"), -starts_with("roc"), -starts_with("iroc"))

int_brier %>%
  dplyr::select(-splits) %>%
  gather() %>%
  ggplot(aes(x = statistic, col = model, fill = model)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(legend.position = "top")

int_brier <- perf_mod(int_brier, seed = 6507, iter = 5000, transform = logit_trans)

pdf <- ggplot(tidy(int_brier)) +
  theme_bw() +
  ylab("Posterior probability for BS")
pdf
ibrier_tab <- tidy(int_brier) %>%
  group_by(model) %>%
  summarise(mean = mean(posterior),
            lower = quantile(posterior, 0.05),
            upper = quantile(posterior, 0.95))
as.data.frame(ibrier_tab) %>% mutate_all(my_round)

comparisons <- contrast_models(
  int_brier,
  list_1 = rep("iClust2", 2),
  list_2 = c( "ER+/HER2-", "Reference"),
  seed = 2
)

compare <- ggplot(comparisons, size =  0.05) +
  theme_bw()
compare

diff_tab <- summary(comparisons, size = 0.05) %>%
  dplyr::select(contrast, starts_with("pract"))
diff_tab

ibrier_Tab <- iclust2prog::post_tab(diff_tab, ibrier_tab)
ibrier_Tab <- ibrier_Tab %>% mutate_all(my_round)
#ibrier_Tab[,c(1,3,2,4,5,6,7,8)]

require(stargazer)
stargazer(ibrier_Tab[,c(1,3,2,4,5,6,7)] , type = "latex", summary = FALSE, digits.extra = 3,
          digits = 3, digit.separator = ".")

library(gridExtra)
grid.arrange(pdf, compare, nrow = 1)



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



