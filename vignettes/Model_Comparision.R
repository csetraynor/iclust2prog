if(!require("devtools")) install.packages("devtools")
if(!require("iclust2prog")) devtools::install_github("iclust2prog")

library(glmnet)
library(purrr)
library(dplyr)
library(tidyr)
library(rsample)
library(tidyposterior)
library(gridExtra)
library(bindrcpp)
devtools::document()
theme_set(theme_bw())

###Load data
data("eph2n_glmnet")
data("iclust2_glmnet")

#Plot glmnet fits
glmnet::plot.cv.glmnet(eph2n_glmnet)
glmnet::plot.cv.glmnet(iclust2_glmnet)

#Extract features, see my_replace in utils
eph2n_features <- extract_features(eph2n_glmnet)
eph2n_features$feature <-my_replace(eph2n_features$feature)

iclust2_features <- extract_features(iclust2_glmnet)
iclust2_features$feature <-my_replace(iclust2_features$feature)

############ Survival analysis vignette

data("ic2dat")

set.seed(9666)
mc_samp <- mc_cv(ic2dat, strata = "status", times = 100)

cens_rate <- function(x) mean(analysis(x)$status == 1)
summary(map_dbl(mc_samp$splits, cens_rate))

############### Create models
mc_samp$mod_iclust2 <- pmap(list(mc_samp$splits),
                            function(data){
                            mod_fit(x = data,
                                      form = iclust2_features,
                                    iter = 1)
                            })
mc_samp$mod_eph2n <- pmap(list(mc_samp$splits),
                          function(data){
                            mod_fit(x = data,
                                    form = eph2n_features,
                                    iter = 1)
                          })


############### Get Brier
mc_samp$brier_eph2n <- pmap(list(mc_samp$splits, mc_samp$mod_eph2n),
                            function(data, model){
                              get_tdbrier(data = data,
                                          mod = model,
                                          form = eph2n_features
                              )
                            })
mc_samp$brier_iclust2 <- pmap(list(mc_samp$splits, mc_samp$mod_iclust2),
                              function(data, model){
                                get_tdbrier(data = data,
                                            mod = model,
                                            form = iclust2_features)
                              })

###integrate Brier
mc_samp$'iC-2' <- map_dbl(mc_samp$brier_iclust2, integrate_tdbrier)
mc_samp$'ER+/HER2-' <- map_dbl(mc_samp$brier_eph2n, integrate_tdbrier)
mc_samp$Reference <- map_dbl(mc_samp$brier_eph2n, integrate_tdbrier_reference)

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

pdf <- ggplot(tidy(int_brier)%>%
                filter(model != "Reference")) +
  theme_bw() +
  ylab("Posterior probability of BS")
pdf
ibrier_tab <- tidy(int_brier) %>%
  filter(model != "Reference") %>%
  group_by(model) %>%
  summarise(mean = mean(posterior),
            lower = quantile(posterior, 0.05),
            upper = quantile(posterior, 0.95))
as.data.frame(ibrier_tab) %>% mutate_all(my_round)

comparisons <- contrast_models(
  int_brier,
  list_1 = rep("iC-2", 1),
  list_2 = c( "ER+/HER2-"),
  seed = 20
)

compare <- ggplot(comparisons, size =  0.01) +
  theme_bw()
compare <- compare + ylab("") + xlab("DeltaBS")

diff_tab <- summary(comparisons, size = 0.01) %>%
  dplyr::select(contrast, starts_with("pract"))
diff_tab

ibrier_Tab <- post_tab(diff_tab, ibrier_tab)
ibrier_Tab <- ibrier_Tab %>% mutate_all(my_round)
pdf("dens1.pdf", 7 ,5)
grid.arrange(pdf, compare, nrow = 1)
dev.off()
