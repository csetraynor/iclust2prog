devtools::document()
# library(iclust2prog)
library(glmnet)
library(purrr)
library(dplyr)
data("ic2dat")

set.seed(9666)
mc_samp <- bootstraps(ic2dat, strata = "status", times = 100)

iclust2_features <- extract_features(iclust2_glmnet)
iclust2_features$feature <- my_replace(iclust2_features$feature)
colnames(ic2dat) <- my_replace(colnames(ic2dat))

# Test with npi, and age -----------------------
ic2dat <- ic2dat %>% dplyr::select(time, status, npi, age_std)

mc_samp$mod_iclust2 <- pmap(list(mc_samp$splits), function(data) {
  mod_fit(x = data, form = iclust2_features, iter = 1)
})

mc_samp$CoxTab_iclust2 <- pmap(list(mc_samp$mod_iclust2), function(mod) {
  par.table(mod, log_scale = TRUE)
})

mc_samp$tab_iclust2 <- pmap(list(mc_samp$CoxTab_iclust2), function(tab) {
  get_table(tab)
})

coeff_tab <- get_coeff(mc_samp, tab = "tab_iclust2")

gene_mod <- coeff_Tab %>% dplyr::mutate(kras = rowSums(.[c("NGF", "GPR4")])) %>%
  dplyr::rename(lef1 = MAP1B) %>% dplyr::select(splits, id, kras, lef1)


int_gene <- tidyposterior::perf_mod(gene_mod, seed = 6507, iter = 5000)

int_gene_tab <- tidy(int_gene) %>% dplyr::group_by(model) %>% dplyr::summarise(HR = exp(mean(posterior)),
                                                                               lower = exp(quantile(posterior, 0.05)), upper = exp(quantile(posterior,
                                                                                                                                            0.95)))
as.data.frame(int_gene_tab) %>% dplyr::mutate_all(my_round)

# Going Bayesian -----------

library("rstanarm")
library("bayesplot")
library("ggplot2")
library("broom")

setwd("~/work/RSTANARM/Tutorial")
data("sampdata")
dat <- read.csv("sampledata.csv",
                header=TRUE)

library(rethinking)
m1.1 <- rethinking::map(
  alist(
    hr ~ dnorm( mu, sigma) ,
    mu ~ dnorm( 0, 10) ,
    sigma ~ dunif( 0, 50)
  ), data = as.data.frame(coeff_tab %>%
                            dplyr::mutate(hr = npi))
)

precis(m1.1)

m1.2 <- rethinking::map(
  alist(
    hr ~ dnorm( mu, sigma) ,
    mu ~ dnorm( 0, 1) ,
    sigma ~ dunif( 0, 50)
  ), data = as.data.frame(coeff_tab %>%
                            dplyr::mutate(hr = npi))
)

precis(m1.2, prob = 0.95)


# Poisson Model ----------------------------

#prepare long format dataset
long_ic2dat <- gen_stan_dat(ic2dat)

m1stan <- rethinking::map2stan(
  alist(
    status ~ dpois(hazard),
    log(hazard) <- log_baseline_hazard[t_id] + log_dtime,
    log_baseline_hazard[t_id] ~ dnorm(0, 1)
  ), data = list( t_id = as.integer(unlist(long_ic2dat$t_id ) ),
                  log_dtime = as.double(unlist(long_ic2dat$log_dtime ) ),
                  status = as.integer(unlist(long_ic2dat$status ) )
  )
)

m2stan <- rethinking::map2stan(
  alist(
    status ~ dpois(hazard),
    log(hazard) <- log_baseline_hazard_mu[t_id] + log_dtime + log_baseline_hazard_raw[t_id],
    log_baseline_hazard_mu[t_id] ~ dnorm(0, 1),
    log_baseline_hazard_raw[t_id] ~ dnorm(0, sigma_hazard),
    sigma_hazard ~ dcauchy(0, 2)
  ), data = list( t_id = as.integer(unlist(long_ic2dat$t_id ) ),
                  log_dtime = as.double(unlist(long_ic2dat$log_dtime ) ),
                  status = as.integer(unlist(long_ic2dat$status ) )
  ) , chains = 2 , iter = 4000 , warmup = 1000
)
show(m2stan)
plot(m2stan)
#devtools::install_github("rmcelreath/rethinking", force = TRUE) workaround for plots error
post <- extract.samples(m2stan)
log_baseline_hazard <- post$log_baseline_hazard_mu + post$log_baseline_hazard_raw
str(post)
dens(post$log_baseline_hazard[,1])

times <- ic2dat %>%
  filter(status) %>% select(time) %>% unique %>% ungroup %>% arrange(time) %>% unlist %>% as.double

exp.link <- function(x) {
  exp( x)
}
haz <- exp.link(log_baseline_hazard)
cumhaz <- apply( haz, 1,  cumsum )
surv <- as.data.frame( exp(- cumhaz ) )
surv <- lapply(seq_along(surv), function(x) {
  data.frame( surv = surv[[x]],
              sample = x,
              times = times)
}
)
surv <- do.call(rbind, surv)
p <- plot_km(ic2dat)
p + geom_line(data = surv, aes(x = times, y = surv, group = sample))

# ------ Generalised Additive Models
m1.gam <- gamm4::gamm4(status~1+offset(log_dtime)+s(time),
                 data = long_ic2dat , family='poisson')

#---- The ML approach
print(m1.gam)

time1 <- long_ic2dat %>% filter(chemotherapy == "YES") %>%
  dplyr::select(time, patient_id) %>%
  arrange(time) %>%
  mutate( log_dtime1 =  log( diff(c(0, time) ) ) )


new_dat <- data.frame(log_dtime = log_dtime, time = time)
ic2dat$log_dtime <- log(ic2dat$time)
p1 <- predict(mgcv::gam(test ~ 1 + offset(log_dtime) + s(time) + factor(chemotherapy), long_ic2dat, family='poisson'), new_dat)
p1
S1<-exp(-cumsum(exp(p1)))

p <- plot_km(ic2surv, strata = "chemotherapy")
p <- p + geom_line(data = data.frame(times = time,
                                     surv = S1,
                                     chemo = ic2surv$chemotherapy), aes(x = times, y = surv, group = chemo))
p

#---------- With bayes csetraynor/tidybayes

#---- The Bayesian model

library(magrittr)
library(dplyr)
library(forcats)
library(tidyr)
library(purrr)
library(modelr)
library(tidybayes)
library(ggplot2)
library(ggstance)
library(ggridges)
library(rstan)
library(rstanarm)
import::from(LaplacesDemon, invlogit)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

m1.stan_gam <- rstanarm::stan_gamm4(status~1+offset(log_dtime)+s(time)+factor(chemotherapy), data = long_ic2dat , family='poisson')
y_rep <- posterior_predict(m1.stan_gam)
dim(y_rep)
summary(m1.stan_gam)

parameters(m1.stan_gam)
post <- as.matrix(m1.stan_gam)
plot_nonlinear(m1.stan_gam)

long_ic2dat$haz <- predict(m1.stan_gam, newdata = long_ic2dat)

plot.frame <- long_ic2dat %>%
  group_by(sample_id) %>%
  mutate(surv = exp( - cumsum(exp(haz))) ) %>%
  filter(row_number()==n()) %>%
  ungroup() %>%
  arrange(time)

p <- plot_km(ic2dat)
p + geom_line(data = plot.frame, aes(x = time, y = surv), linetype = 2)


# Calculta Brier Score ----------------
set.seed(10)
newdata <- ic2dat[sample(nrow(ic2dat), 15 ), ] %>% arrange(time)
timepoints <-  seq(0, max(newdata$time),
                   length.out = 100L)
longdatanew <- gen_stan_dat(newdata, timepoints = timepoints)
longdatanew$testml <- predict(m1.stan_gam, newdata = longdatanew)
test <- longdatanew %>% select(time, testml, sample_id)
test_matrix <- split(test, as.factor(test$sample_id)) %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="time"), .) %>%
  select(contains("testml")) %>%
  as.matrix() %>%
  exp

cumhaz <- apply( test_matrix, 1,  cumsum )

probs <- exp( - cumhaz )

brier <- pec::pec(probs, Surv(time, status) ~ 1,
                  data = newdata,
                  maxtime = max(timepoints),
                  exact = FALSE,
                  exactness = nrow(test2) - 1)

#--- Full Bayesian
post <- posterior_predict(m1.stan_gam)
long_ic2dat$test <- post[1, ]
