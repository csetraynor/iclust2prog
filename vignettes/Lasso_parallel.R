#mc-cross-validation.R

#Load libraries
library(dplyr)
require(doMC)

#Load data
ERpos = readRDS("ERpos_data.RDS")
iclust2 = readRDS("iclust2_data.RDS")

#Get brca iclust2

iclust2$intclust <- NULL
iclust2$patient_id <- NULL
ERpos$intclust <- NULL
ERpos$patient_id <- NULL


model_fit <- function(dat){

# create predictor matrix
x <- dat %>% dplyr::select(-os_months, -os_deceased)

###### Fit models
##################################################
#grouped enet
p.fac = rep(1, ncol(x))
p.fac[match(c("npi","age_std"), colnames(x))] = 0
#prepare
x <- as.matrix(x)
y <- as.matrix(dat %>%
                 dplyr::select(time = os_months,
                               status = os_deceased), ncol = 2)
require(doMC)
registerDoMC(cores=8)
### Apply elastic net with a=1 to perform Lasso
mod <-  glmnet::cv.glmnet(x, y, family = "cox",
                          grouped = TRUE,
                          alpha = 1,
                          parallel = TRUE, penalty.factor = p.fac)
return(mod)
}

mod <- model_fit(iclust2)

saveRDS(mod, "iclust2_lasso.RDS")
rm(mod)
rm(iclust2)
modE <- model_fit(ERpos)
saveRDS(modE, "erpos_lasso.RDS")


