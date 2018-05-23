# Package iclust2prog
iclust2prog



This package includes relevant code and figures of various prognostic models in iClust2-patients from the METABRIC trial.
A recent breakthrough in breast cancer has stratified the METABRIC cohort (Pereira et al) using integrative clustering confirming
that there are at least 10 subtypes of the disease which have different clinical outcome. Depending on
the category patients may be spare unnecessary chemotherapy or benefit of targeted antibody HER-2 therapy
(intClust-3 and intClust-5, respectively). Moreover, it is established that intClust-2 (currently into ER-positive)
has complicated medical management and prognosis. In this study, we revisit the (Pereira et al) dataset, downloaded via the
cBioPortal platform  (Gao et al), to describe a novel gene signature predictive of prognosis in intClust-2 patients, with the
potential to develop new tailored therapies by mapping genes to pathological pathways.


## Getting Started

This vignette serves as an approach to the modeling of survival data in medical statistics.
```
library(iclust2prog)
library(dplyr)
library(tidyr)
library(tidyposterior)
library(DOSE)
theme_set(theme_bw())

###Load data

data("intclustdat")


intclustdat <- intclustdat %>%
  rename(time = os_months,
         status = os_deceased,
         "Sep_1" = "1-Sep") %>%
  mutate(status = status == 1)


mc_samp <- mc_cv(intclustdat, strata = "status", times = 100)

mc_samp$mod_pooled <- pmap(list(mc_samp$splits),
                            function(data){
                              mod_fit(x = data,
                                      form = pooled,
                                      inits = pooled_inits,
                                      iter = 0)
                            })
mc_samp$brier_pooled <- pmap(list(mc_samp$splits, mc_samp$mod_pooled),
                              function(data, model){
                                get_tdbrier(data = data,
                                            mod = model,
                                            inits = pooled_inits,
                                            iters = 0
                                            )
                              })

```

### Prerequisites

This package is build in R to install R follow R-CRAN.

In addition, it's advisable to also install Stan, Installation and documentation can be found here: http://mc-stan.org/users/documentation/


### Installing
```
devtools::install_github("csetraynor/iclust2prog")
```

## Acknowledgment

I would like to thank my supervisors Prof Michael Chappell,
Dr Neil Evans, Dr Tarj Sahota and Ms Helen Tomkinson for giving me the opportunity to study for this PhD at Warwick. 

## References

  Pereira, B. et al. Nat. Commun. 7, 2016.
  
  Gao et al. Sci. Signal. 2013 & Cerami et al. Cancer Discov. 2012

Noah Simon, Jerome Friedman, Trevor Hastie, Rob Tibshirani (2011).
  Regularization Paths for Cox's Proportional Hazards Model via
  Coordinate Descent. Journal of Statistical Software, 39(5), 1-13. URL
  http://www.jstatsoft.org/v39/i05/.

Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012). Evaluating
  Random Forests for Survival Analysis Using Prediction Error Curves.
  Journal of Statistical Software, 50(11), 1-23. URL
  http://www.jstatsoft.org/v50/i11/.

Max Kuhn and Hadley Wickham (2017). rsample: General Resampling
  Infrastructure. R package version 0.0.2.
  https://CRAN.R-project.org/package=rsample
  
Liang Li, Cai Wu Department of Biostatistics and The University of Texas MD
  Anderson Cancer Center (2016). tdROC: Nonparametric Estimation of
  Time-Dependent ROC Curve from Right Censored Survival Data. R package version
  1.0. https://CRAN.R-project.org/package=tdROC

Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He. clusterProfiler: an R
  package for comparing biological themes among gene clusters. OMICS: A Journal
  of Integrative Biology 2012, 16(5):284-287



