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

### Prerequisites

This package is build in R to install R follow R-CRAN.

In addition, it's advisable to also install Stan, Installation and documentation can be found here: http://mc-stan.org/users/documentation/

### Installing
The package has to be installed in your R library either via git clone . or using devtools for example:
```
if(!require("devtools")) install.packages("devtools")
devtools::install_github("csetraynor/iclust2prog")
```

## Getting Started

The vignette explains how to perform MC-cross-validation in survival analysis , which is very appealing in Machine Learning (McLachlan, G., 2004) and has been proposed as an alternative to classical statistical analysis ( Benavoli et al, 2017)
The dataset we are going to work with contains curated clinical outcomes and genomic (gene expression and CNA) features from the METABRIC trial.
You will be able to download the dataset via the cBioPortal and this is made easy with the introduction of the cgdsr package.
For example:

```
library(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

the_study_list = getCancerStudies(mycgds)[25,1]
case_list = getCaseLists(mycgds, the_study_list)[2,1]
clinical_data <-  getClinicalData(mycgds, case_list)
```

Otherwise, the data which will be relevant for this study is also included in this package:
```
library(iclust2prog)

###Load data

data("intclustdat")
?intclustdat

```
The vignette Model_comparision serves as an approach to the modeling of high-dimensional genomic data in survival analysis. We can for example learn which model selection strategy has superior performance by comparing the integrated brier score. In addition, we can perform a post-hoc analysis and by HMC sample draws from the posterior distribution of BS. The figure bellow show the empirical distribution of 100-mc-cross-validation experiment for iclust2, and ER+/HER2-.

![posterior](https://user-images.githubusercontent.com/33321804/40587915-b5b389da-61cd-11e8-92b4-a8e1ad39a22f.png)

We can quantify the differences via the rope statistic for their practical statistical significance, e.g. whether the model is better or worse than the competitors and plot the results.

![contrasts](https://user-images.githubusercontent.com/33321804/40587912-aa0e0ef2-61cd-11e8-978c-c58a8717a334.png)

Finally, we can also make a simple plot of the gene pathways involved in our model by using clusterProfiler (Guangchuang Yu et al).

![enrichmap](https://user-images.githubusercontent.com/33321804/40587914-b21c17f6-61cd-11e8-9577-e6c995592ea1.png)

## Acknowledgment

I would like to thank my supervisors Prof Michael Chappell,
Dr Neil Evans, Dr Tarj Sahota and Ms Helen Tomkinson for giving me the opportunity to study for this PhD at Warwick. 

In addition, many thanks to the original authors of the study METABRIC (Pereira et al) and the creators of cBioPortal (Gao et al) for makin easier to share knowledge in biology and promote the development of science that may find cures for the difficult cancerous diseases.

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

McLachlan, G. Discriminant analysis and statistical pattern recognition. 544, (John Wiley & Sons, 2004).
