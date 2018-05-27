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
The vignette Model_comparision serves as an approach to the modeling high-dimensional genomic data in survival analysis. We can for example learn which model has greater performance by comparing the integrated brier score. The figure belllow is the empirical distribution of 100-mc-cross-validation experiment for iclust2, erpos, and pooled analyses.

<p>
    <img src="https://github.com/csetraynor/iclust2prog/Plots/IBrier_pdf.pdf" width="220" height="240" />
</p>

In addition, we can perform a post-hoc analysis and by HMC sample draws from the posterior of the IBRIER distributions. The figure bellow shos the results:

[IBRIER_POSTERIOR.pdf](https://github.com/csetraynor/iclust2prog/Plots/IBRIER_POSTERIOR.pdf)

Finally, we are going to consider if the differences are significant, whether the model is better or worse than the competitors. We can clarify the differences by the plot bellow.

[diff_plot.pdf](https://github.com/csetraynor/iclust2prog/Plots/diff_plot.pdf)

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

McLachlan, G. Discriminant analysis and statistical pattern recognition. 544, (John Wiley & Sons, 2004).
