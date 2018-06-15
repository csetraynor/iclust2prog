brca <- readRDS("/home/mtr/rfactory/brca_data.RDS")
cna <- readRDS("/home/mtr/rfactory/cna_expression.RDS")

brca <- cbind(brca, cna)

iclust2 <- brca[brca$intclust == 2, ]
saveRDS(iclust2, "iclust2_data.RDS")

#####

devtools::document()
data("clinical_data")
gene_expr <- readRDS("/home/mtr/rfactory/brca_metabric/brca_data.RDS")
gene_cna <- readRDS("/home/mtr/rfactory/brca_metabric/cna_expression.RDS")
combined <- cbind(gene_expr, gene_cna)

library(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

the_study_list = getCancerStudies(mycgds)[25, 1]
case_list = getCaseLists(mycgds, the_study_list)[2, 1]
clinical_data <- getClinicalData(mycgds, case_list)
colnames(clinical_data) <- tolower(colnames(clinical_data))
clinical_data <- tibble::rownames_to_column(clinical_data, var = "patient_id")
dplyr::glimpse(clinical_data)
treatment_effect <- clinical_data %>%
  select(patient_id, chemotherapy, radio_therapy)
devtools::use_data(treatment_effect)


str(clinical_data)
dplyr::glimpse(clinical_data)
colnames(clinical_data) <- tolower(colnames(clinical_data))
clinical_data <- tibble::rownames_to_column(clinical_data, var = "patient_id")
table(clinical_data$er_status)

er_status <- clinical_data[clinical_data$er_status == "+",]$patient_id
her2 <- clinical_data[clinical_data$her2_status == "-",]$patient_id
three_gene <- clinical_data[clinical_data$threegene == "ER+/HER2- High Prolif" | clinical_data$threegene == "ER+/HER2- Low Prolif",]$patient_id


er_status <-gsub("\\.", "-", er_status)
her2 <- gsub("\\.", "-", her2 )
three_gene <- gsub("\\.", "-", three_gene )

er <- clinical_data[ (clinical_data$patient_id %in% er_status & clinical_data$patient_id %in% her2) |  clinical_data$patient_id %in% three_gene,]
iclust2 <- clinical_data[clinical_data$intclust == 2,]


combined <- combined[ (combined$patient_id %in% er_status & combined$patient_id %in% her2) |  combined$patient_id %in% three_gene,]
iclust2 <- combined[combined$intclust == 2,]


saveRDS(combined, "combined.RDS")

combined <- readRDS("combined.RDS")

gene_expr <- combined[, c(colnames(gene_expr))]
gene_cna <- combined[, c(colnames(gene_cna))]

gene_expr$patient_id <- NULL
gene_expr$intclust <- NULL
gene_expr$npi <- NULL
gene_expr$os_months <- NULL
gene_expr$os_deceased <- NULL


library(foreach)
library(parallel)

require(doMC)
registerDoMC(cores=4)

correlations =  lapply(seq_along(gene_expr), function(i) {
    correlations_i = lapply(seq_along(gene_cna),
                              function(x)  {
      ltm::biserial.cor(x = as.double(gene_expr[,i]),
                        y = as.factor(gene_cna[,x]))
    })
    as.data.frame(do.call(cbind, correlations_i) )
})
as.data.frame(do.call(rbind, correlations) )
saveRDS(correlations, "correlations.RDS")



# correlations
row_indic <- apply(correlations, 1, function(x) sum(x > 0.3 | x < -0.3) > 1)
correlations<- correlations[row_indic ,row_indic ]
pdf('corrplot.pdf')
corrplot::corrplot(correlations, method="square")
dev.off()
