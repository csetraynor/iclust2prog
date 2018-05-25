brca <- readRDS("C:/RFactory/Rdata_brca/brca_data.RDS")
cna <- readRDS("C:/RFactory/Rdata_brca/cna_expression.RDS")
brca <- cbind(brca, cna)


library(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

the_study_list = getCancerStudies(mycgds)[25,1]
case_list = getCaseLists(mycgds, the_study_list)[2,1]
clinical_data <-  getClinicalData(mycgds, case_list)

str(clinical_data)
colnames(clinical_data) <- tolower(colnames(clinical_data))
clinical_data <- tibble::rownames_to_column(clinical_data, var = "patient_id")
table(clinical_data$er_status)



er_status <- clinical_data[clinical_data$er_status == "+",]$patient_id
er_ihc <- clinical_data[clinical_data$er_ihc == "pos",]$patient_id
her2 <- clinical_data[clinical_data$her2_status == "+",]$patient_id
her2_3 <- clinical_data[clinical_data$threegene == "HER2+",]$patient_id

er_status <-gsub("\\.", "-", er_status)
er_ihc <-gsub("\\.", "-", er_ihc )
her2 <- gsub("\\.", "-", her2 )
her2_3 <- gsub("\\.", "-", her2_3 )

brca <- brca[brca$patient_id %in% er_status | brca$patient_id %in% er_ihc, ]
brca <- brca[!(brca$patient_id %in% her2),]
brca <- brca[!(brca$patient_id %in% her2_3),]
iclust2 <- brca[brca$intclust == 2,]


nrow( brca[!(brca$patient_id %in% her2_3),])

saveRDS(brca, "brca.RDS")
brca$os

unique(clinical_data$her2_status)
unique(clinical_data$threegene)
