##### Relaxed Cox
library(iclust2prog)
library(glmnet)
data("intclustdat")
intclustdat <- intclustdat %>%
  dplyr::rename(time = os_months,
         status = os_deceased) %>%
  dplyr::mutate(status = status == 1)
data("iclust2_glmnet")

iclust2_features <- extract_features(iclust2_glmnet)
iclust2_features$feature <-my_replace(iclust2_features$feature)
colnames(intclustdat) <- my_replace(colnames(intclustdat))

X <- intclustdat[,iclust2_features$feature]
relaxed_enet_fit <- suppressWarnings(coxph(Surv(intclustdat$time, intclustdat$status)~ . ,
  data = X, init = iclust2_features$coef, control = coxph.control(iter.max = 5) ))


CoxTab <- par.table(relaxed_enet_fit)
require(stargazer)
stargazer(CoxTab  , type = "latex", summary = FALSE, digits.extra = 3,
          digits = 3, digit.separator = ".",
          rownames = F)


genedata <- data.frame(Hugo_Symbol = CoxTab$Parameter,
                       coef = CoxTab$HR)
if("age_std" %in% genedata$Hugo_Symbol){
  geneTable <- genedata[(-match(c("age_std", "npi"),genedata$Hugo_Symbol)),]
}else{
  geneTable <- genedata
}

geneTable$Hugo_Symbol_cna <- as.character(geneTable$Hugo_Symbol)
geneTable$Hugo_Symbol <- gsub("_cna.*", "", geneTable$Hugo_Symbol_cna)

####Load data target names
data("target_names")

geneTable$Entrez_Gene_Id <- target_names$Entrez_Gene_Id[match(geneTable$Hugo_Symbol, target_names$Hugo_Symbol)]
geneTable

#Optional Check other Hugo Symbol synonims in databases
# gene <-  geneTable$Entrez_Gene_Id
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("clusterProfiler","org.Hs.eg.db")
# library(clusterProfiler)
# library(org.Hs.eg.db)
#
# gene.df <- clusterProfiler::bitr(gene,
#                                  fromType = "ENTREZID",
#                                  toType = c("ENSEMBL", "SYMBOL"),
#                                  OrgDb = org.Hs.eg.db)
# gene.df

geneList <- geneTable[!is.na(geneTable$Entrez_Gene_Id),]$coef
names(geneList) <- geneTable[!is.na(geneTable$Entrez_Gene_Id),]$Entrez_Gene_Id
geneList

ontology_search(ont = "c1", gene_list = geneTable)
ontology_search(ont = "c2", gene_list = geneTable)
ontology_search(ont = "c3", gene_list = geneTable)
ontology_search(ont = "c4", gene_list = geneTable)
ontology_search(ont = "c5", gene_list = geneTable)
ontology_search(ont = "c6", gene_list = geneTable)
ontology_search(ont = "H", gene_list = geneTable)

