##### Relaxed Cox
library(iclust2prog)
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
relaxed_enet_fit <-   coxph(Surv(intclustdat$time, intclustdat$status)~ . ,
  data = X, init = iclust2_features$coef, control = coxph.control(iter.max = 5) )


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

# genedata <- data.frame(Hugo_Symbol = iclust2_features$feature,
#                        coef = iclust2_features$coef)

geneTable$Hugo_Symbol <- as.character(geneTable$Hugo_Symbol)
geneTable$Hugo_Symbol <- gsub("_cna.*", "", geneTable$Hugo_Symbol)
############ Create Gene List
# gene_names <- readr::read_tsv("brca_metabric/data_expression.txt")
#cna_names <- readr::read_tsv("brca_metabric/data_CNA.txt")
# cna_names <- data.frame(Hugo_Symbol = cna_names$Hugo_Symbol,
#                         Entrez_Gene_Id = cna_names$Entrez_Gene_Id)
# gene_names$Hugo_Symbol <-my_replace(gene_names$Hugo_Symbol)
# cna_names$Hugo_Symbol <-my_replace(cna_names$Hugo_Symbol)
# match(geneTable$Hugo_Symbol, cna_names$Hugo_Symbol)
# gene_names <- data.frame(Hugo_Symbol = gene_names$Hugo_Symbol,
#                         Entrez_Gene_Id = gene_names$Entrez_Gene_Id)
# target_names <- rbind(cna_names, gene_names)
# target_names <- target_names[complete.cases(target_names),]
# target_names <- unique(target_names)
# devtools::use_data(target_names, overwrite = T)

data("target_names")

geneTable$Entrez_Gene_Id <- target_names$Entrez_Gene_Id[match(geneTable$Hugo_Symbol, target_names$Hugo_Symbol)]
geneTable

#Check other Hugo Symbol synonims in databases
gene <-  geneTable$Entrez_Gene_Id
source("https://bioconductor.org/biocLite.R")
#biocLite(c("clusterProfiler","org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)

gene.df <- clusterProfiler::bitr(gene,
                                 fromType = "ENTREZID",
                                 toType = c("ENSEMBL", "SYMBOL"),
                                 OrgDb = org.Hs.eg.db)
gene.df

geneList <- geneTable[!is.na(geneTable$Entrez_Gene_Id),]$coef
names(geneList) <- geneTable[!is.na(geneTable$Entrez_Gene_Id),]$Entrez_Gene_Id
geneList

#devtools::install_github("csetraynor/predsurv")
library(predsurv)

ontology_search(ont = "c1", gene_list = geneTable)
ontology_search(ont = "c2", gene_list = geneTable)
ontology_search(ont = "c3", gene_list = geneTable)
ontology_search(ont = "c4", gene_list = geneTable)
ontology_search(ont = "c5", gene_list = geneTable)
ontology_search(ont = "c6", gene_list = geneTable)
ontology_search(ont = "H", gene_list = geneTable)















#Targets from Microarray
test <- data.frame( Hugo_Symbol = c("CMC4", "NMRK1", "MPA1B", "CSAG1", "NDUFA4L2", "NGF", "BOLA1",
          "CAMKK1", "GRP4", "DB451841", "AK055323"))
test$Entrez_Gene_Id <- target_names$Entrez_Gene_Id[match(test$Hugo_Symbol, target_names$Hugo_Symbol)]
test$coef <- seq(1:11)

testList <- test[!is.na(test$Entrez_Gene_Id),]$coef
names(testList) <- test[!is.na(test$Entrez_Gene_Id),]$Entrez_Gene_Id

###Enrichment analysis
library(DOSE)

gsecc <- gseGO(geneList=sort(geneList, decreasing = TRUE),
               minGSSize     =  1,
               maxGSSize =  length(geneList),
               pvalueCutoff = 0.1,
               OrgDb=org.Hs.eg.db,
               verbose=F)

head(summary(gsecc))

plotGOgraph(gsecc)
cnetplot(gsecc, foldChange=geneList)
enrichMap(gsecc, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
dotplot(gsecc, showCategory=50)

gsecc <- enrichGO(names(sort(geneList, decreasing = TRUE)),
                minGSSize     =  1,
                maxGSSize =  length(geneList),
                OrgDb = org.Hs.eg.db,
                pAdjustMethod = "DO")
head(summary(gsecc))

plotGOgraph(gsecc)
cnetplot(gsecc, foldChange=geneList)
enrichMap(gsecc, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
dotplot(gsecc, showCategory=30)

