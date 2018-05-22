##### Relaxed Cox


enet_data <- readRDS("C:/RFactory/Rdata_brca/enet_fit_data.RDS")

relaxed_enet_fit <- coxph(Surv(time, status) ~ npi+MTCP1NB+C9orf95+DB451841+MAP1B+CSAG1+NDUFA4L2+NGF+BOLA1+CAMKK1+AK055323+GPR4+age_std, data = enet_data)


CoxTab <- predsurv::par.table(relaxed_enet_fit)

genedata <- data.frame(Hugo_Symbol = CoxTab$Parameter,
                       coef = CoxTab$HR)
if("age_std" %in% genedata$Hugo_Symbol){
  geneTable <- genedata[(-match(c("age_std", "npi"),genedata$Hugo_Symbol)),]
}else{
  geneTable <- genedata
}


############ Create Gene List

devtools::install_github("GuangchuangYu/clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)
brcadata <- readr::read_tsv("C:/RFactory/Downloads/brca_metabric/data_expression.txt", col_names = TRUE)

gene_names <- data.frame(Entrez_Gene_Id = brcadata$Entrez_Gene_Id,
                         Hugo_Symbol = brcadata$Hugo_Symbol)

saveRDS(gene_names, "gene_names.RDS")

geneTable$Entrez_Gene_Id <- brcadata$Entrez_Gene_Id[match(geneTable$Hugo_Symbol, brcadata$Hugo_Symbol)]

geneTable
gene <-  geneTable$Entrez_Gene_Id

gene.df <- clusterProfiler::bitr(gene, fromType = "ENTREZID",
                                 toType = c("ENSEMBL", "SYMBOL"),
                                 OrgDb = org.Hs.eg.db)
gene.df$SYMBOL

geneList <- geneTable[!is.na(geneTable$Entrez_Gene_Id),]$coef
names(geneList) <- geneTable[!is.na(geneTable$Entrez_Gene_Id),]$Entrez_Gene_Id
geneList


###Enrichment analysis
library(DOSE)
ncg <- gseNCG(sort(geneList, decreasing = TRUE),
              nPerm         = 100,
              minGSSize = 1,
              maxGSSize = 9,
              pvalueCutoff  = 1,
              pAdjustMethod = "BH",
              verbose       = TRUE)

ncg <- setReadable(ncg, 'org.Hs.eg.db')
head(ncg, 3)

gsecc <- gseGO(geneList=sort(geneList, decreasing = TRUE),
               ont="BP",
               minGSSize = 1,
               maxGSSize = 9,
               OrgDb=org.Hs.eg.db,
               verbose=F)
head(summary(gsecc))

plotGOgraph(gsecc)
cnetplot(gsecc, foldChange=geneList)
enrichMap(gsecc, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
dotplot(gsecc, showCategory=30)
