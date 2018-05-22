##### Relaxed Cox
data("intclustdat")

relaxed_enet_fit <- coxph(Surv(time, status) ~ npi+MTCP1NB+C9orf95+DB451841+MAP1B+CSAG1+NDUFA4L2+NGF+BOLA1+CAMKK1+AK055323+GPR4+age_std, data = enet_data)


CoxTab <- par.table(relaxed_enet_fit)

genedata <- data.frame(Hugo_Symbol = CoxTab$Parameter,
                       coef = CoxTab$HR)
if("age_std" %in% genedata$Hugo_Symbol){
  geneTable <- genedata[(-match(c("age_std", "npi"),genedata$Hugo_Symbol)),]
}else{
  geneTable <- genedata
}

############ Create Gene List
data("gene_names")
geneTable$Entrez_Gene_Id <- brcadata$Entrez_Gene_Id[match(geneTable$Hugo_Symbol, brcadata$Hugo_Symbol)]

geneTable
gene <-  geneTable$Entrez_Gene_Id

devtools::install_github("GuangchuangYu/clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)

gene.df <- clusterProfiler::bitr(gene, fromType = "ENTREZID",
                                 toType = c("ENSEMBL", "SYMBOL"),
                                 OrgDb = org.Hs.eg.db)
gene.df$SYMBOL

geneList <- geneTable[!is.na(geneTable$Entrez_Gene_Id),]$coef
names(geneList) <- geneTable[!is.na(geneTable$Entrez_Gene_Id),]$Entrez_Gene_Id
geneList

###Enrichment analysis
library(DOSE)

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
