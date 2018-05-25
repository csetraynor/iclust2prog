##### Relaxed Cox
library(iclust2prog)
data("intclustdat")
data("iclust2_glmnet")

iclust2_features <- extract_features(iclust2_glmnet)
iclust2_features$feature[match("`GCSAML-AS1_cna`1", iclust2_features$feature)] <- "GCSAML_AS1_cna1"
colnames(intclustdat)[match("GCSAML-AS1_cna1", colnames(intclustdat))] <- "GCSAML_AS1_cna1"

X <- intclustdat[,iclust2_features$feature]
relaxed_enet_fit <-   coxph(Surv(intclustdat$time, intclustdat$status)~ . ,
  data = X, init = iclust2_features$coef, control = coxph.control(iter.max = 1) )


CoxTab <- par.table(relaxed_enet_fit)

genedata <- data.frame(Hugo_Symbol = CoxTab$Parameter,
                       coef = CoxTab$HR)
if("age_std" %in% genedata$Hugo_Symbol){
  geneTable <- genedata[(-match(c("age_std", "npi"),genedata$Hugo_Symbol)),]
}else{
  geneTable <- genedata
}
geneTable$Hugo_Symbol <- as.character(geneTable$Hugo_Symbol)

############ Create Gene List
data("gene_names")
gene_names <- readr::read_tsv("C:/RFactory/bymetabric_files/metabricdata/brca_metabric/data_expression.txt")
cna_names <- readr::read_tsv("C:/RFactory/bymetabric_files/metabricdata/brca_metabric/data_CNA.txt")
cna_names <- data.frame(Hugo_Symbol = cna_names$Hugo_Symbol,
                        Entrez_Gene_Id = cna_names$Entrez_Gene_Id)
gene_names <- data.frame(Hugo_Symbol = gene_names$Hugo_Symbol,
                        Entrez_Gene_Id = gene_names$Entrez_Gene_Id)

target_names <- rbind(cna_names, gene_names)
target_names <- unique(target_names)

target_names$Hugo_Symbol[duplicated(target_names$Hugo_Symbol) == T]

target_names$Entrez_Gene_Id <- gene_names$Entrez_Gene_Id[match(target_names$Hugo_Symbol, gene_names$Hugo_Symbol)]

target_names$Entrez_Gene_Id <- cna_names$Entrez_Gene_Id[match(target_names$Hugo_Symbol, cna_names$Hugo_Symbol)]

geneTable$Hugo_Symbol <- gsub("_cna*", "", geneTable$Hugo_Symbol)

geneTable$Entrez_Gene_Id <- target_names$Entrez_Gene_Id[match(geneTable$Hugo_Symbol, target_names$Hugo_Symbol)]

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

gsecc <- gseDO(geneList=sort(geneList, decreasing = TRUE),
               minGSSize = 1,
               maxGSSize = 9,
               verbose=F)
head(summary(gsecc))

plotGOgraph(gsecc)
cnetplot(gsecc, foldChange=geneList)
enrichMap(gsecc, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
dotplot(gsecc, showCategory=30)


