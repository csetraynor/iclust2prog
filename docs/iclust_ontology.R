C1 <- qusage::read.gmt("C:/RFactory/predsurv/vignettes/Gene sets/c1.all.v6.1.entrez.gmt")

C7 <- qusage::read.gmt("C:/RFactory/clusterProfiler/vignettes/Gene sets/c7.all.v6.1.entrez.gmt")

gene_names <- readRDS("C:/RFactory/Rdata_brca/gene_names.RDS")

match("GATA3", gene_names$Hugo_Symbol)

gene_names$Entrez_Gene_Id[match("GATA3", gene_names$Hugo_Symbol)]

mg <- function(x){
  gene_names$Entrez_Gene_Id[match(x, gene_names$Hugo_Symbol)]
}


intClust1 <- unique(c(C1$chr17q23, mg("GATA3")  , mg("RPS6KB1") , mg("PPM1D,"), mg("PTRH2"), mg("APPBP2") ) )

intClust1 <-intClust1[!is.na(intClust1)]

intClust2 <- unique(c( C1$chr11q13, C1$chr11q14, mg("CCND1")  ,mg("EMSY")  ,mg("PAK1") ) )
intClust2 <-intClust2[!is.na(intClust2)]

intClust3 <- unique(c(mg("TP53")  ,mg("PIK3CA")  ,mg("CDH1"), mg("RUNX1") ) )
intClust3 <-intClust3[!is.na(intClust3)]

intClust4 <- unique(c( C7[grep("TCR", names(C7))]  ) )
intClust4 <-intClust4[!is.na(intClust4)]

intClust5 <- unique(c(C1$chr17q12, mg("TP53") ))
intClust5 <-intClust5[!is.na(intClust5)]

intClust6 <- unique(c(C1$chr8p12 , mg("ZNF703")))
intClust6 <-intClust6[!is.na(intClust6)]

intClust7 <- unique(c(C1$chr16p, C1$chr8q, mg("MAP3K1"), mg("CTCF")))
intClust7 <-intClust7[!is.na(intClust7)]

intClust8 <- unique(c(C1$chr1q, mg("PIK3CA"), mg("GATA3"), mg("MAP2K4") ) )
intClust8 <-intClust8[!is.na(intClust8)]

intClust9 <- unique(c( C1$chr8q, mg("TP53") , mg("PPP2R2A")   ))
intClust9 <- intClust9[!is.na(intClust9)]

intClust10 <- unique(c( mg("TP53") , C1$chr5q, C1$chr8q, C1$chr12p, 
                        mg("AURKB"),  mg("BCL2"),mg("BUB1"),mg("CDCA3"),mg("CDCA4"),mg("CDC20"),mg("CDC45"),mg("CHEK1"),mg("FOXM1"),mg("HDAC2"), mg("IGF1R"),mg("KIF2C"),mg("KIFC1"),mg("MTHFD1L"),mg("RAD51AP1"),mg("TTK"), mg("UBE2C") ))
intClust10 <- intClust10[!is.na(intClust10)]


intClustGO <- list(intClust1,intClust2,
                   intClust3,intClust4,
                   intClust5,intClust6,
                   intClust7,intClust8,
                   intClust9,intClust10)


