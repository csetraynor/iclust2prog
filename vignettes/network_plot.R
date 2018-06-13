devtools::document()
data("int_coeff_tab")

library(ggnet)
library(network)
library(sna)

#example
# random graph
net = rgraph(10, mode = "graph", tprob = 0.5)
net = network(net, directed = FALSE)

# vertex names
network.vertex.names(net) = letters[1:10]
ggnet2(net)

oncosig <- ontology_search(ont = "c6", gene_list = geneTable)
names(oncosig)
kras <- oncosig[grep("KRAS", names(oncosig))]
entrez <- unique(unlist(kras))

names(kras) <- tolower(names(kras))
kras_genes <- c(rep("TBK1", 3), rep("KRAS", 7))


hnet = network(hmat, directed = FALSE)
ggnet2(hnet)$data

ggnet2(hnet, label = TRUE)


oncomodule <- ontology_search(ont = "c6", gene_list = geneTable)
module <- names(oncomodule)
entrez <- unique(unlist(oncomodule))
oncomat <- get_mat(oncomodule)
onconet = network(oncomat, directed = FALSE)
ggnet2(onconet)$data
ggnet2(onconet,  color = "color")

onconet %v% "regulate" = ifelse(grepl("DN", names(oncomodule)), "down", "up")
onconet %v% "regulation" = ifelse(onconet %v% "regulate" == "up", "red", "green")



hrsize = sapply(oncomodule, function(x){
  round(sum(geneTable$HR[match(geneTable$Entrez_Gene_Id, unlist(x) )],
      na.rm = T), 3)
})

onconet %v% "Hazard ratio" = hrsize

ggnet2(onconet,  color = "regulation", size = "Hazard ratio",
       label = T)

