#calculate the degree of graph
library(igraph)
graph<-graph_from_data_frame(reference,directed = T)
de<-degree(graph,mode="total")
de<-as.matrix(sort(de,decreasing = T))
de<-cbind(de,de[,1])
de[,2]<-row.names(de)
#select the top 595 genes
ss<-x[c(de[1:595,2]),]