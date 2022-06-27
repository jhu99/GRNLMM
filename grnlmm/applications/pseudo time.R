#calculate pseudo time information
library(monocle3)
x<-read.csv("tumor_全exp.csv",header = TRUE, sep = ",", quote = "\"",
            dec = ".", fill = TRUE, comment.char = "!")
row.names(x)<-x[,1]
x<-x[,-c(1)]
g<-g[,-c(1)]
g<-as.matrix(g)
g<-g[1:100]
g<-as.matrix(row.names(x))
row.names(g)<-g[,1]
colnames(g)<-"gene_short_name"
c<-matrix(0,1568,1)
row.names(c)<-colnames(x)
colnames(c)<-"cell_type"
c[,1]<-"B cells"
x<-as.matrix(x)
cds <- new_cell_data_set(x, cell_metadata = c, gene_metadata = g)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds<-order_cells(cds,root_cells = "ATCATCTGTAGCGCTC_DLBCL002B")
plot_cells(cds,color_cells_by = "pseudotime")
write.table(time,col.names = F,row.names = T,file = "time.txt",sep="\t")