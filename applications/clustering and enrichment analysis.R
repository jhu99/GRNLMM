#k-means
netw<-read.csv("network1-5.csv",header = TRUE, sep = ",", quote = "\"",
               dec = ".", fill = TRUE, comment.char = "!")
cluster1<-kmeans(x,centers = 5,nstart = 20)
rescluster1<-as.matrix(cluster1[["cluster"]])
rescluster1<-rescluster1[order(rescluster1[,1],decreasing = F),]
rescluster1<-as.matrix(rescluster1)
class1<-row.names(rescluster1)[1:15]
class2<-row.names(rescluster1)[36:64]
class3<-row.names(rescluster1)[85:100]
class4<-row.names(rescluster1)[65:84]
class5<-row.names(rescluster1)[16:35]

#enrichment analysis
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
output1<-bitr(class1,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
output2<-bitr(class2,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
output3<-bitr(class3,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
output4<-bitr(class4,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
output5<-bitr(class5,fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = 'org.Hs.eg.db')
ego1<-enrichGO(gene = output1[,2],OrgDb = org.Hs.eg.db,keyType = 'ENSEMBL',ont = "BP",readable = T)
ego2<-enrichGO(gene = output2[,2],OrgDb = org.Hs.eg.db,keyType = 'ENSEMBL',ont = "BP",readable = T)
ego3<-enrichGO(gene = output3[,2],OrgDb = org.Hs.eg.db,keyType = 'ENSEMBL',ont = "BP",readable = T)
ego4<-enrichGO(gene = output4[,2],OrgDb = org.Hs.eg.db,keyType = 'ENSEMBL',ont = "BP",readable = T)
ego5<-enrichGO(gene = output5[,2],OrgDb = org.Hs.eg.db,keyType = 'ENSEMBL',ont = "BP",readable = T)
#drawing
barplot(ego1,showCategory = 15)
C1<-matrix(0,15,2)
C1[,1]<-ego1@result[["Description"]][1:15]
C1[,2]<--log(ego1@result[["p.adjust"]][1:15])
C1<-as.data.frame(C1)
C1[,2]<-as.numeric(C1[,2])
colnames(C1)<-c("fun","p")
C1<-within(C1,{fun<-factor(fun,levels=c(ego1@result[["Description"]][15:1]))})
ggplot(data=C1,mapping = aes(x=fun,y=p,group=factor(1)))+geom_bar(stat = "identity",fill="#FF2D2D")+coord_flip()+
  theme(legend.position = "none")+theme(axis.text.x = element_text(size = 16,face = "plain",color="black"))+
  theme(axis.title.x = element_text(size = 15,face = "plain",color="black"))+
  theme(axis.text.y = element_text(size = 13,face = "plain",color = "black"))+
  xlab(NULL)+ylab("- log (p.adjust)")+scale_x_discrete(labels=function(x) str_wrap(x, width=50))