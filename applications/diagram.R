#draw the scatter diagram of two genes' expression at several time point
library(ggplot2)
library(gcookbook)
library(ggpubr)
y<-matrix(0,188,2)
colnames(y)<-c("POU5F1","TERF1")
y[,1]<-x[34,571:758]
y[,2]<-x[59,571:758]
y<-as.data.frame(y)

g1<-ggplot(y1,aes(x=POU5F1,y=TERF1)) + geom_point(shape = 3,color="red") + ylim(0,10) + xlim (0,8)+
  annotate("text", x=7.2, y=9.5, label="0.59",size=5, fontface="plain",color = "red")
g2<-ggplot(y2,aes(x=POU5F1,y=TERF1)) + geom_point(shape = 3,color="green")  + ylim(0,10) + xlim (0,8)+
  annotate("text", x=7.2, y=9.5, label="0.50",size=5, fontface="plain",color = "green")
g3<-ggplot(y3,aes(x=POU5F1,y=TERF1)) + geom_point(shape = 3,color="blue")  + ylim(0,10) + xlim (0,8)+
  annotate("text", x=7.2, y=9.5, label="0.39",size=5, fontface="plain",color = "blue")
g4<-ggplot(y4,aes(x=POU5F1,y=TERF1)) + geom_point(shape = 3,color="orange")  + ylim(0,10) + xlim (0,8)+
  annotate("text", x=7.2, y=9.5, label="0.09",size=5, fontface="plain",color = "orange")
g5<-ggplot(y5,aes(x=POU5F1,y=TERF1)) + geom_point(shape = 3,color="pink")  + ylim(0,10) + xlim (0,8)+
  annotate("text", x=7.2, y=9.5, label="-0.35",size=5, fontface="plain",color = "pink")
ggarrange(g1,g2,g3,g4,g5,ncol=3,nrow=2,labels=c("0h","12h","24h","72h","96h"),label.x= 0,label.y = 0.17)

#line chart
library(ggplot2)
library(reshape2)
library(gridExtra)
y1<-matrix(0,5,3)
y1[,1]<-c(0,12,24,72,96)
y1[,2]<-c("ZFX","ZFX","ZFX","ZFX","ZFX")
y1[,3]<-c(mean(x["ZFX",1:92]),mean(x["ZFX",93:194]),mean(x["ZFX",195:260]),mean(x["ZFX",433:570]),mean(x["ZFX",571:758]))
y2<-matrix(0,5,3)
y2[,1]<-c(0,12,24,72,96)
y2[,2]<-c("NANOG","NANOG","NANOG","NANOG","NANOG")
y2[,3]<-c(mean(x["NANOG",1:92]),mean(x["NANOG",93:194]),mean(x["NANOG",195:260]),mean(x["NANOG",433:570]),mean(x["NANOG",571:758]))
y3<-matrix(0,5,3)
y3[,1]<-c(0,12,24,72,96)
y3[,2]<-c("SOX2","SOX2","SOX2","SOX2","SOX2")
y3[,3]<-c(mean(x["SOX2",1:92]),mean(x["SOX2",93:194]),mean(x["SOX2",195:260]),mean(x["SOX2",433:570]),mean(x["SOX2",571:758]))
y4<-matrix(0,5,3)
y4[,1]<-c(0,12,24,72,96)
y4[,2]<-c("TUB","TUB","TUB","TUB","TUB")
y4[,3]<-c(mean(x["TUB",1:92]),mean(x["TUB",93:194]),mean(x["TUB",195:260]),mean(x["TUB",433:570]),mean(x["TUB",571:758]))
y5<-matrix(0,5,3)
y5[,1]<-c(0,12,24,72,96)
y5[,2]<-c("CEBPZ","CEBPZ","CEBPZ","CEBPZ","CEBPZ")
y5[,3]<-c(mean(x["CEBPZ",1:92]),mean(x["CEBPZ",93:194]),mean(x["CEBPZ",195:260]),mean(x["CEBPZ",433:570]),mean(x["CEBPZ",571:758]))
y6<-matrix(0,5,3)
y6[,1]<-c(0,12,24,72,96)
y6[,2]<-c("ID4","ID4","ID4","ID4","ID4")
y6[,3]<-c(mean(x["ID4",1:92]),mean(x["ID4",93:194]),mean(x["ID4",195:260]),mean(x["ID4",433:570]),mean(x["ID4",571:758]))
y<-rbind(y1,y2,y3,y4,y5,y6)
y<-as.data.frame(y)
colnames(y)<-c("time","genes","expression")
y<-y[order(y[,1],decreasing = F),]
y[,3]<-as.numeric(y[,3])
ggplot(data=y,aes(x=time,y=expression,group=genes,color=genes)) +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),legend.position = )+
  geom_line(lwd=0.8) +scale_x_discrete(expand = c(0,0))+
  annotate("text", x=4.6, y=3.4, label="PITX2",size=4, fontface="plain")+
  annotate("text", x=4.6, y=1.9, label="GATA3",size=4, fontface="plain")+
  annotate("text", x=4.6, y=1.15, label="LEF1",size=4, fontface="plain")+
  annotate("text", x=4.6, y=0.45, label="ID4",size=4, fontface="plain")+
  xlab("time")+theme(axis.title.x = element_text(size = 13,face = "bold"))+
  ylab("expression")+theme(axis.title.y = element_text(size = 13,face = "bold"))

# Figure 3B
x<-as.matrix(read.table("data.txt", sep="\t"))
ti<-as.matrix(read.table("time.txt", sep="\t"))
tf<-as.matrix(read.table("tf.txt", sep="\t"))
row.names(x)<-tf
da<-x[c("SOX2","NANOG","ZFX"),]
da<-rbind(da,ti[,2])
da<-t(da)
da<-rbind(da,da,da)
da[759:1516,1]<-da[1:758,2]
da[1517:2274,1]<-da[1:758,3]
da[1517:2274,2]<-"ZFX"
da[759:1516,2]<-"NANOG"
da[1:758,2]<-"SOX2"
da<-da[,-c(3)]
colnames(da)<-c("expression","gene_name","pseudo_time")
da[,1]<-as.numeric(da[,1])
da[,3]<-as.numeric(da[,3])
da<-as.data.frame(da)
library(ggthemes)
ggplot(data = da,aes(x=pseudo_time,y=expression,colour=gene_name))+geom_point(shape=19,size=1,alpha=0.8)+
  theme_wsj()+
  scale_colour_wsj()+
  guides(size=guide_legend(title=NULL),colour=guide_legend(title=NULL))+
  theme(legend.text = element_text(size=14,face="plain"))+
  theme(plot.margin = unit(rep(2,4),'lines'))+geom_smooth(method = "loess")

#histogram (Figure 2C)
library(ggplot2)
library(gcookbook)
library(ggpubr)
#the number of edges needs to be calculated
y<-matrix(c(82,65,65,68,"GRNLMM","SCODE","GENIE3","scLink"),4,2)
colnames(y)<-c("edge","methods")
y<-as.data.frame(y)
y<-within(y,{methods<-factor(methods,levels=c("GRNLMM","scLink","GENIE3","SCODE"))})
y[,1]<-as.numeric(y[,1])
ggplot(data=y,mapping=aes(x=methods,y=edge,fill="#FF6633",group=factor(1)))+geom_bar(stat = "identity",width = 0.5)+
  theme_bw() + theme(panel.grid=element_blank(),axis.text.x = element_text( vjust=-2,hjust = 0.5,size = 16,face = "bold")) +
  xlab(NULL)+ylab("# true predicted edges")+theme(axis.text.y = element_text(size = 14,face = "bold"),axis.title.y = element_text(size=17,face = "plain",vjust = 3))+
  theme(legend.position = "none",legend.text = element_text(size=14,face="plain"),legend.title = element_blank())+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  theme(legend.background=element_rect(fill="white", colour="black"))+
  theme(plot.margin = unit(rep(1,4),'lines'))+ylim(c(0,90))
ggarrange(g1,g2,common.legend = T,legend = "right")

#correlation heatmap
A<-A[c(class1,class2,class3,class4,class5),]
A<-A[,c(class1,class2,class3,class4,class5)]
annotation_col = data.frame(GeneCluster = factor(rep(c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5"), c(15,29,16,20,20))))
row.names(annotation_col)<-row.names(A)
annotation_row = data.frame(GeneCluster = factor(rep(c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5"), c(15,29,16,20,20))))
row.names(annotation_row)<-row.names(A)
ann_colors=list(GeneCluster=c(Cluster1="#8CEA00",Cluster2="#FFC78E",Cluster3="#66B3FF",
                              Cluster4="#ff7575",Cluster5="#FFFF37"))
A<-apply(A, 1, as.numeric)
pheatmap(A,cluster_cols = F,cluster_rows = F,show_rownames = F,show_colnames = F,
         annotation_col = annotation_col,annotation_row = annotation_row,annotation_legend = FALSE,
         annotation_names_col = F,annotation_names_row = F,annotation_colors = ann_colors)

#violin diagram
my_comparisons<-list(c("NFtumor","NFnormal"),c("ATtumor","ATnormal"))
ggviolin(ceshi,x = "comp",y="expression",fill = "type")+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+xlab(NULL)+
  theme(axis.text.x = element_blank(),axis.text.y = element_text(size = 14),axis.title.y = element_text(size = 14))+
  theme(legend.text = element_text(size=14),legend.title = element_text(size = 14))