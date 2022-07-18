# Filter ground truth based on tf set
reference<-reference[order(reference[,3],decreasing = T),]
re<-matrix(0,nrow(res),3)
k<-1
for (i in 1:nrow(reference)) {
  if(reference[i,1] %in% tf && reference[i,2] %in% tf){
    re[k,1]<-reference[i,1]
    re[k,2]<-reference[i,2]
    k<-k+1
  }
}

# add labels 
re<-re[1:(k-1),]
re<-paste(re[,1],re[,2],sep = "+")
re<-re[1:1000]
ge<-paste(res[,1],res[,2],sep = "+")
for (i in 1:nrow(res)) {
  if(ge[i] %in% re){
    res[i,4]<-1
  }
} 

#calculate the AUROC
res<-read.csv("vg_res.csv",header = TRUE, sep = ",", quote = "\"",
              dec = ".", fill = TRUE, comment.char = "!")
res<-res[order(abs(res[,3]),decreasing = T),]
# res<-res[-c(1:100),]
res[,3]<-abs(res[,3])
pred1<-prediction(res[,3],res[,4])
performance(pred1,"auc")@y.values[[1]]

#draw the PRC
plot(performance(pred1,"prec","rec"),lwd=1,col="#A6A600",ylim=c(0,0.2),cex.lab=1.3,add=T)
legend("topright",legend = c("GRNLMM","GENIE3","SCODE","scLink"),
       col=c("#FF5809","#46A3FF", "#FF44FF", "#A6A600"),lty = 1,lwd = 2)

#draw the ROC
library(ggplot2)
library(pROC)
#pay attention to the replacement of res
rocscode<-roc(res[,4],res[,3])
rocgenie3<-roc(res[,4],res[,3])
rocscLink<-roc(res[,4],res[,3])
rocour<-roc(res[,4],res[,3])
g2<-ggroc(list(GRNLMM=rocour,GENIE3=rocgenie3,SCODE=rocscode,scLink=rocscLink),lwd=0.8,legacy.axes=TRUE)
g2 +
  theme(legend.position=c(0.88,0.135))+
  theme(legend.text=element_text(size=14, face="plain"),legend.title = element_blank())+
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  theme(legend.background=element_rect(fill="white", colour="black"))+
  annotate("text", x=0.129, y=1.05, label="GRNLMM: 0.561",size=5, fontface="plain",) +
  annotate("text", x=0.143, y=1, label="GENIE3: 0.550",size=5, fontface="plain",) +
  annotate("text", x=0.145, y=0.95, label="SCODE: 0.541",size=5, fontface="plain",) +
  annotate("text", x=0.162, y=0.90, label="scLink: 0.501",size=5, fontface="plain",) +
  scale_colour_manual("methods",values=c("#FF5809","#46A3FF", "#FF44FF", "#A6A600"))+
  xlab("false positive rate")+theme(axis.title.x = element_text(size = 17,face = "plain",vjust = -5))+
  ylab("true positive rate")+theme(axis.title.y = element_text(size = 17,face = "plain",vjust = 8))+
  theme(axis.text.y = element_text(size = 14,face = "bold",angle = 90,vjust=3,hjust=0.5))+
  theme(axis.text.x = element_text(size = 14,face = "bold",vjust=-2,hjust=0.5))+
  theme(plot.margin = unit(rep(2,4),'lines'))

