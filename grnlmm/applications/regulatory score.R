#calculate the regulatory score
i<-1
RC<-matrix(0,100,2)
rname<-tf
while(i<=sqrt(nrow(res))){
  j<-0
  while(j<sqrt(nrow(res))){
    if(res[(i-1)*100+1+j,3]>0.4 || res[(i-1)*100+1+j,3]<(-0.4)){
      RC[i,1]<-RC[i,1]+1
    }
    if(res[(i-1)*100+1+j,4]>0.4 || res[(i-1)*100+1+j,4]<(-0.4)){
      RC[i,2]<-RC[i,2]+1
    }
    j<-j+1
  }
  rname[i]<-res[(i-1)*100+1,1]
  i<-i+1
}
row.names(RC)<-rname
colnames(RC)<-c("normal","tumor") 
RC<-RC/100
RC<-RC[c(output1[,1],output2[,1],output3[,1],output4[,1],output5[,1]),]
RC<-RC[order(RC[,1],decreasing = T),]
RC<-RC[order(RC[,2],decreasing = T),]
RC<-RC[order(RC[,3],decreasing = T),]
RC<-RC[order(RC[,4],decreasing = T),]
RC<-RC[order(RC[,5],decreasing = T),]
RCrate<-matrix(0,414,1)
for (i in 1:414) {
  if(RC[i,2]>RC[i,1]){
    RCrate[i,1]<-RCrate[i,1]+(RC[i,2]-RC[i,1])
  }
  if(RC[i,3]>RC[i,2]){
    RCrate[i,1]<-RCrate[i,1]+(RC[i,3]-RC[i,2])
  }
  if(RC[i,4]>RC[i,3]){
    RCrate[i,1]<-RCrate[i,1]+(RC[i,4]-RC[i,3])
  }
}
row.names(RCrate)<-row.names(RC)
RCrate<-RCrate[order(RCrate[,1],decreasing = T),]
RCrate<-as.matrix(RCrate)
RC<-RC[c(row.names(RCrate)),]
library(pheatmap)
pheatmap(RC,cluster_rows = F,cluster_cols = F)