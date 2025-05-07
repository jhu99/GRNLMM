#confusion matrix
res<-res[,1:4]
la<-replicate(nrow(res),0)
res<-cbind(res,la)
for (i in 1:100) {
  res[i,5]<-1
}
tp=NA
fp=NA
fn=NA
tp = sum(res[,5]==1 & res[,4]==1)
fp = sum(res[,5]==1 & res[,4]==0)
fn = sum(res[,5]==0 & res[,4]==1)
#precision
P = tp/(tp+fp)
#recall
R = tp/(tp+fn)
#F1-score
f1 = 2*(P*R)/(P+R)