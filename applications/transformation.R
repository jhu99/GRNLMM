# It converts the covariance matrix of gcnvda's output to the correlation matrix
# Di is the inverse matrix of D
A<-read.csv("v_g.csv",header = TRUE, sep = ",", quote = "\"",
            dec = ".", fill = TRUE, comment.char = "!")
tf<-as.matrix(read.table("tf.txt", sep="\t"))
D<-matrix(0,nrow=nrow(A),ncol = ncol(A))
for(i in 1:nrow(A)){
  D[i,i]<-sqrt(A[i,i])  
}
Di<-solve(D)
A<-as.matrix(A)
R<-Di %*% A %*% Di
rownames(R)<-tf
colnames(R)<-tf