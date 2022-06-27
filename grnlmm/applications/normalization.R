#normalization
x<-read.csv("expression.csv",header = TRUE, sep = ",", quote = "\"",
                dec = ".", fill = TRUE, comment.char = "!")
g<-read.csv("g.csv",header = TRUE, sep = ",", quote = "\"",
            dec = ".", fill = TRUE, comment.char = "!")
x<-x[,-c(1)]
x<-t(x)
row.names(x)<-x[,1]
x<-x[,-c(1)]
x<-apply(x,2,as.numeric)
x<-sweep(x, MARGIN = 1, 1e6/rowSums(x), FUN = "*") 
x<-log10(x + 1)
x<-t(x)