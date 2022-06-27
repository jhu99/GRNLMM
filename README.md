# GRNLMM: Constructing gene co-expression networks from single-cell expression data using linear mixed model


## Installation
For installation please use the following codes in R

```
install_github("jhu99/GRNLMM/grnlmm")
```
## Example
```
library(grnlmm)
load('data/expressiondata.rda')
load('data/vg.rda')
load('data/ve.rda')
Vg <- grnlmm (x, vg, ve)
```

### Input of GRNLMM
 <li type="disc">x&nbsp;:&nbsp;G x C matrix of expression data, where G is the number of genes and C is the number of cells</li>
 <li type="disc">V_g&nbsp;:&nbsp;G x G symmetric matrix of initial value of genetic covariance matrix</li>
 <li type="disc">V_e&nbsp;:&nbsp;G x G symmetric matrix of initial value of error covariance matrix</li>

### Output of GRNLMM 
The output of GRNLMM is a G x G covariance symmetric matrix V_g with the format :
```
 0.31  0.15 -0.43  ...
 0.15  1.50  0.60  ...
-0.43  0.60  1.13  ...
 .
 .
 .
```
where V_g [&nbsp;i&nbsp;,&nbsp;j&nbsp;] represents the correlation between gene i and gene j.

Further, we can transform it to a matrix R of correlation coefficient by
```
for(i in 1:nrow(A)){
    D[i,i] <- sqrt(A[i,i])  
}

Di <- solve(D)

R <- Di %*% A %*% Di
```

## Applications
The experimental code implementation in the paper can be viewed in applications in grnlmm folder.

## Citation
