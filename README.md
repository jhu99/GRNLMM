# GCNVDA: Inference of gene co-expression networks from single-cell transcriptome data based on variance decomposition analysis

![](https://img.shields.io/github/r-package/v/jhu99/GCNVDA)
![](https://img.shields.io/github/license/jhu99/GCNVDA)
[![](https://img.shields.io/badge/downloads-108-green)](https://github.com/jhu99/GCNVDA/graphs/traffic)
![](https://img.shields.io/github/stars/jhu99/gcnvda?style=social)

&emsp;&emsp;We development a new method, GCNVDA, which models single-cell expression data by variance decomposition and uses the covariance matrix of random effect terms to characterize the correlation between genes. To overcome the influence of randomness of intercellular expression and improve the accuracy of the predicted GCNs, we use a known correlation matrix to reflect the relationship between cells and add a noise term to the model. Our results show that GCNVDA has advantages in accurately identifying the co-expressive relationships between genes and can explore genes and gene function modules that play an essential role in biological processes.

## Installation
For installation please use the following codes in R

```
install_github("jhu99/GCNVDA")
```
## Example
```
library(gcnvda)
load('data/expressiondata.rda')
load('data/vg.rda')
load('data/ve.rda')
Vg <- gcnvda (x, vg, ve)
```

### Input of GCNVDA
 <li type="disc">x&nbsp;:&nbsp;G x C matrix of expression data, where G is the number of genes and C is the number of cells</li>
 <li type="disc">V_g&nbsp;:&nbsp;G x G symmetric matrix of initial value of genetic covariance matrix</li>
 <li type="disc">V_e&nbsp;:&nbsp;G x G symmetric matrix of initial value of error covariance matrix</li>

### Output of GCNVDA 
The output of GCNVDA is a G x G covariance symmetric matrix V_g with the format :
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
The experimental code implementation in the paper can be viewed in applications folder.

## Citation
Jialu Hu, Bin Lian, Haohui Zhang, Tao Wang, Yongtian Wang, Xuequn Shang, Inference of gene co-expression networks from single-cell transcriptome data based on variance decomposition analysis (unpublished).
