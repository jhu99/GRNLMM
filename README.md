# GRNLMM
Constructing gene co-expression networks from single-cell expression data using linear mixed model.
# Running GRNLMM
## Loading required functions
`source('grnlmm.R')`

## Usage
`V_g <- grnlmm(x,V_g,V_e)`

## Input data
x: G x C matrix of expression data, where G is the number of genes and C is the number of cells

V_g: G x G symmetric matrix of initial value of genetic covariance matrix

V_e: G x G symmetric matrix of initial value of error covariance matrix

## Example of input matrix x
` 2.1&nbsp;&nbsp;&nbsp;1.1&nbsp;&nbsp;&nbsp;0&nbsp;&nbsp;&nbsp;...

 3.6&nbsp&nbsp&nbsp0&nbsp&nbsp&nbsp0&nbsp&nbsp&nbsp...

.&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp

.&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp

.&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp`

## Example of input matrix V_g
` 1   0   0   ...`

` 0   1   0   ...`

`.`

`.`

`.           1`

## Example of input matrix V_e
` 1   0   0   ...`

` 0   1   0   ...`

`.`

`.`

`.           1`

#Output of GRNLMM
The output of GRNLMM is a G x G covariance symmetric matrix V_g, where V_g[1,2] represents the correlation of gene1 and gene2.

Further, we can transform it to a matrix R of correlation coefficient by

R = Di %*% V_g %*% Di, where Di is the inverse matrix of D, D = sqrt(V_g).
