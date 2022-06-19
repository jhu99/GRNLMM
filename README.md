# GRNLMM
Constructing gene co-expression networks from single-cell expression data using linear mixed model.
# Running GRNLMM

## Installation
`install_github("jhu99/GRNLMM/grnlmm")`


## Input data
x: G x C matrix of expression data, where G is the number of genes and C is the number of cells

V_g: G x G symmetric matrix of initial value of genetic covariance matrix

V_e: G x G symmetric matrix of initial value of error covariance matrix


## Example
`library(grnlmm)`
`load('data/expressiondata.rda')`
`load('data/vg.rda')`
`load('data/ve.rda)`
`Vg <- grnlmm (x, vg, ve)`


# Output of GRNLMM
The output of GRNLMM is a G x G covariance symmetric matrix V_g, where V_g[1,2] represents the correlation of gene1 and gene2.

Further, we can transform it to a matrix R of correlation coefficient by

R = Di · V_g · Di, where Di is the inverse matrix of D, D = sqrt (V_g).
