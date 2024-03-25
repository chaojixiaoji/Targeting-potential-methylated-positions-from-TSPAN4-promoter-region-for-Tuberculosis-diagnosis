rm(list = ls())
setwd("E:\\wenzhang\\TB meth blood\\完整流程")
options(stringsAsFactors = F)

require(GEOquery)
require(Biobase)
library("impute")

info=read.table("group.txt",sep="\t",header=T)
library(data.table)
b=info
rownames(b)=b[,1]

a=fread("genematrix.txt",data.table = F )
a[1:4,1:4]
rownames(a)=a[,1]
a=a[,-1]
beta=as.matrix(a)
beta=impute.knn(beta)
betaData=beta$data
betaData=betaData+0.00001
a=betaData
a[1:4,1:4]
identical(colnames(a),rownames(b))

library(ChAMP)

myLoad=champ.filter(beta = a,pd = b, arraytype = "450K")
myLoad
save(myLoad,file = 'step1-output.Rdata')
