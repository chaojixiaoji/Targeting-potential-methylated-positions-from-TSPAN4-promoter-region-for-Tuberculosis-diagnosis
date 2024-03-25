rm(list=ls())
a <- read.table("expreset.txt",header = T)
a[1:4,1:4]
class(a)
rownames(a) <- a[,1]
a <- a[,-1]
a <- as.matrix(a)
class(a)
a['GAPDH',]#查看GAPDH在样本中的表达量
a['ACTB',]#查看β-actin在样本中的表达量
boxplot(a)#画样本表达量的箱型图
pd <- read.table("group.txt",header=T)
b <- pd
rownames(b) <- b[,1]
myLoad <- list(beta = a,pd = b)
save(myLoad,file = 'step1-output.Rdata')
