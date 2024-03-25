rm(list=ls())
library(limma)              
expFile="output.txt"    
conFile="sample1.txt"           
treatFile="sample2.txt"           
setwd("E:\\wenzhang\\TB meth blood\\完整流程\\lasso")    


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

s1=read.table(conFile, header=F, sep="\t", check.names=F)
sampleName1=as.vector(s1[,1])
conData=data[,sampleName1]


s2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName2=as.vector(s2[,1])
treatData=data[,sampleName2]


rt=cbind(conData, treatData)



treatNum=ncol(treatData)
conNum=ncol(conData)
Type=c(rep("Control",conNum),rep("Treat",treatNum))
outData=rbind(id=paste0(colnames(rt),"_",Type),rt)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)


