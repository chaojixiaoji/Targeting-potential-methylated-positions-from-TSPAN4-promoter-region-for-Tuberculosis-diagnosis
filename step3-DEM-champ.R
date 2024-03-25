rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load(file = 'step1-output.Rdata')
myLoad    
load(file = 'step2-champ_myNorm.Rdata')
group_list=myLoad$pd$group
table(group_list)
myDMP <- champ.DMP(beta = myNorm,pheno=group_list, adjPVal = 1)
write.table(myDMP,file="myDMP.txt",sep="\t")
write.csv(myDMP,file = "myDMP.csv")
head(myDMP[[1]])
save(myDMP,file = 'step3-output-myDMP.Rdata')

myDMR <- champ.DMR(beta = myLoad$beta,pheno=group_list,method="Bumphunter",minProbes=7,cores=6,maxGap=200)
write.table(myDMR,file="myDMR min2.txt",sep = "\t")
write.csv(myDMR,file = "myDMR min2.csv")
DMR.GUI(beta = myLoad$beta,pheno = myLoad$pd$group)


