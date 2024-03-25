
rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load(file = 'step1-output.Rdata')

myLoad  
# 耗时步骤，运行一次后，就注释掉
if(T){
  myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
  dim(myNorm) 
  pD=myLoad$pd
  save(myNorm,pD,file = 'step2-champ_myNorm.Rdata')
}
load(file = 'step2-champ_myNorm.Rdata')
# 原来的450K经过质控过滤后是400K啦
beta.m=myNorm
group_list=myLoad$pd$group
dim(beta.m) 

