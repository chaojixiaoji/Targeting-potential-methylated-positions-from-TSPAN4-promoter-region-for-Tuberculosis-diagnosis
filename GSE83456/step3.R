rm(list=ls())
load(file = 'step1-output.Rdata')
exprSet <- myLoad$beta
library(pheatmap)
library(AnnoProbe)
group_list=myLoad$pd$group
group_list
#差异分析
library(limma)#就是给表达矩阵、分组矩阵、比较矩阵，limma给你分析差异结果
design=model.matrix(~0+factor(group_list))#分组矩阵
design#分组矩阵制作好了，用0和1代替
colnames(design)=levels(factor(group_list))#给分组矩阵的列名用分组信息来表示
rownames(design)=colnames(exprSet)#行名用样本名来表示

contrast.matrix=makeContrasts(TB-HC,levels = design)
contrast.matrix #比较矩阵，makeContrasts用来告诉谁比谁
fit=lmFit(exprSet,design)#lmFit函数为线性模型拟合，给表达矩阵和分组矩阵，输出
fit2 <- contrasts.fit(fit, contrast.matrix) #把上一步的线性拟合矩阵矩阵和比较矩阵给它，计算差异和标准误
fit2 <- eBayes(fit2)#通过eBayes将标准误调整到一个相同的水平再来计算t、logFC等值
tempOutput = topTable(fit2, coef=1, n=Inf)#用topTable提取
nrDEG = na.omit(tempOutput) #该函数是删去缺少某些值的对象

nrDEG <- nrDEG[nrDEG$P.Value < 0.05, ]#提取P<0.05基因的差异分析结果

head(nrDEG)
write.csv(nrDEG,"limma_notrend.results.csv",quote = F)#将差异分析结果保存到本地
save(nrDEG,file="step2-DEG.Rdata")
