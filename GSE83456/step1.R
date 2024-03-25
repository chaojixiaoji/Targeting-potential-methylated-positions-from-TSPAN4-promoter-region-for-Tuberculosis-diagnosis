
rm(list=ls())
library(limma)
x <- read.ilmn("GSE83456 non_normalized.txt",
               expr = "sample",
               probeid = "ID_REF")
y <- neqc(x)
exp=y$E

dat=rbind(ID=colnames(exp),exp)
write.table(dat,file="norexp.txt",sep="\t",quote = F,col.names = F)


cols=rainbow(ncol(exp))
pdf(file="nor.pdf",width = 6,height = 4.5)
par(cex=0.3,mar=c(8,8,8,8))
boxplot(exp[,91:106],col=cols)
dev.off()

#对探针进行注释
library(AnnoProbe)
ids=idmap('GPL10558',type = 'bioc')
head(ids)
colnames(ids)=c('probe_id','symbol')  
ids=ids[ids$symbol != '',]
ids=ids[ids$probe_id %in%  rownames(dat),]

dat[1:4,1:4]   
dat=dat[ids$probe_id,] 

###去除重复探针，保留中位数最大的探针与基因相对应
ids$median=apply(dat,1,median) #ids新建median这一列，列名为median，同时对dat这个矩阵按行操作，取每一行的中位数，将结果给到median这一列的每一行
ids=ids[order(ids$symbol,ids$median,decreasing = T),]
#对ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
ids=ids[!duplicated(ids$symbol),]#将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
dat=dat[ids$probe_id,] #新的ids取出probe_id这一列，将dat按照取出的这一列中的每一行组成一个新的dat
rownames(dat)=ids$symbol#把ids的symbol这一列中的每一行给dat作为dat的行名
dat[1:4,1:4]  #保留每个基因ID第一次出现的信息(median最大的那个)
dat=rbind(ID=colnames(dat),dat)#手动添加一行列名，写出文件的时候，选择不写出列名，就不会出现列名错位的情况啦
write.table(dat,file="expreset.txt",sep = "\t",col.names = FALSE,quote = FALSE)
