
library(e1071)
library(caret)
library(pROC)
library(ggplot2)

rm(list=ls())

mydata <- read.table("output SVM.txt",row.names = 1,header = T)


table(mydata$status)
X <- mydata[,2:57] 
Y <- as.numeric(as.factor(mydata$status))

control <- trainControl(method = "repeatedcv",   
                       
                        number = 5,     
                       
                        repeats = 1,     
                       
                        search = "random"   
                          
                        )
  

set.seed(12345)  
svm_rfe <- rfe(X,
               Y,
               sizes = c(1:50),
               rfeControl = rfeControl(functions = caretFuncs,
                                      
                                       method = "repeatedcv",
                                      
                                       number = 5,
                                      
                                       repeats = 1,
                                      
                                       verbose = FALSE),
                                    
               method = "svmLinear",
        
               trControl = control,
           
               preProc = c("center", "scale")
      
               )


svm_rfe_ranking <- svm_rfe$variables
head(svm_rfe_ranking)



varImp(svm_rfe)

varImp_dataframe <- data.frame(Gene = row.names(varImp(svm_rfe))[1:10],
                               importance = varImp(svm_rfe)[1:10, 1])
                          
write.table(svm_rfe_ranking,file="SVM gene.txt",sep = "\t")

varImp_dataframe <- na.omit(varImp_dataframe)


mycolors <- c('#D4E2A7','#88D7A4','#A136A1','#BAE8BC','#C757AF',
              '#DF9FCE','#D5E1F1','#305691','#B6C2E7','#E8EFF7',
              '#9FDFDF','#EEE0F5','#267336','#98CEDD','#CDE2EE',
              '#DAD490','#372E8A','#4C862D','#81D5B0','#BAE8C9',
              '#A7DCE2','#AFDE9C')

ggplot(varImp_dataframe, aes(x = reorder(Gene, -importance), y = importance , fill = Gene)) + 

  geom_col() +

  ggtitle("Hub Genes") +
  
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 12, color = "black",angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.title = element_text(margin = margin(b = 20)),
        panel.grid.major = element_line(color = "grey", size = 0.2)) +
  xlab("position") + ylab("Importance") +
  scale_fill_manual(values = mycolors)



top_10_vars <- svm_rfe_ranking$var[1:10]
top_10_vars


top10_SVM_data <- mydata[,top_10_vars]


X_plot = svm_rfe$results$Variables
Y_plot = svm_rfe$results$RMSE
plot(X_plot, Y_plot,  
     xlab="Variable Number",  
     ylab="RMSE (Cross-Validation)",  
     col="#81D5B0",  
     pch=16,                
     cex=1.5,          
     lwd=3,            
     type="b",         
     ylim=c(0.3, 0.7)) 

abline(h=min(Y_plot), col="skyblue")  
grid(col="grey",lwd=1,lty=3)  
 
legend("topright",c("Cross-Validation RMSE"),
       col=c("#81D5B0"),pch=c(16,NA),lwd=2,bg="white")   # 添加图例


wmin <- which.min(svm_rfe$results$RMSE)
wmin


points(wmin, svm_rfe$results$RMSE[wmin], col = "orange", pch = 16, cex=1.5)  
text(wmin, svm_rfe$results$RMSE[wmin],  
     paste0("N=", wmin),
     pos = 4, col = "orange", cex=1.5)  


Target_Genes <- svm_rfe$optVariables
Target_Genes


Best_SVM_data <- mydata[,Target_Genes]

