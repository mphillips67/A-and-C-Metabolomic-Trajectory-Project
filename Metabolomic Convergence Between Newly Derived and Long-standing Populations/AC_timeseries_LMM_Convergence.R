
library(lme4)
library(car)
library(dplyr)
library(qvalue)
library(MuMIn)

# Clear the global environment, but not the plots
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

#read in data table
data <- read.csv("A-C Normalized Metabolomics Time Series Data.csv", header =T)
data$Day <- as.numeric(data$Day)
data$SimpleSel <- gsub("-Type", "", data$SimpleSel)

#A type comparison
A_type <- subset(data, SimpleSel =="A")

#make empty table
result_table_Atype <- data.frame(
  Compound = character(),
  Sel_pval = numeric(),
  Day_pval = numeric(),
  Interaction_pval=numeric()
)

for (i in 8:ncol(A_type)){
  
  temp <- A_type[,c(1:4,i)]
  model <-  suppressMessages(lmer(temp[,5] ~ Sel * Day + (1 | Pop), data = temp))
  
  test <- suppressMessages(Anova(model))
  
  result_table_Atype <- rbind(result_table_Atype, data.frame(
    Compound = colnames(temp[5]),
    SubPop_pval = test$`Pr(>Chisq)`[1],
    Day_pval = test$`Pr(>Chisq)`[2],
    Interaction_pval = test$`Pr(>Chisq)`[3]
  )
  )
}

#adjust for multiple comparisons
result_table_Atype$SubPop_fdr <- p.adjust(result_table_Atype$SubPop_pval, method = "fdr")
result_table_Atype$Day_fdr <- p.adjust(result_table_Atype$Day_pval, method = "fdr")
result_table_Atype$Interaction_fdr <- p.adjust(result_table_Atype$Interaction_pval, method = "fdr")


#C type comparison
#A type comparison
C_type <- subset(data, SimpleSel =="C")

#make empty table
result_table_Ctype <- data.frame(
  Compound = character(),
  Sel_pval = numeric(),
  Day_pval = numeric(),
  Interaction_pval=numeric()
)

for (i in 8:ncol(C_type)){
  
  temp <- C_type[,c(1:4,i)]
  model <-  suppressMessages(lmer(temp[,5] ~ Sel * Day + (1 | Pop), data = temp))
  
  test <- suppressMessages(Anova(model))
  
  result_table_Ctype <- rbind(result_table_Ctype, data.frame(
    Compound = colnames(temp[5]),
    SubPop_pval = test$`Pr(>Chisq)`[1],
    Day_pval = test$`Pr(>Chisq)`[2],
    Interaction_pval = test$`Pr(>Chisq)`[3]
  )
  )
}

#adjust for multiple comparisons
result_table_Ctype$SubPop_fdr <- p.adjust(result_table_Ctype$SubPop_pval, method = "fdr")
result_table_Ctype$Day_fdr <- p.adjust(result_table_Ctype$Day_pval, method = "fdr")
result_table_Ctype$Interaction_fdr <- p.adjust(result_table_Ctype$Interaction_pval, method = "fdr")


#write.csv(result_table_Ctype, file = "LMM_C-Type_results.csv")
#write.csv(result_table_Atype, file = "LMM_A-Type_results.csv")
