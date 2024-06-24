setwd("C:/Users/huber/Box/Phillips Lab/Projects/A and C Metabolomic Trajectory Project/Manuscript/Githubcode/")

library(lme4)
library(car)
library(dplyr)

# Clear the global environment, but not the plots
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)


#read in data table
data <- read.csv("A-C Normalized Metabolomics Time Series Data.csv", header =T)
data$Day <- as.numeric(data$Day)
data$SimpleSel <- gsub("-Type", "", data$SimpleSel)


#make empty table
result_table <- data.frame(
  Compound = character(),
  Sel_pval = numeric(),
  Day_pval = numeric(),
  Interaction_pval=numeric()
)

#run loop for mixed effect model 
for (i in 8:ncol(data)){
  
  temp <- data[,c(1:4,i)]
  model <-  suppressMessages(lmer(temp[,5] ~ SimpleSel * Day + (1 | Pop), data = temp))
  
  test <- suppressMessages(Anova(model))
  
  
  result_table <- rbind(result_table, data.frame(
    Compound = colnames(temp[5]),
    Sel_pval = test$`Pr(>Chisq)`[1],
    Day_pval = test$`Pr(>Chisq)`[2],
    Interaction_pval = test$`Pr(>Chisq)`[3]
  )
  )
}

#adjust for multiple comparisons
result_table$Sel_fdr <- p.adjust(result_table$Sel_pval, method = "fdr")
result_table$Day_fdr <- p.adjust(result_table$Day_pval, method = "fdr")
result_table$Interaction_fdr <- p.adjust(result_table$Interaction_pval, method = "fdr")

write.csv(result_table, "Mixed_effect_results_AC_metabolomics_First80samples.csv", row.names = FALSE)


### Figure 5A ###

#load packages
library(ggplot2)
library(ggvenn)
library(dbplyr)


#read in data
df <- read.csv("A-C Normalized Metabolomics Time Series Data.csv", na = "N/A") #read in csv with relevant data

# Replace "9" with "09" in the combined Day
df$Day <- gsub("9", "09", df$Day)
df$Day <- as.factor(df$Day)
#rename column heading
colnames(df)[colnames(df) == "SimpleSel"] <- "Selection"
#change capitalization of "type"
df$Selection <- gsub("Type", "type", df$Selection)



df$combined <- paste(df$Selection, df$Day, sep = "_") #make new column combining type and day
# Use dplyr to move the 'combined' column to the 8th position
df <- df %>%
  select(-combined) %>%
  mutate(combined = paste(Selection, Day, sep = "_")) %>%
  select(1:7, combined, everything())


#pca
pca_df <-as.matrix(na.omit(df[,9:ncol(df)]))
pca_df <- prcomp(pca_df, scale=FALSE) 

pca.var_df <- pca_df$sdev^2
pca.var.per_df <- round(pca.var_df/sum(pca.var_df)*100, 1)

# Set labels for PCA
Type <- df$Day

# gradient
colors <- c("09" = "yellow3","21" = "orange", "28" = "orangered", "35"= "red3", "70" = "firebrick4")

# Create ggplot object and print for Selection
Run_PCA <- function() {
  pca.data_df <- data.frame(Sample=rownames(pca_df$x),
                            X=pca_df$x[,1],
                            Y=pca_df$x[,2])
  
  
  p_df <-ggplot(pca.data_df,aes(x=X, y=Y,color=Type, label=Type, shape=df$Selection)) +
    scale_color_manual(values = colors) +
    guides(shape = guide_legend(title = "Selection"),color = guide_legend(title = "Age from Egg (days)")) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
          legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
          panel.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_line(color = "gray90"),
          legend.key = element_rect(fill = "transparent", color = NA),
          axis.line = element_line(color = "black"))
  
  p_df<-p_df+geom_point() + 
    xlab(paste("PC1 - ", pca.var.per_df[1], "%", sep="")) + 
    ylab(paste("PC2 - ", pca.var.per_df[2], "%", sep="")) + 
    geom_point(size = 4) + 
    theme(legend.key.size = unit(1, 'cm'))
  
  p_df
}

Run_PCA()


### Create PCA including Subpopulation data

# Create ggplot object and print for Selection
Run_PCA_2 <- function() {
  pca.data_df <- data.frame(Sample=rownames(pca_df$x),
                            X=pca_df$x[,1],
                            Y=pca_df$x[,2])
  
  
  p_df2 <-ggplot(pca.data_df,aes(x=X, y=Y,color=Type, label=Type, shape=df$Sel)) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = c(19, 1, 17, 2)) +
    guides(shape = guide_legend(title = "Selection"),color = guide_legend(title = "Age from Egg (days)")) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
          legend.text = element_text(size = 12), legend.title = element_text(size = 12), 
          panel.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_line(color = "gray90"),
          legend.key = element_rect(fill = "transparent", color = NA),
          axis.line = element_line(color = "black"))
  
  p_df2<-p_df2+geom_point() + 
    xlab(paste("PC1 - ", pca.var.per_df[1], "%", sep="")) + 
    ylab(paste("PC2 - ", pca.var.per_df[2], "%", sep="")) + 
    geom_point(size = 4) + 
    theme(legend.key.size = unit(1, 'cm'))
  
  p_df2
}

Run_PCA_2()



### Figure 5B ###

library(ggplot2)
library(tidyverse)
library(ggvenn)
library(dplyr)
library(vegan)
library(UpSetR)

# Clear the global environment, but not the plots
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)


#### UpsetPlot ####
data <- read.csv("Mixed_effect_results_AC_metabolomics_First80samples.csv", header = TRUE)

all <- data[,2]
interaction <- subset(data , Interaction_fdr < 0.01 )[,2]
selection <- subset(data , Sel_fdr < 0.01 )[,2]
time <- subset(data , Day_fdr < 0.01 )[,2]

ven <- list(All = all, Selection =  selection, Age = time, Interaction = interaction)

# create upset plot with better colors
upset(fromList(ven), nintersects = NA, mainbar.y.max = 40, 
      sets.x.label = "Number of Metabolites", 
      mainbar.y.label = "Intersections",
      sets = c("Selection", "Age", "Interaction"), 
      sets.bar.color = "#56B4E9",
      main.bar.color = "black",
      matrix.color = "Orange3",
      order.by = "degree",
      decreasing = F, 
      empty.intersections = NULL, 
      text.scale = c(3, 3, 1.5, 2, 2, 3),
      point.size = 6, line.size = 2)

