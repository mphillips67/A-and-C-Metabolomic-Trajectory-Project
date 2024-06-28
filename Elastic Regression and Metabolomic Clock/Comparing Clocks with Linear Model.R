###
###   Zachary Greenspan
###   David Hubert
###



# Clear ggplot2# Clear the global environment, but not the plots
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

#load packages
library(lme4)

df_A <- read.csv("df_A.csv")
df_C <- read.csv("df_C.csv")

se <- function(x) sqrt(var(x) / length(x))

r2_score <- function(actual, predicted) {
  residSS <- sum((actual - predicted) ** 2)
  totalSS <- sum((actual - mean(actual)) ** 2)
  1 - (residSS / totalSS)
}

TSS <- function(actual)
  sum((actual - mean(actual)) ** 2)

RSS <- function(actual, predicted)
  sum((actual - predicted) ** 2)

dfs_comb <- function(df_CtoC, df_CtoA) {
  lm_CtoC <- lm(df_CtoC$pred_day ~ df_CtoC$actual_day)
  intercept_CtoC <- lm_CtoC$coefficients[1]
  slope_CtoC <- lm_CtoC$coefficients[2]
  r2_CtoC <- r2_score(df_CtoC$pred_day,
                      df_CtoC$actual_day)
  
  df_CtoA <- df_CtoA[df_CtoA$SimpleSel == "A-Type",]
  lm_CtoA <- lm(df_CtoA$pred_day ~ df_CtoA$actual_day)
  intercept_CtoA <- lm_CtoA$coefficients[1]
  slope_CtoA <- lm_CtoA$coefficients[2]
  r2_CtoA <- r2_score(df_CtoA$actual_day,
                      df_CtoA$pred_day)
  
  combined_Df <- rbind(df_CtoC, df_CtoA)
}

dfs_comb_A <- function(df_AtoA, df_AtoC) {
  lm_AtoA <- lm(df_AtoA$pred_day ~ df_AtoA$actual_day)
  intercept_AtoA <- lm_AtoA$coefficients[1]
  slope_AtoA <- lm_AtoA$coefficients[2]
  r2_AtoA <- r2_score(df_AtoA$pred_day,
                      df_AtoA$actual_day)
  
  df_AtoC <- df_AtoC[df_AtoC$SimpleSel == "C-Type",]
  lm_AtoC <- lm(df_AtoC$pred_day ~ df_AtoC$actual_day)
  intercept_AtoC <- lm_AtoC$coefficients[1]
  slope_AtoC <- lm_AtoC$coefficients[2]
  r2_AtoC <- r2_score(df_AtoC$actual_day,
                      df_AtoC$pred_day)
  
  combined_Df <- rbind(df_AtoA, df_AtoC)
}




predictingAgeSelectionData <- read.csv("predictingAgeSelectionData.csv")

c<- predictingAgeSelectionData$SimpleSel == "A-Type"
temp <- dfs_comb(df_C, predictingAgeSelectionData[c,1:4])

linearModel_C <- lm(pred_day ~ actual_day * SimpleSel2, data = temp)

c<- predictingAgeSelectionData$SimpleSel == "C-Type"
temp <- dfs_comb_A(df_A, predictingAgeSelectionData[c,1:4])

linearModel_A <- lm(pred_day ~ actual_day * SimpleSel2, data = temp)

summary(linearModel_A)$coefficients
summary(linearModel_C)$coefficients


summary(linearModel_A)
summary(linearModel_C)