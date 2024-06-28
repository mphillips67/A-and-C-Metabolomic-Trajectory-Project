###
###   Zachary Greenspan
###

library(utils)
library(ggplot2)
library(glmnet)
library(caret)
library(dplyr)
library(foreach)
library(doParallel)
library(gridExtra)
library(ggpubr)
library(ggtext)
library(gridExtra)
library(ggpubr)


#no_cores <- parallel::detectCores() - 1
#cluster <- makePSOCKcluster(no_cores)
#registerDoParallel(cluster)

#unregister_dopar <- function() {
#  env <- foreach:::.foreachGlobals
#  rm(list=ls(name=env), pos=env)
#}

metaData <- read.csv("A-C Normalized Metabolomics Time Series Data.csv")


####
####
#### Basic Functions Start
####
####


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

stratifiedFolds <- function(numSamples, kFolds, stratas) {
  numEachType <- numSamples / stratas
  result <- sample(1:numEachType %% kFolds + 1)
  if (stratas > 1)
    for (i in 1:(stratas - 1)) {
      result <- c(result, sample(1:numEachType %% kFolds + 1))
    }
  result
}


####
####
#### Basic Functions End
####
####



####
####
#### Machine Learning Functions Start
####
####


runStratifiedML <- function(metaData, phenoData, phenoDay, populations, selection, method = "glmnet", folds = c(21, 22), numPops = 20) {
  result <- NULL
  phenos <- NULL
  sels <- NULL
  days <- NULL
  test_pops <- NULL
  result_train <- NULL
  phenos_train <- NULL
  sels_train <- NULL
  days_train <- NULL
  pops_train <- NULL
  models <- NULL
  if (folds[1] == 21)
    x <- stratifiedFolds(numPops, 5, 2)
  else
    x <- folds
  pops <- sort(unique(populations))
  
  
  tuneGrid = expand.grid(
    alpha=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
    lambda=c(seq(0.000001,1,length=100),seq(1,200,length=201))
  )
  
  for (i in 1:5) {
    Omics_train <- metaData[populations %in% pops[x != i],]
    Pheno_train <- phenoData[populations %in% pops[x != i]]
    
    model<-train(
      Omics_train, 
      Pheno_train,
      tuneGrid = tuneGrid,
      method = method
    )
    
    Omics_test <- metaData[populations %in% pops[x == i],]
    Pheno_test <- phenoData[populations %in% pops[x == i]]
    
    pred <- predict(model, Omics_test)
    pred_train <- predict(model, Omics_train)
    
    result <- c(result, pred)
    phenos <- c(phenos, Pheno_test)
    sels <- c(sels, selection[populations %in% pops[x == i]])
    days <- c(days, phenoDay[populations %in% pops[x == i]])
    test_pops <- c(test_pops, populations[populations %in% pops[x == i]])
    
    result_train <- c(result_train, pred_train)
    phenos_train <- c(phenos_train, Pheno_train)
    sels_train <- c(sels_train, selection[populations %in% pops[x != i]])
    days_train <- c(phenoDay, selection[populations %in% pops[x != i]])
    pops_train <- c(pops_train, populations[populations %in% pops[x != i]])
  
    models <- c(models, model)
  }
  end_train <- list()
  end_train[[1]] <- result_train
  end_train[[2]] <- phenos_train
  end_train[[3]] <- sels_train
  end_train[[4]] <- days_train
  end_train[[5]] <- pops_train
  end_test <- list()
  end_test[[1]] <- result
  end_test[[2]] <- phenos
  end_test[[3]] <- sels
  end_test[[4]] <- days
  end_test[[5]] <- test_pops
  end <- list()
  end[[1]] <- end_train
  end[[2]] <- end_test
  end[[3]] <- models
  end
}

runStratifiedML.multi <- function(metaData, phenoData, phenoDay, populations, selection, method = "glmnet", folds = c(21, 22), nPer = 5, time = FALSE) {
  result <- NULL
  phenos <- NULL
  sels <- NULL
  days <- NULL
  test_pops <- NULL
  result_train <- NULL
  phenos_train <- NULL
  sels_train <- NULL
  days_train <- NULL
  pops_train <- NULL
  if (folds[1] == 21)
    x <- stratifiedFolds(20, 5, 2)
  else
    x <- folds
  pops <- sort(unique(populations))
  
  tuneGrid = expand.grid(
    alpha=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
    lambda=c(seq(0.000001,1,length=100),seq(1,200,length=201))
  )
  
  if (time)
    start_time <- Sys.time()
  
  models.list <- foreach(i=1:nPer,.errorhandling = c('pass')) %dopar%{ 
    Omics_train <- metaData[populations %in% pops[x != i],]
    Pheno_train <- phenoData[populations %in% pops[x != i]]
    
    model<-train(
      Omics_train, 
      Pheno_train,
      tuneGrid = tuneGrid,
      method = method
    )
  }
  if (time) {
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  
  for (i in 1:5) {
    Omics_train <- metaData[populations %in% pops[x != i],]
    Pheno_train <- phenoData[populations %in% pops[x != i]]
    
    Omics_test <- metaData[populations %in% pops[x == i],]
    Pheno_test <- phenoData[populations %in% pops[x == i]]
    
    pred <- predict(models.list[[i]], Omics_test)
    pred_train <- predict(models.list[[i]], Omics_train)
    
    #result[i] <- r2(pred, Pheno_test)
    result <- c(result, pred)
    phenos <- c(phenos, Pheno_test)
    sels <- c(sels, selection[populations %in% pops[x == i]])
    days <- c(days, phenoDay[populations %in% pops[x == i]])
    test_pops <- c(test_pops, populations[populations %in% pops[x == i]])
    
    result_train <- c(result_train, pred_train)
    phenos_train <- c(phenos_train, Pheno_train)
    sels_train <- c(sels_train, selection[populations %in% pops[x != i]])
    days_train <- c(phenoDay, selection[populations %in% pops[x != i]])
    pops_train <- c(pops_train, populations[populations %in% pops[x != i]])
  }
  
  end_train <- list()
  end_train[[1]] <- result_train
  end_train[[2]] <- phenos_train
  end_train[[3]] <- sels_train
  end_train[[4]] <- days_train
  end_train[[5]] <- pops_train
  end_test <- list()
  end_test[[1]] <- result
  end_test[[2]] <- phenos
  end_test[[3]] <- sels
  end_test[[4]] <- days
  end_test[[5]] <- test_pops
  end <- list()
  end[[1]] <- end_train
  end[[2]] <- end_test
  end
}

runPermutationStratified <- function(metaData, phenoData, phenoDay, populations, selection, method = "glmnet", folds = c(21, 22), nPer = 5, time = FALSE, nSim = 100) {
  result <- list()
  startTime <- Sys.time()
  for (i in 1:nSim) {
    if(i %% 10 == 0) {
      print(i)
      print(curTime - startTime)
    }
    result[[i]] <- runStratifiedML.multi(metaData, phenoData, phenoDay, populations, selection, method = "glmnet", folds = c(21, 22), nPer = 5, time = FALSE)
    curTime <- Sys.time()
  }
  result
}

plotPredictionAgeAll <- function(using, theme = theme_classic, title = "", r2 = FALSE) {
  tempDf <- data.frame(pred = using[[2]][[1]])
  tempDf$Pheno_test <- using[[2]][[2]]
  tempDf$SimpleSel <- using[[2]][[3]]
  tempDf$Days <- using[[2]][[4]]
  
  tempDf$combined <- paste(tempDf$SimpleSel, tempDf$Days, sep = "_")
  tempDf$combined <- gsub("9", "09", tempDf$combined)
  
  r2_reg <- r2_score(tempDf$Pheno_test, tempDf$pred)
  
  #colors <- c("A-Type_09" = "firebrick4","A-Type_21" = "red3", "A-Type_28" = "firebrick1","A-Type_35" = "orange", 
  #            "C-Type_21" = "midnightblue","C-Type_28" = "blue","C-Type_35" = "deepskyblue3", "C-Type_70" = "lightskyblue" )
  
  g1 <- ggplot(tempDf) +
    theme_light() +
    geom_point(mapping = aes(x = Pheno_test, y = pred,
                             color = SimpleSel,
                             shape = SimpleSel),
               #position = position_jitterdodge()
               size = 3) +
    geom_abline(slope = 1, intercept = 0) +
    #scale_color_manual(values = colors) +
    scale_colour_manual(name = "Legend", values = c("#BB271A", "#0000F3")) +
    theme_light() +
    theme(legend.position = c(.75,.85),
          #legend.key.size = unit(1.1, "cm"),
          legend.box.background = element_rect(colour = "black")) +
    labs(title = title, 
         x = "Actual Age", 
         y = "Predicted Age", 
         color = "Legend", 
         shape = "Legend", 
         group = "Legend") +
    scale_x_continuous(breaks = seq(10, 70, 10)) +
    scale_y_continuous(breaks = seq(10, 70, 10)) +
    theme
  
  if (r2)
    g1 <- g1 + annotate(geom = "label", x = 16, y = 70, 
                        label = paste('R^2 ==', round(r2_reg, 3), sep = ""), color  = "black",
                        size = 5, parse = TRUE)
  
  g1
}

makeDualModels <- function(Omics_A, Omics_C, pheno_A, pheno_C, method) {
  tuneGrid = expand.grid(
    alpha=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
    lambda=c(seq(0.000001,1,length=100),seq(1,200,length=201))
  )
  
  modelA<-train(
    Omics_A, 
    pheno_A,
    tuneGrid = tuneGrid,
    method = method,
    trControl = trainControl(
      method = "cv", 
      number = 5, 
      verboseIter = TRUE
    )
  )
  
  modelC<-train(
    Omics_C, 
    pheno_C,
    tuneGrid = tuneGrid,
    method = method,
    trControl = trainControl(
      method = "cv", 
      number = 5, 
      verboseIter = TRUE
    )
  )
  
  result <- list()
  result[[1]] <- modelA
  result[[2]] <- modelC
  result
}

makePredictiveDf <- function(Omics_A, Omics_C, pheno_A, pheno_C, method = "glmnet") {
  models <- makeDualModels(Omics_A, Omics_C, pheno_A, pheno_C, method)
  model_A <- models[[1]]
  model_C <- models[[2]]
  pred_A <- predict(model_C, Omics_A)
  pred_C <- predict(model_A, Omics_C)
  
  df_A <- data.frame(pred_day = pred_A)
  df_A$actual_day <- pheno_A
  df_A$SimpleSel <- "A-Type"
  df_A$SimpleSel2 <- "Predicting A-Type Age\nfrom C-Type Data"
  df_C <- data.frame(pred_day = pred_C)
  df_C$actual_day <- pheno_C
  df_C$SimpleSel <- "C-Type"
  df_C$SimpleSel2 <- "Predicting C-Type Age\nfrom A-Type Data"
  df_Both_A_C <- rbind(df_A, df_C)
  
  df_Both_A_C
}

makePredictiveDf_from_models <- function(model_A, model_C, Omics_A, Omics_C, pheno_A, pheno_C) {
  pred_A <- predict(model_C, Omics_A)
  pred_C <- predict(model_A, Omics_C)
  
  df_A <- data.frame(pred_day = pred_A)
  df_A$actual_day <- pheno_A
  df_A$SimpleSel <- "A-Type"
  df_A$SimpleSel2 <- "Predicting A-Type Age\nfrom C-Type Data"
  df_C <- data.frame(pred_day = pred_C)
  df_C$actual_day <- pheno_C
  df_C$SimpleSel <- "C-Type"
  df_C$SimpleSel2 <- "Predicting C-Type Age\nfrom A-Type Data"
  df_Both_A_C <- rbind(df_A, df_C)
  
  df_Both_A_C
}

makePredictiveDf2 <- function(using) {
  tempDf <- data.frame(pred_day = using[[2]][[1]])
  tempDf$actual_day <- using[[2]][[2]]
  tempDf$SimpleSel <- using[[2]][[3]]
  tempDf$SimpleSel2 <- using[[2]][[3]]
  tempDf
}

plot_Trained_On_CType <- function(df_CtoC, df_CtoA, theme = theme_classic, title = "", r2 = FALSE, r2_b = "No") {
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
  
  g1 <- ggplot(combined_Df, mapping = aes(x = actual_day, y = pred_day,
                                        color = SimpleSel,
                                        shape = SimpleSel,
                                        group = SimpleSel)) +
    #geom_abline(slope = 1, intercept = 0, linewidth = 1.5) +
    geom_point(size = 3.5) +
    geom_abline(slope = slope_CtoA, intercept = intercept_CtoA, color = "#BB271A") +
    geom_abline(slope = slope_CtoC, intercept = intercept_CtoC, color = "#0000F3") +
    scale_colour_manual(name = "Legend", values = c("#BB271A", "#0000F3")) +
    #scale_shape_manual(name = "Group", values = c("#BB271A", "#0000F3")) +
    theme_light() +
    theme(legend.position = c(.75,.85),
          #legend.key.size = unit(1.1, "cm"),
          legend.box.background = element_rect(colour = "black")) +
    labs(title = title, 
         x = "Actual Age", 
         y = "Predicted Age", 
         color = "Legend",
         shape = "Legend",
         group = "Legend") +
    scale_x_continuous(breaks = seq(10, 70, 10)) +
    scale_y_continuous(breaks = seq(10, 80, 10)) +
    ylim(1, 80) +
    xlim(1, 80) +
    theme #+ 
  
  if (r2) {
    label1.r <- paste('R<sup>2</sup> = ', round(r2_CtoA, 3))#, "<br>")
    label1.a <- paste('α = ', round(slope_CtoA, 3), "<br>")
    label1.b <- paste('β = ', round(intercept_CtoA, 3))
    label1 <- paste0(label1.r, label1.a, label1.b)
    label2.r <- paste('R<sup>2</sup> = ', round(r2_CtoC, 3))#, "<br>")
    label2.a <- paste('α = ', round(slope_CtoC, 3), "<br>")
    label2.b <- paste('β = ', round(intercept_CtoC, 3))
    label2 <- paste0(label2.r, label2.a, label2.b)
    g1 <- g1 + 
      ggtext::geom_richtext(aes(x,y, label = label),
                            hjust = 0,
                            color  = "#BB271A",
                            data = data.frame(x = 16, y = 70,
                                              SimpleSel = "A-Type",
                                              label = label1.r)) + 
      ggtext::geom_richtext(aes(x,y, label = label),
                            hjust = 0,
                            color  = "#0000F3",
                            data = data.frame(x = 50, y = 30,
                                              SimpleSel = "C-Type",
                                              label = label2.r))
  }
  if (r2_b == "Half") {
    label1.r <- paste('R<sup>2</sup> = ', round(r2_CtoA, 3))#, "<br>")
    label1.a <- paste('α = ', round(slope_CtoA, 3), "<br>")
    label1.b <- paste('β = ', round(intercept_CtoA, 3))
    label1 <- paste0(label1.r, label1.a, label1.b)
    label2.r <- paste('R<sup>2</sup> = ', round(r2_CtoC, 3))#, "<br>")
    label2.a <- paste('α = ', round(slope_CtoC, 3), "<br>")
    label2.b <- paste('β = ', round(intercept_CtoC, 3))
    label2 <- paste0(label2.r, label2.a, label2.b)
    g1 <- g1 + 
      ggtext::geom_richtext(aes(x,y, label = label),
                            hjust = 0,
                            color  = "#0000F3",
                            data = data.frame(x = 50, y = 30,
                                              SimpleSel = "C-Type",
                                              label = label2.r))
  }
  
  g1
}

plot_Trained_On_AType <- function(df_AtoA, df_AtoC, theme = theme_classic, title = "", r2 = FALSE, r2_b = "No") {
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
  
  g1 <- ggplot(combined_Df, mapping = aes(x = actual_day, y = pred_day,
                                          color = SimpleSel,
                                          shape = SimpleSel,
                                          group = SimpleSel)) +
    #geom_abline(slope = 1, intercept = 0, linewidth = 1.5) +
    geom_point(size = 3.5) +
    geom_abline(slope = slope_AtoA, intercept = intercept_AtoA, color = "#BB271A") +
    geom_abline(slope = slope_AtoC, intercept = intercept_AtoC, color = "#0000F3") +
    scale_colour_manual(name = "Legend", values = c("#BB271A", "#0000F3")) +
    #scale_shape_manual(name = "Group", values = c("#BB271A", "#0000F3")) +
    theme_light() +
    theme(legend.position = c(.75,.85),
          #legend.key.size = unit(1.1, "cm"),
          legend.box.background = element_rect(colour = "black")) +
    labs(title = title, 
         x = "Actual Age", 
         y = "Predicted Age", 
         color = "Legend",
         shape = "Legend",
         group = "Legend") +
    scale_x_continuous(breaks = seq(10, 70, 10)) +
    scale_y_continuous(breaks = seq(10, 80, 10)) +
    ylim(1, 80) +
    xlim(1, 80) +
    theme #+ 
  
  if (r2) {
    label1.r <- paste('R<sup>2</sup> = ', round(r2_AtoA, 3))#, "<br>")
    label1.a <- paste('α = ', round(slope_AtoA, 3), "<br>")
    label1.b <- paste('β = ', round(intercept_AtoA, 3))
    label1 <- paste0(label1.r, label1.a, label1.b)
    label2.r <- paste('R<sup>2</sup> = ', round(r2_AtoC, 3))#, "<br>")
    label2.a <- paste('α = ', round(slope_AtoC, 3), "<br>")
    label2.b <- paste('β = ', round(intercept_AtoC, 3))
    label2 <- paste0(label2.r, label2.a, label2.b)
    g1 <- g1 + 
      ggtext::geom_richtext(aes(x,y, label = label),
                            hjust = 0,
                            color  = "#BB271A",
                            data = data.frame(x = 16, y = 40,
                                              SimpleSel = "A-Type",
                                              label = label1.r)) + 
      ggtext::geom_richtext(aes(x,y, label = label),
                            hjust = 0,
                            color  = "#0000F3",
                            data = data.frame(x = 50, y = 20,
                                              SimpleSel = "C-Type",
                                              label = label2.r))
  }
  if (r2_b == "Half") {
    label1.r <- paste('R<sup>2</sup> = ', round(r2_AtoA, 3))#, "<br>")
    label1.a <- paste('α = ', round(slope_AtoA, 3), "<br>")
    label1.b <- paste('β = ', round(intercept_AtoA, 3))
    label1 <- paste0(label1.r, label1.a, label1.b)
    label2.r <- paste('R<sup>2</sup> = ', round(r2_AtoC, 3))#, "<br>")
    label2.a <- paste('α = ', round(slope_AtoC, 3), "<br>")
    label2.b <- paste('β = ', round(intercept_AtoC, 3))
    label2 <- paste0(label2.r, label2.a, label2.b)
    g1 <- g1 + 
      ggtext::geom_richtext(aes(x,y, label = label),
                            hjust = 0,
                            color  = "#BB271A",
                            data = data.frame(x = 16, y = 40,
                                              SimpleSel = "A-Type",
                                              label = label1.r))
  }
  
  g1
}

plotPredictionAgeSelection <- function(dataFrame, theme = theme_classic, title = "", r2 = FALSE) {
  cAfromC <- dataFrame$SimpleSel2 == "Predicting A-Type Age\nfrom C-Type Data"
  cCfromA <- dataFrame$SimpleSel2 == "Predicting C-Type Age\nfrom A-Type Data"
  
  r2_AfromC <- r2_score(dataFrame$actual_day[cAfromC],
                        dataFrame$pred_day[cAfromC])
  
  r2_CfromA <- r2_score(dataFrame$actual_day[cCfromA],
                        dataFrame$pred_day[cCfromA])
  
  g1 <- ggplot(dataFrame, mapping = aes(x = actual_day, y = pred_day,
                                  color = SimpleSel,
                                  shape = SimpleSel,
                                  group = SimpleSel)) +
    geom_point(size = 3.5) +
    geom_abline(slope = 1, intercept = 0) +
    scale_colour_manual(name = "Legend", values = c("#BB271A", "#0000F3")) +
    #scale_shape_manual(name = "Group", values = c("#BB271A", "#0000F3")) +
    theme_light() +
    theme(legend.position = c(.75,.85),
          #legend.key.size = unit(1.1, "cm"),
          legend.box.background = element_rect(colour = "black")) +
    labs(title = title, 
         x = "Actual Age", 
         y = "Predicted Age", 
         color = "Legend",
         shape = "Legend",
         group = "Legend") +
    scale_x_continuous(breaks = seq(10, 70, 10)) +
    scale_y_continuous(breaks = seq(10, 70, 10)) +
    theme #+ 
  
  if (r2) {
    g1 <- g1 + 
      annotate(geom = "label", x = 16,y = 70,
               label = paste('R^2 ==', round(r2_AfromC, 3), sep = ""),
               color  = "#BB271A", size = 5, parse = TRUE) + 
      annotate(geom = "label", x = 60,y = 30, 
               label = paste('R^2 ==', round(r2_CfromA, 3), sep = ""),
               color = "#0000F3", size = 5, parse = TRUE)
  }
  
  g1
}

plotPredictionMortality <- function(using) {
  tempDf <- data.frame(pred = using[[2]][[1]])
  tempDf$Pheno_test <- using[[2]][[2]]
  tempDf$SimpleSel <- using[[2]][[3]]
  tempDf$Days <- using[[2]][[4]]
  
  tempDf$combined <- paste(tempDf$SimpleSel, tempDf$Days, sep = "_")
  tempDf$combined <- gsub("9", "09", tempDf$combined)
  
  r2_reg <- r2_score(tempDf$Pheno_test, tempDf$pred)
  
  colors <- c("A-Type_09" = "firebrick4","A-Type_21" = "red3", "A-Type_28" = "firebrick1","A-Type_35" = "orange", 
              "C-Type_21" = "midnightblue","C-Type_28" = "blue","C-Type_35" = "deepskyblue3", "C-Type_70" = "lightskyblue" )
  
  g1 <- ggplot(tempDf) +
    theme_light() +
    geom_point(mapping = aes(x = Pheno_test, y = pred,
                             color = combined,
                             shape = SimpleSel),
               #position = position_jitterdodge()
               size = 3) +
    geom_abline(slope = 1, intercept = 0) +
    scale_color_manual(values = colors) +
    #scale_colour_manual(name = "Legend", values = c("#BB271A", "#0000F3")) +
    theme_light() +
    theme(#legend.position = c(.175,.85),
      #legend.key.size = unit(1.1, "cm"),
      legend.box.background = element_rect(colour = "black")) +
    labs(title = "", 
         x = "Actual Age-Specific Mortality", 
         y = "Predicted Age-Specific Mortality", 
         color = "Legend", 
         shape = "Legend", 
         group = "Legend") +
    machineLearningTheme# + 
  #annotate(geom = "label", x = 16,y = 58, 
  #         label = paste("R2 = ", round(r2_reg, 3), sep = ""), color  = "#BB271A",
  #         size = 5)
  
  g1
}

####
####
#### Machine Learning Functions End
####
####



####
####
#### Variables Start
####
####


Omics_A_Day <- metaData[metaData$SimpleSel == "A-Type", 8:209]
Omics_C_Day <- metaData[metaData$SimpleSel == "C-Type", 8:209]

Pheno_A_Day <- metaData[metaData$SimpleSel == "A-Type",]$Day
Pheno_C_Day <- metaData[metaData$SimpleSel == "C-Type",]$Day

Pop_A_Day <- metaData[metaData$SimpleSel == "A-Type",]$Pop
Pop_C_Day <- metaData[metaData$SimpleSel == "C-Type",]$Pop

Sel_A_Day <- metaData[metaData$SimpleSel == "A-Type",]$SimpleSel
Sel_C_Day <- metaData[metaData$SimpleSel == "C-Type",]$SimpleSel

textSize <- 18
machineLearningTheme = theme(
  axis.title.x = element_text(size = textSize, face = "bold"),
  axis.text.x = element_text(size = textSize - 2),
  axis.title.y = element_text(size = textSize, face = "bold"),
  axis.text.y = element_text(size = textSize - 2),
  legend.title = element_text(size = textSize - 4), 
  legend.text = element_text(size = textSize - 6),
  #legend.position = c(.85,.4),
  #legend.box.background = element_rect(colour = "black")), 
  legend.background = element_rect(fill="white",
                                   linewidth=0.5, linetype="solid", 
                                   colour ="black"),
  plot.title = element_text(size=18, hjust = 0.5, face = "bold"))


####
####
#### Variables End
####
####



####
####
#### Machine Learning Runs and Graphs Start
####
####



set.seed(22222)

predictingAgeFullData <- runStratifiedML(metaData[,8:209], metaData[,2], metaData$Day, 
                                         metaData$Pop, metaData$SimpleSel)

predictingAgeCFullData <- runStratifiedML(Omics_C_Day, Pheno_C_Day, Pheno_C_Day, 
                                          Pop_C_Day, Sel_C_Day, numPops = 10,
                                          folds =  stratifiedFolds(10, 10, 1))

predictingAgeAFullData <- runStratifiedML(Omics_A_Day, Pheno_A_Day, Pheno_A_Day, 
                                          Pop_A_Day, Sel_A_Day, numPops = 10,
                                          folds =  stratifiedFolds(10, 10, 1))

predictingAgeSelectionData <- makePredictiveDf(Omics_A_Day, Omics_C_Day, 
                                               Pheno_A_Day, Pheno_C_Day,
                                               method = "glmnet")

df_C <- makePredictiveDf2(predictingAgeCFullData)

df_A <- makePredictiveDf2(predictingAgeAFullData)

g1 <- plotPredictionAgeAll(predictingAgeFullData, machineLearningTheme, "", r2 = TRUE)

g2 <- plotPredictionAgeSelection(predictingAgeSelectionData, machineLearningTheme, "" ,r2 = TRUE)

ggarrange(g1, g2, nrow=2, common.legend = FALSE, legend="right")

g3 <- g1 + geom_label(label = "A", x = 3.5, y = 82, size = 13) + coord_cartesian(clip="off")
g4 <- g2 + geom_label(label = "B", x = 3.5, y = 81, size = 13, color = "black") + coord_cartesian(clip="off")

ggarrange(g3, g4, nrow=2, common.legend = FALSE, legend="right")

theme_panelTitle <- theme(
  plot.title = element_text(size=14, hjust = 0, face = "plain"))
g5 <- g1 + ggtitle("A. Prediction on Full Data Set") + theme_panelTitle
g6 <- g2 + ggtitle("B. Prediction on Selection Regime Data Sets") + theme_panelTitle

ggarrange(g5, g6, nrow=2, common.legend = FALSE, legend="right")

g7 <- plot_Trained_On_CType(df_C, predictingAgeSelectionData, theme = machineLearningTheme, r2 = FALSE)
g8 <- plot_Trained_On_AType(df_A, predictingAgeSelectionData, theme = machineLearningTheme, r2 = FALSE)

#final <- ggarrange(g7, g8, nrow=1, ncol=2, common.legend = FALSE, legend="right")
final <- ggarrange(g7, g8, nrow=1, ncol=2, common.legend = TRUE, legend="bottom")
ggsave(final, filename = "2024_06_03_Machine_Learning_Graph_No_R2_V2.pdf", device = "pdf", width = 15, height = 7.5, units = "in")
ggsave(final, filename = "2024_06_03_Machine_Learning_Graph_No_R2_V4.pdf", device = "pdf", width = 8, height = 4, units = "in")

figure <- ggarrange(g7 + rremove("ylab"), g8 + rremove("ylab"), ncol=2, common.legend = TRUE, legend="bottom")

annotate_figure(figure, left = ggpubr::text_Grob("Predicted Age", rot = 90, vjust = 1, gp = gpar(cex = 2)))


g7b <- plot_Trained_On_CType(df_C, predictingAgeSelectionData, theme = machineLearningTheme, r2 = FALSE, r2_b = "Half")
g8b <- plot_Trained_On_AType(df_A, predictingAgeSelectionData, theme = machineLearningTheme, r2 = FALSE, r2_b = "Half")

final <- ggarrange(g7b, g8b, nrow=2, common.legend = FALSE, legend="right")
ggsave(final, filename = "2024_05_24_Machine_Learning_Graph_w_Half_R2.pdf", device = "pdf", width = 7, height = 10, units = "in")

g7c <- plot_Trained_On_CType(df_C, predictingAgeSelectionData, theme = machineLearningTheme, r2 = TRUE)
g8c <- plot_Trained_On_AType(df_A, predictingAgeSelectionData, theme = machineLearningTheme, r2 = TRUE)

final <- ggarrange(g7c, g8c, nrow=2, common.legend = FALSE, legend="right")
ggsave(final, filename = "2024_05_24_Machine_Learning_Graph_w_All_R2.pdf", device = "pdf", width = 7, height = 10, units = "in")

theme_panelTitle <- theme(
  plot.title = element_text(size=14, hjust = 0, face = "plain"))
g09 <- g7 + ggtitle("A. Predictions from C-Type Data") + theme_panelTitle
g10 <- g8 + ggtitle("B. Predictions from A-Type Data") + theme_panelTitle

ggarrange(g09, g10, nrow=2, common.legend = FALSE, legend="right")

#7x10 inches

####
####
#### Machine Learning Runs and Graphs End
####
####

unregister_dopar()
