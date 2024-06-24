library(readxl)
library(bit)

### Reading in File and Setting up

metabolomics <- read_excel("Metabolomics Raw Data.xlsx", sheet = "Relative Quant Data")
metabolomicsSampleInfo <- read_excel("Metabolomics Raw Data.xlsx", sheet = "Sample Info")

metabolomicsPruned <- metabolomics[,4:89]
metabolomicsNames <- metabolomics[,1]

### Reading in File and Setting up



### Data Formatting - Removing Partial and Empty Rows

metabolomicsByIself <- data.frame(metabolomicsPruned)
for (i in 1:393) {
  metabolomicsByIself[i,] <- as.integer(metabolomicsByIself[i,])
}

whichNumsBad <- 2:394
for (i in 1:393) {
  if (sum(is.na(metabolomicsByIself[i,])) == 0) {
    whichNumsBad[i] <- TRUE
  }
  else
    whichNumsBad[i] <- FALSE
}
whichNumsBad <- as.booltype(whichNumsBad)

metabolomicsMod <- data.frame(t(metabolomicsPruned[whichNumsBad,])[,1:202])

### Data Formating - Removing Partial and Empty Rows



### Scaling by Sample - 
### Removing Spiked Stable Isotope Labeled Internal Standards (SILISs)

for (i in 1:86) {
  temp <- as.double(metabolomicsMod[i,1:202])
  temp <- log10(temp)
  temp <- ( temp - mean(temp) ) / sd(temp)
  metabolomicsMod[i,] <- temp
}

### Scaling by Sample - 
### Removing Spiked Stable Isotope Labeled Internal Standards (SILISs)



### Name Formating to Remove Special Characters for Final Table
### Handling Duplicates after Said Process

inBetweenNames <- metabolomicsNames$COMPOUND[whichNumsBad]
nameVector <- 1:202
for (i in 1:202) {
  nameVector[i] <- gsub('[[:digit:]]', '', inBetweenNames[i])
  nameVector[i] <- gsub('[[:punct:]]', '', nameVector[i])
  nameVector[i] <- gsub(' ', '', nameVector[i])
  nameVector[duplicated(nameVector)] <- paste(nameVector[duplicated(nameVector)], 2, sep = '')
}

### Name Formating to Remove Special Characters for Final Table
### Handling Duplicates after Said Process



### Final Table Experimental Variables Organization

colnames(metabolomicsMod) <- nameVector

nameTable <- data.frame(SimplifiedNames = nameVector)
nameTable$ExpandedNames <- inBetweenNames[1:202]

finalTable <- data.frame('Pop' = metabolomicsSampleInfo$Population)
finalTable$Day <- metabolomicsSampleInfo$Day
finalTable$Sel <- 1
finalTable$SimpleSel <- 1
finalTable$RunOrder <- metabolomicsSampleInfo$`NWMRC ID`
finalTable$Batch <- 1
finalTable$Batch <- 3
finalTable$Batch[finalTable$RunOrder < 61] <- 2
finalTable$Batch[finalTable$RunOrder < 31] <- 1
finalTable$NumberOfFliesInSample <- metabolomicsSampleInfo$`number of flies`

### Final Table Experimental Variables Organization



### Adding KEGG IDs to File

#metabolomics$COMPOUND %in% nameTable$ExpandedNames
ourKEGG <- metabolomics$`KEGG ID`[metabolomics$COMPOUND %in% 
                                       nameTable$ExpandedNames]
nameTable$KEGG_ID <- ourKEGG
nameTable$First_KEGG_Id <- substr(ourKEGG, 1, 6)

### Adding KEGG IDs to File



### Adding Extra Population Structure Labeling

for (i in 1:86) {
  temp <- substr(finalTable$Pop[i], 1, 2)
  if (temp == "CO") {
    finalTable$Sel[i] <- "CO"
    finalTable$SimpleSel[i] <- "C-Type"
  }
  else if (temp == "NC") {
    finalTable$Sel[i] <- "nCO"
    finalTable$SimpleSel[i] <- "C-Type"
  }
  else if (temp == "AO") {
    finalTable$Sel[i] <- "AO"
    finalTable$SimpleSel[i] <- "A-Type"
  }
  else if (temp == "AC") {
    finalTable$Sel[i] <- "ACO"
    finalTable$SimpleSel[i] <- "A-Type"
  }
}

finalTable <- cbind(finalTable, metabolomicsMod)

for (i in 8:209) {
  finalTable[,i] <- as.double(finalTable[,i])
}

### Adding Extra Population Structure Labeling



### Normalize for Batch Effects

finalTable$Batch <- as.factor(finalTable$Batch)
finalTableRaw <- finalTable
for (i in 1:202) {
  formula <- paste(nameVector[i], ' ~ Batch', sep = '')
  formula <- as.formula(formula)
  x <- lm(formula, finalTable)
  y <- residuals(x)
  finalTable[,i + 7] <- y
}

### Normalize for Batch Effects



### Remove Samples Not in Experiment

finalTable2 <- finalTable[metabolomicsSampleInfo$`Sample ID` < 81,]

### Remove Samples Not in Experiment



### Write to File

write.csv(finalTable2, "A-C Normalized Metabolomics Time Series Data.csv", row.names=FALSE)
write.csv(nameTable, "A-C Metabolites Names.csv", row.names=FALSE)

### Write to File
