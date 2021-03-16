library(caret)

setwd("describeBIC")

classif = read.table("tcga_vs_myPAM50.txt", header=T)
nona = classif[which(! is.na(classif$TCGAprediction)),]
diff = nona[which(nona$TCGAprediction != nona$myPrediction),]

a = confusionMatrix(nona$TCGAprediction, nona$myPrediction)
a[["byClass"]][ , "F1"] 