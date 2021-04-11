# bladder-cancer-tumour-cell-phenotype-classifiers
IHC based classifier to predict major subtype of bladder cancer:

Using available dataset on ~30 proteins in ~300 samples from advanced bladder cancer tissue, a decision tree classifer model was made. The general steps that I followed to create these classifers are as:
### 1- data exploratory analysis (in R)
### 2- feature selection by ROC analysis (in R)
### 3- Creating a decision tree multiclass classifer (by sklearn library in Python) 
___________________________________________________________________________________________________________________________________________________________________

## 1- Data exploratory analysis

```R
# Loading libraries

library(factoextra)
library(FactoMineR)
library(missMDA)
library(corrplot)
library(ggplot2)
library(reshape)

# reading data
data2017 <- data.frame(data.table::fread("4.JPath2017_PtLevel.csv", header=T, sep=",", ))

# Removing those cases with unknow IHC subtype
data2017 <- data2017[-which(data2017$Subtype_IHC == "N/A"),]


# converting N/A to na
for (i in 1:41){
  data2017[,i][data2017[,i] == "N/A"] <- NA
}
#______Data distribution plots______________________________##
#Color by groups
new_df <- data2017[,c(9,11:41)]
# counting na
for (i in 1:length(colnames(new_df))){
  print(sum(is.na(new_df[,i])))
}
# convert to numeric
for (i in 2:length(colnames(new_df))){
  new_df[,i] <- as.numeric(new_df[,i])
}


plotDf <- melt(new_df, id.vars = 'Subtype_IHC')

# dropping non-common subtypes
table(plotDf$Subtype_IHC)

# BaSq    GU Mes-L Sc/NE   Uro 
# 1767  3379   465   558  6045 

plotDf <- plotDf[-which(plotDf$Subtype_IHC == "Mes-L" |plotDf$Subtype_IHC == "Sc/NE") ,]

# density plot to see how data distribution look like
ggplot(data = plotDf, aes(x=value, fill=Subtype_IHC))+
  geom_density(aes(x=value, y=..scaled.., fill=Subtype_IHC), alpha=0.5)+
  facet_wrap(~variable)+
  theme(legend.position="none")
```
