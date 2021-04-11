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
![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/histograms.JPG)

```R
###__PCA analysis to see how samples cluster togather_____##


# Imputing data for missing ones
new_df <- imputePCA(new_df[,-1], ncp = 2)
new_df <- new_df$completeObs
new_df <- data.frame(cbind(data2017[,c(1:9)], new_df))

# doing PCA
new.pca <- PCA(new_df[,c(10:40)], graph = FALSE)


fviz_pca_ind(new.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = new_df$Subtype_IHC, # color by groups
             palette = "jco",
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups")
```
![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/pca1.JPG)

```R

# making a subset: uro + basal
binary_df <- new_df[which(new_df$Subtype_IHC == "BaSq" | new_df$Subtype_IHC == "Uro"),]
new.pca <- PCA(binary_df[,c(10:40)], graph = FALSE)

fviz_pca_ind(new.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = binary_df$Subtype_IHC, # color by groups
             palette = "jco",
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups")
```
![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/pca2.JPG)
