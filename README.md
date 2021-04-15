# bladder-cancer-tumour-cell-phenotype-classifiers
IHC based classifier to predict major subtype of bladder cancer:

Using available dataset on ~30 proteins in ~300 samples from advanced bladder cancer tissue, a decision tree classifer model was made. The general steps that I followed to create these classifers are as:
#### 1- Data exploratory analysis (in R)
#### 2- Feature selection by ROC analysis (in R)
#### 3- Creating a decision tree multiclass classifer (by sklearn library in Python) 
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

```R
# to see contribution of each variables to PCs
fviz_eig(new.pca, addlabels = TRUE, ylim = c(0, 30))
#
```
![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/scree_plot.JPG) 
```R
# correlation plot to see how features contributed to PCs
var <- get_pca_var(new.pca)
corrplot(var$cos2, is.corr=FALSE)
```
![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/corrplot.JPG) 

___________________________________________________________________________________________________________________________________________________________________

## 2- Feature selection 

In order to identify which features can potentially be useful to classify samples into diffrent subtypes a set of binay ROC analysis were performed. 

```R
#
data2017 <- data.frame(data.table::fread("4.JPath2017_PtLevel.csv", header=T, sep=",", ))

#
data2017 <- data2017[-which(data2017$Subtype_IHC == "N/A"),]

table(data2017$Subtype_IHC)

# convert preotein expression value to numeric vectors
for (i in 11:length(colnames(data2017))){
  data2017[,i] <- as.numeric(data2017[,i])
}


UroBasalGu_2017 <- data2017[which(data2017$Subtype_IHC == "BaSq" | data2017$Subtype_IHC == "GU" |data2017$Subtype_IHC == "Uro"),]

# Uro+GU = Luminal
LumBas_2017 <-UroBasalGu_2017
LumBas_2017$Subtype_IHC <- ifelse(LumBas_2017$Subtype_IHC == "BaSq", "basal", "luminal")

table(LumBas_2017$Subtype_IHC)

#basal luminal 
#57     304 
prop.table(table(LumBas_2017$Subtype_IHC))
#basal   luminal 
#0.1578947 0.8421053



library(pROC)

proteins <- names(LumBas_2017)[11:41]
# define empty vectors for auc, protein name
auc <-c()
protein_name <- c()
auc_ci_ll <- c()
auc_ci_ul <- c()


for (i in 1:length(proteins)){
  df <- data.frame(subtype = as.factor(LumBas_2017$Subtype_IHC), 
                   prot = as.numeric(LumBas_2017[,proteins[i]]))
  names(df)[2] <- proteins[i]
  d <- roc(subtype ~ df[,2], data=df, percent = TRUE, ci = TRUE, auc= TRUE)
  protein_name[i] <- proteins[i]
  auc[i] <- d$auc
  auc_ci_ll[i] <- round(data.frame(d$ci)[1,],2)
  auc_ci_ul[i] <- round(data.frame(d$ci)[3,],2)
}

LumBas_2017_auc <- data.frame(protein = protein_name,
                              AUC = auc,
                              LL_CI =auc_ci_ll ,
                              UL_CI = auc_ci_ul)
```
![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/luminal_basal_auc.JPG)
```R


UroGu_2017 <- UroBasalGu_2017[-which(UroBasalGu_2017$Subtype_IHC == "BaSq"),]

proteins_UroGu <- names(UroGu_2017)[11:41]

# for loop
# define empty vectors for auc, protein name
auc <-c()
protein_name <- c()
auc_ci_ll <- c()
auc_ci_ul <- c()


for (i in 1:length(proteins_UroGu)){
  df <- data.frame(subtype = as.factor(UroGu_2017$Subtype_IHC), 
                   prot = as.numeric(UroGu_2017[,proteins_UroGu[i]]))
  names(df)[2] <- proteins_UroGu[i]
  d <- roc(subtype ~ df[,2], data=df, percent = TRUE, ci = TRUE, auc= TRUE)
  protein_name[i] <- proteins_UroGu[i]
  auc[i] <- d$auc
  auc_ci_ll[i] <- round(data.frame(d$ci)[1,],2)
  auc_ci_ul[i] <- round(data.frame(d$ci)[3,],2)
}

UroGu_2017_auc <- data.frame(protein = protein_name,
                             AUC = auc,
                             LL_CI =auc_ci_ll ,
                             UL_CI = auc_ci_ul)
```
![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/uroGuAUC.JPG)

```R
GuBas_2017 <- UroBasalGu_2017[-which(UroBasalGu_2017$Subtype_IHC == "Uro"),]


proteins_GuBas <- names(GuBas_2017)[11:41]

# for loop
# define empty vectors for auc, protein name
auc <-c()
protein_name <- c()
auc_ci_ll <- c()
auc_ci_ul <- c()


for (i in 1:length(proteins_GuBas)){
  df <- data.frame(subtype = as.factor(GuBas_2017$Subtype_IHC), 
                   prot = as.numeric(GuBas_2017[,proteins_GuBas[i]]))
  names(df)[2] <- proteins_GuBas[i]
  d <- roc(subtype ~ df[,2], data=df, percent = TRUE, ci = TRUE, auc= TRUE)
  protein_name[i] <- proteins_GuBas[i]
  auc[i] <- d$auc
  auc_ci_ll[i] <- round(data.frame(d$ci)[1,],2)
  auc_ci_ul[i] <- round(data.frame(d$ci)[3,],2)
}

GuBas_2017_auc <- data.frame(protein = protein_name,
                             AUC = auc,
                             LL_CI =auc_ci_ll ,
                             UL_CI = auc_ci_ul)


```

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/GU_Basal_AUC.JPG)

___________________________________________________________________________________________________________________________________________________________________
## 3- Creating decision tree multiclass classifers
Considering AUC analysis these markers found to do good job in classifying subtypes: 

"KRT14", "KRT5", "CDH3", "FOXA1", "GATA3", "PPARG", "RB1" , "CCND1", "CDKN2Ap16", "FGFR3", "TP63"


```python
# Import libraries
import numpy as np
import pandas as pd
import os as os 

#reading dataset
data = pd.read_csv("4.JPath2017_PtLevel.csv")
data.head()
```
![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/dataset_1.JPG)

```python
# Review missing data in the demographic information
print(data.isnull().sum())

#TMA#                      0
#Age                       0
#Sex                       0
#Grade                     6
#Stage                     1
#Clinical stage           28
#Hyb#                    118
#Subtype_LundTax_RNA     120
#Subtype_IHC              31
#CK5_Pattern_Bernardo    175
#CCNB1                     9
#CCND1                     4
#CD3density                2
#CD68density               0
#CDH1                      5
#CDH3                      8
#CDKN2Ap16                 8
#CHGA                     15
#E2F3                      6
#EGFR                     25
#EPCAM                    16
#ERBB2                     8
#ERBB3                     7
#FGFR3                    24
#FOXA1                    21
#FOXM1                    11
#GATA3                     8
#KRT14                     5
#KRT20                    25
#KRT5                      6
#NCAM1                     7
#PPARG                     3
#pSTAT3                    7
#RB1                       8
#RXRA                      6
#SYP                      12
#TP63                      5
#TUBB2B                   20
#UPK3                     22
#VIM                      18
#ZEB2                     30
#dtype: int64

# dropping NAs from Subtype_IHC
data = data.dropna(subset=['Subtype_IHC'])
#print(data.isnull().sum())

print(data.shape)
#(394, 41)

# slicing based on columns
d = data.iloc[:, 10:41]

# imputing missing values by mean of columns
d.fillna(d.mean(), inplace = True)

#print(d.isnull().sum())

# replace a part of data farme by othe dataframe
data.iloc[:, 0:9]

data = pd.concat([data.iloc[:, 0:9], d], axis=1)

#print(data.isnull().sum())

# drop values from subtype_IHC column
data = data.drop(data[(data['Subtype_IHC'] == 'Mes-L') | (data['Subtype_IHC'] == 'Sc/NE')].index)

# feature selection based on ROC analysis:
## KRT14, KRT5, CDH3, FOXA1, GATA3, PPARG for Uro vs. Bas
## RB1 , CCND1, CDKN2Ap16, FGFR3

corrmat = data[["Subtype_IHC", "KRT14", "KRT5", "CDH3", "FOXA1", "GATA3", "PPARG", "RB1" , "CCND1", "CDKN2Ap16", "FGFR3", "TP63"]].corr()
# visualization 
import seaborn as sn
import matplotlib.pyplot as plt
sn.heatmap(corrmat, annot=True)
plt.show()
```
![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/corrmat.JPG)

```python
#Making classifiers based on selected features
X = data[["KRT14", "KRT5", "RB1"]].values
y = data.iloc[:,8].values
#y
# encode IHC subtype as a dummy variable
y,class_names = pd.factorize(y)

# Splitting the dataset into the Training set and Test set
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.25, stratify=y, random_state = 42)

# Fitting Classifier to the Training Set
from sklearn.tree import DecisionTreeClassifier
classifier = DecisionTreeClassifier(criterion='entropy',max_depth=3, random_state=42)
classifier.fit(X_train, y_train)

#DecisionTreeClassifier(criterion='entropy', max_depth=3, random_state=42)

# Model performance on training set
y_pred_train =classifier.predict(X_train)

from sklearn import metrics
from sklearn.metrics import confusion_matrix, classification_report

accuracy = metrics.accuracy_score(y_train, y_pred_train)
print("Accuracy: {:.2f}".format(accuracy))
cm=confusion_matrix(y_train,y_pred_train)
print('Confusion Matrix: \n', cm)
print(classification_report(y_train, y_pred_train))
```
```python
#___________________Model performance on training set_________________________#
#Accuracy: 0.92
#Confusion Matrix: 
# [[ 38   5   0]
# [  0 143   3]
# [  0  13  68]]
#              precision    recall  f1-score   support
#
#           0       1.00      0.88      0.94        43
#           1       0.89      0.98      0.93       146
#           2       0.96      0.84      0.89        81
#
#    accuracy                           0.92       270
#   macro avg       0.95      0.90      0.92       270
#weighted avg       0.93      0.92      0.92       270
```

```python
# Predicting the test results
y_pred=classifier.predict(X_test)

# Classification results on test set
from sklearn import metrics
accuracy = metrics.accuracy_score(y_test, y_pred)
print("Accuracy: {:.2f}".format(accuracy))

from sklearn.metrics import confusion_matrix, classification_report
cm=confusion_matrix(y_test,y_pred)
print('Confusion Matrix: \n', cm)
print(classification_report(y_test, y_pred))
```
```python
#___________________Model performance on test set_________________________#
Accuracy: 0.89
Confusion Matrix: 
 [[11  3  0]
 [ 0 49  0]
 [ 2  5 21]]
              precision    recall  f1-score   support

           0       0.85      0.79      0.81        14
           1       0.86      1.00      0.92        49
           2       1.00      0.75      0.86        28

    accuracy                           0.89        91
   macro avg       0.90      0.85      0.87        91
weighted avg       0.90      0.89      0.89        91
```


```python
# Visualize the tree by graphiz
import graphviz
from sklearn import tree
from IPython.display import Image
from IPython.display import Image  

feature_names = ["KRT14", "KRT5","RB1" ]
dot_data = tree.export_graphviz(classifier, out_file=None, filled=True,feature_names = feature_names, rounded = True, class_names=class_names)
graph = graphviz.Source(dot_data)
graph
#graph.render('round-table.pdf', view=True)  
