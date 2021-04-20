# bladder-cancer-tumour-cell-phenotype-classifiers
IHC based classifiers to predict major subtype of bladder cancer:
Using the available dataset on ~30 proteins in ~300 samples from advanced bladder cancer tissue, a decision tree classifier model was made. The general steps that I followed to create these classifiers are as follows:
#### 1- Data exploratory analysis 
#### 2- Feature selection by ROC analysis 
#### 3- Creating a decision tree multiclass classifier 
___________________________________________________________________________________________________________________________________________________________________

## 1- Data exploratory analysis

#### variables data distribution/shape

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/histograms.JPG)

#### How samples cluster together
##### all subtypes
![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/pca1.JPG)

##### uro vs basal

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/pca2.JPG)

##### How PCs contribute to the observed variance

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/scree_plot.JPG) 

##### What features contributed to each PCs

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/corrplot.JPG) 

___________________________________________________________________________________________________________________________________________________________________

## 2- Feature selection 

In order to identify which features can potentially be useful to classify samples into different subtypes a set of binary ROC analyses were performed. 


##### Lum vs. Bas AUC result

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/luminal_basal_auc.JPG)

##### Uro vs. Gu AUC result

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/uroGuAUC.JPG)


##### Gu vs. Bas AUC result

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/GU_Basal_AUC.JPG)

___________________________________________________________________________________________________________________________________________________________________
## 3- Creating decision tree multiclass classifiers

Considering AUC analysis,  these markers found to do a good job in classifying subtypes: 

"KRT14", "KRT5", "CDH3", "FOXA1", "GATA3", "PPARG", "RB1" , "CCND1", "CDKN2Ap16", "FGFR3", "TP63"

##### a look at the dataset


![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/dataset_1.JPG)

##### correlation matrix

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/corrmat.JPG)

##### Model result
```python
#___________________Model performance on training set_________________________#
## sample encoding : 0 = Basal, 1= Uro, 2 = GU

#Accuracy: 0.92
#Confusion Matrix: 

# [[ 38   5   0]
# [  0 143   3]
# [  0  13  68]]

#              precision    recall  f1-score   support
#
#           0       1.00      0.88      0.94        43
#           1       0.89      0.98      0.93       146
#           2       0.96      0.84      0.89        81
#
#    accuracy                           0.92       270
#   macro avg       0.95      0.90      0.92       270
#weighted avg       0.93      0.92      0.92       270
```


```python
#___________________Model performance on test set_________________________#
#Accuracy: 0.89
#Confusion Matrix: 
# [[11  3  0]
# [ 0 49  0]
# [ 2  5 21]]
#              precision    recall  f1-score   support
#
#           0       0.85      0.79      0.81        14
#           1       0.86      1.00      0.92        49
#           2       1.00      0.75      0.86        28
#
#    accuracy                           0.89        91
#   macro avg       0.90      0.85      0.87        91
# weighted avg       0.90      0.89      0.89        91
```

````python
# Tuning the depth of the trees,
# accuracy for the model from depth 1 to 5:
#[0.7912087912087912,
# 0.8681318681318682,
# 0.8461538461538461,
# 0.8461538461538461,
# 0.8461538461538461]
# best accuracy achived with depth 2! 
````
#### feature importance:

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/feature_importance.JPG)

#### tree visualization [to be added]


*
