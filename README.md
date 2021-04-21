# bladder-cancer-tumour-cell-phenotype-classifiers

#### Decision tree classifier to predict major subtypes of muscle invasive bladder cancer using IHC data

#### Steps:
##### 1- Data exploratory analysis 
##### 2- Feature selection by ROC analysis 
##### 3- Creating a decision tree multiclass classifier 

___________________________________________________________________________________________________________________________________________________________________

## 1- Data exploratory analysis

#### variables data distribution/shape. 
Those markers which shows distinct pattern for each subgroups might worth to have them as feature in making the classifier. 

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/histograms.JPG)

#### How samples cluster together

All sample                 |  Uro vs. BaSq
:-------------------------:|:-------------------------:
<img src="https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/pca1.JPG" width="800" height="500">  |  <img src="https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/pca2.JPG" width="800" height="500">

##### How PCs contribute to the observed variance

<img src="https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/scree_plot.JPG" width="700" height="500">

##### What features contributed to each PCs

<img src="https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/corrplot.JPG" width="300" height="800">
___________________________________________________________________________________________________________________________________________________________________

## 2- Feature selection 

In order to identify which features can potentially be useful to classify samples into different subtypes a set of binary ROC analyses were performed. 


##### AUC result


Lum vs. Bas                |  Uro vs. Gu               |  Gu vs. Bas
:-------------------------:|:-------------------------:|:-------------------------:
![](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/luminal_basal_auc.JPG)  |  ![](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/uroGuAUC.JPG) | ![](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/GU_Basal_AUC.JPG) 

__________________________________________________________________________________________________________________________________________________________________
## 3- Creating decision tree multiclass classifiers

Considering AUC analysis,  these markers found to do a good job in classifying subtypes: 

"KRT14", "KRT5", "CDH3", "FOXA1", "GATA3", "PPARG", "RB1" , "CCND1", "CDKN2Ap16", "FGFR3", "TP63"

##### A look at the dataset


![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/dataset_1.JPG)

##### Correlation matrix
This would be helpful to reduce the feature number by keeping one feature from a set of highly correlated (r >= .8) features.

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/corrmat.JPG)

##### Model result
```python
#___________________Model performance on training set_________________________#
subtypes codes
  Code  Subtype
------  ---------
     0  BaSq
     1  Uro
     2  GU

Classification results on training set
Accuracy: 0.93
Confusion Matrix: 
 [[ 37   6   0]
 [  0 146   0]
 [  0  14  67]]
              precision    recall  f1-score   support

           0       1.00      0.86      0.92        43
           1       0.88      1.00      0.94       146
           2       1.00      0.83      0.91        81

    accuracy                           0.93       270
   macro avg       0.96      0.90      0.92       270
weighted avg       0.93      0.93      0.93       270
```


```python

#___________________Model performance on test set_________________________#
Classification results on test set
Accuracy: 0.86
Confusion Matrix: 
 [[10  3  1]
 [ 0 49  0]
 [ 1  8 19]]
              precision    recall  f1-score   support

           0       0.91      0.71      0.80        14
           1       0.82      1.00      0.90        49
           2       0.95      0.68      0.79        28

    accuracy                           0.86        91
   macro avg       0.89      0.80      0.83        91
weighted avg       0.87      0.86      0.85        91

```
### Visualization of the decision tree
<img src="https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/decision_tree_1.JPG">

#### Finding the otimal tree depth by cross validation
<img src="https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/tree%20depth.JPG" width="400" height="800">

#### feature importance:

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/feature_importance.JPG)
