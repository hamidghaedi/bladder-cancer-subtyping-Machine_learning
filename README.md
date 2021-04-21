# bladder-cancer-tumour-cell-phenotype-classifiers

#### Decision tree classifier to predict major subtypes of muscle invasive bladder cancer using IHC data

#### Steps:
##### 1- Data exploratory analysis 
##### 2- Feature selection by ROC analysis 
##### 3- Creating a decision tree multiclass classifier 

_____________________________________________________________________________________________________________________________________

## 1- Data exploratory analysis

#### variables data distribution/shape. 
Those markers which shows distinct pattern for each subgroups might worth to have them as feature in making the classifier. 

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/histograms.JPG)

#### How samples cluster together

All sample                 |  Uro vs. BaSq
:-------------------------:|:-------------------------:
<img src="https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/pca1.JPG" width="800" height="500">  |  <img src="https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/pca2.JPG" width="800" height="500">

___________________________________________________________________________________________________________________________________
## 2- Feature selection 

In order to identify which features can potentially be useful to classify samples into different subtypes a set of binary ROC analyses were performed. 


##### AUC result


Lum vs. Bas                |  Uro vs. Gu               |  Gu vs. Bas
:-------------------------:|:-------------------------:|:-------------------------:
![](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/luminal_basal_auc.JPG)  |  ![](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/uroGuAUC.JPG) | ![](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/GU_Basal_AUC.JPG) 

________________________________________________________________________________________________________________________________
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

#### Finding the optimal tree depth by cross validation
Best depth based on the cross validation accuracy apears to be three. 
<img src="https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/tree%20depth.JPG" width="800" height="400">

####  Multiclass AUC analysis:
 To assess the performance of the classifier another measure that I wanted to use is AUC. 

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/AUC_1.JPG)



### Classifier with reduced features

In order to make classifiers to work on NMIBC dataset like Sjödahl et al. Int. J. Cancer: 146, 2636–2647 (2019), the classifiers should be modified in the way that they work with the common features between the two dataset.Since both dataset are coming from the same laboratory the classifres should work on NMIBC dataset also without being affected by serious biases like those arising from data preparation.

Fortunately important features are common (7 out of 10) between two datasets.

Subtypes of samples in the NMIBC data set is as follows:


Ba/Sq  |  GU    | Mes | Sc/NE |  Uro 
:-----:|:------:|:---:|:-----:|:---:
8      | 55     |  7  | 6     | 192   

Common features between NMIBC dataset and the selected markers of MIBC:

"CCND1","CDKN2Ap16","FGFR3","FOXA1","GATA3","KRT14","KRT5","RB1"

```python
#___________________Model performance on the training set_________________________#

subtypes codes
  Code  Subtype
------  ---------
     0  BaSq
     1  Uro
     2  GU

Classification results on training set
Accuracy: 0.92
Confusion Matrix: 
 [[ 37   5   1]
 [  0 143   3]
 [  0  13  68]]
              precision    recall  f1-score   support

           0       1.00      0.86      0.92        43
           1       0.89      0.98      0.93       146
           2       0.94      0.84      0.89        81

    accuracy                           0.92       270
   macro avg       0.94      0.89      0.92       270
weighted avg       0.92      0.92      0.92       270

#___________________Model performance on the test set_________________________#

Classification results on test set
Accuracy: 0.90
Confusion Matrix: 
 [[10  3  1]
 [ 0 49  0]
 [ 0  5 23]]
              precision    recall  f1-score   support

           0       1.00      0.71      0.83        14
           1       0.86      1.00      0.92        49
           2       0.96      0.82      0.88        28

    accuracy                           0.90        91
   macro avg       0.94      0.85      0.88        91
weighted avg       0.91      0.90      0.90        91
```

### Reduced feature decision tree visualization 
![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/decision%20tree2.JPG)

####  Multiclass AUC analysis for the reduced feature decision tree

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/AUC_2.JPG)

### The model performance of the NMIBC data

```python


Classification results on validation set
Accuracy: 0.90
Confusion Matrix: 
 [[  6   1   1]
 [  0 186   6]
 [  0  17  38]]
              precision    recall  f1-score   support

           0       1.00      0.75      0.86         8
           1       0.91      0.97      0.94       192
           2       0.84      0.69      0.76        55

    accuracy                           0.90       255
   macro avg       0.92      0.80      0.85       255
weighted avg       0.90      0.90      0.90       255

```
