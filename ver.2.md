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

![alt text](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/00_marker_val_dist.png)

#### How samples cluster together


<img src="https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/01_pca.png" width="800" height="500"> 

___________________________________________________________________________________________________________________________________
## 2- Feature selection 

In order to identify which features can potentially be useful to classify samples into different subtypes a set of binary ROC analyses were performed. 


##### AUC result


![](https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/features.PNG) 

________________________________________________________________________________________________________________________________
## 3- Creating decision tree multiclass classifiers

Considering AUC analysis,  these markers found to do a good job (AUC >= 80%) in classifying subtypes, were included.

##### Finding tree depth by K-Fold cross validation
```python
Depth: 1 Accuracy: 0.706
Depth: 2 Accuracy: 0.836
Depth: 3 Accuracy: 0.830
Depth: 4 Accuracy: 0.858
Depth: 5 Accuracy: 0.836
Depth: 6 Accuracy: 0.836
Depth: 7 Accuracy: 0.842
Depth: 8 Accuracy: 0.853
Depth: 9 Accuracy: 0.853
```


##### Model trained by 10-K stratified cross validation
```python
Mean score: 0.856 (+/-0.007)

#___________________Model performance on test set_________________________#
subtypes codes
  Code  Subtype
------  ---------
     0  BaSq
     1  Uro
     2  GU


Classification results on test set
Accuracy: 0.82
Confusion Matrix: 
 [[ 8  2  0]
 [ 1 20  0]
 [ 1  4  9]]
              precision    recall  f1-score   support

           0       0.80      0.80      0.80        10
           1       0.77      0.95      0.85        21
           2       1.00      0.64      0.78        14

    accuracy                           0.82        45
   macro avg       0.86      0.80      0.81        45
weighted avg       0.85      0.82      0.82        45


```
### Visualization of the decision tree
<img src="https://raw.githubusercontent.com/hamidghaedi/bladder-cancer-tumour-cell-phenotype-classification/main/figs/second_tree.PNG">


