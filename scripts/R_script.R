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

```R
# to see contribution of each variables to PCs
fviz_eig(new.pca, addlabels = TRUE, ylim = c(0, 30))
#
```

```R
# correlation plot to see how features contributed to PCs
var <- get_pca_var(new.pca)
corrplot(var$cos2, is.corr=FALSE)
```
#__________________________________________________________________________
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
