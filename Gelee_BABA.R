# Test de specialité Statistiques 2017
# BABA HALIMA IENM2
# 03 juillet 2017

rm(list=ls())
graphics.off()

# Chargement des librairies
library(verification)
library(MASS)
library(class)
library(randomForest)
library(chron)
library(car)
library(ipred)
library(rpart)
library(gbm)
library(e1071)
library(boot)
library(ROCR)

# Fonction calcul du pss
pss=function(prev,obs) {
	table=table(prev,obs)
	H=table[2,2]/(table[2,2]+table[1,2])
	F=table[2,1]/(table[2,1]+table[1,1])
	TSG=sum(diag(table))/sum(table)
	pss=H-F
	print(paste("H =",round(H,3)))
	print(paste("F =",round(F,3)))
	print(paste("TSG=",round(TSG,3)))		
	print(paste("pss=",round(pss,3))) 
	resultat=c(H,F,TSG,pss)
	return(resultat)}

# Fonction qui réalise le graphique des résidus
plot.res=function(x,y,titre="titre") {
	plot(x,y,col="blue",xlim=c(0,1),ylim=c(-1,1),
	ylab="Résidus",xlab="Valeurs predites",main=titre)
	abline(h=0,col="green")}

# Fonction qui réalise le graphique des valeurs ajustées
plot.fit=function(x,y,titre="titre") {
	plot(x,y,col="blue",xlim=c(0,1),ylim=c(-0.1,1),
	ylab="Ajustées",xlab="Observations P60",main=titre)
	abline(0,1,col="green")
	abline(150,0,col="pink")
	abline(v=150,col="pink")}

##############################################
#          1. Chargement des donnees         #
##############################################


data = read.table("DataGelee.txt",sep=";",header=T)
names(data)
summary(data)

# Suppression de la variable "DATE"
data = data[,-1] 

# Division par 10 des température
data[,"T"]=data[,"T"]/10
data[,"TD"]=data[,"TD"]/10

# Transformation des variables GLBp et GLBo en variables catégoriques
data[,"GLBp"]=as.factor(data[,"GLBp"])
data[,"GLBo"]=as.factor(data[,"GLBo"])

# Apreçu des donnees
x11()
par(mfrow=c(4,2))
hist(data[,"N"]);hist(data[,"T"])
hist(data[,"TD"]);hist(data[,"U"])
hist(data[,"FF"])

# Etude de la correlation des vairables
cor(data)
x11()
pairs(data,gap=0,pch='.')

# Puisque T et TD sont corrélées, on peut supprimer l'une des deux variables
data=data[,-4]


##############################################
#          2. Estimation du meilleur         #
#    modele statistique pour l'occurence     #
##############################################


# Extraction d'un échantillon aléatoire
set.seed(222)                   # initialisation du générateur
test.ratio=.2                   # part de l'échantillon test (20%)
npop=nrow(data)                 # nombre de lignes dans les données
nvar=ncol(data)                 # nombre de colonnes
ntest=ceiling(npop*test.ratio)  # taille de l'échantillon test

testi=sample(1:npop,ntest)   # indices (tirés aléatoirement) de l'échantillon test
appri=setdiff(1:npop,testi)  # indices complémentaires de l'échantillon d'apprentissage

# Construction des échantillons test et apprentissage
data_appr = data[appri,] 
data_test = data[testi,]



#### REGRESSION LOGISTIQUE SANS INTERACTION


glb = as.numeric(data_test[,"GLBo"])-1
glbr = as.numeric(data_appr[,"GLBo"])-1
glbp = as.numeric(data_test[,"GLBp"])-1
glbpr = as.numeric(data_appr[,"GLBp"])-1


### MODELE 1 : COMPLET


log.lm=glm(GLBo~.,data=data_appr,family=binomial)
summary(log.lm)

# Matrice de confusion de l’échantillon d’apprentissage
table(log.lm$fitted.values>0.5,data_appr[,"GLBo"])

# Courbe ROC
x11()
ROC=roc.plot(glbr,log.lm$fitted.values)
i_lm_r=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_lm_r=ROC$plot.data[i_lm_r,1,1]
seuil_lm_r

# Calcul des scores de Pierce et de Brier
pss(log.lm$fitted.values>seuil_lm_r,data_appr[,"GLBo"])
brier(glbr,log.lm$fitted.values)

# Diagramme d'attribut
verif=verify(glbr,log.lm$fitted.values)
x11()
attribute(verif,main="Diagramme d'attribut modele complet apprentissage \n Regression log sans interaction")

# Diagramme de vraisemblance
x11()
discrimination.plot(glbr,log.lm$fitted.values,main="Diagramme de vraisemblance modele complet apprentissage \n Regression log sans interaction")


### PREVISIONS


pred.log=predict(log.lm,newdata=data_test,type="response")

# Matrice de confusion pour la prévision
table(pred.log>0.5,data_test[,"GLBo"])

# Courbe ROC
x11()
ROC=roc.plot(glb,pred.log)
i_lm_t=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_lm_t=ROC$plot.data[i_lm_t,1,1]
seuil_lm_t

# Calcul des scores de Pierce et de Brier
pss(pred.log>seuil_lm_t,data_test[,"GLBo"])
brier(glb,pred.log)

# Diagramme d'attribut
verif=verify(glb,pred.log)
x11()
attribute(verif,main="Diagramme d'attribut modele complet test \n Regression log sans interaction")

# Diagramme de vraisemblance
x11()
discrimination.plot(glb,pred.log,main="Diagramme de vraisemblance modele complet test \n Regression log sans interaction")

# Calcul des scores de Pierce et de Brier pour GLBp
pss(glbp>seuil_lm_t,data_test[,"GLBo"])
brier(glb,glbp)


### MODELE 2 : Sélection de modèle par sélection de variables


# AIC backward
log.lm.AIC=stepAIC(log.lm,direction="backward")
summary(log.lm.AIC)

# AIC foreward
mod0=glm(GLBo~1,data=data_appr,family=binomial)
log.lm.AICF=stepAIC(mod0,GLBo~.,trace=TRUE,direction=c("forward"))
summary(log.lm.AICF)

# AIC stepwise
log.lm.AICS=stepAIC(log.lm,~.,trace=TRUE,direction=c("both"))
summary(log.lm.AICS)

# BIC stepwise
# k=log(npop) pour BIC au lieu de k=2 pour AIC
log.lm.BIC=stepAIC(log.lm,~.,trace=TRUE,direction=c("both"),k=log(npop))

# L'application des méthodes AIC et BIC ne change pas le nombre de predicteurs, on elime alors l'utilisation de ces modeles



#### REGRESSION LOGISTIQUE AVEC INTERACTION



### MODELE 1 : COMPLET


log.qm=glm(GLBo~(.)^2,data=data_appr,family=binomial)
summary(log.qm)

# Matrice de confusion de l’échantillon d’apprentissage
table(log.qm$fitted.values>0.5,data_appr[,"GLBo"])

# Courbe ROC
x11()
ROC=roc.plot(glbr,log.qm$fitted.values)
i_qm_r=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_qm_r=ROC$plot.data[i_qm_r,1,1]
seuil_qm_r

# Calcul des scores de Pierce et de Brier
pss(log.qm$fitted.values>seuil_qm_r,data_appr[,"GLBo"])
brier(glbr,log.qm$fitted.values)

# Diagramme d'attribut
verif=verify(glbr,log.qm$fitted.values)
x11()
attribute(verif,main="Diagramme d'attribut modele complet apprentissage \n Regression log avec interaction")

# Diagramme de vraisemblance
x11()
discrimination.plot(glbr,log.qm$fitted.values,main="Diagramme de vraisemblance modele complet apprentissage\n Regression log avec interaction")


### PREVISIONS


pred.log.qm=predict(log.qm,newdata=data_test,type="response")

# Matrice de confusion pour la prévision
table(pred.log.qm>0.5,data_test[,"GLBo"])

# Courbe ROC
x11()
ROC=roc.plot(glb,pred.log.qm)
i_qm_t=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_qm_t=ROC$plot.data[i_qm_t,1,1]
seuil_qm_t

# Calcul des scores de Pierce et de Brier
pss(pred.log.qm>seuil_qm_t,data_test[,"GLBo"])
brier(glb,pred.log.qm)

# Diagramme d'attribut
verif=verify(glb,pred.log.qm)
x11()
attribute(verif,main="Diagramme d'attribut modele complet test \n Regression log avec interaction")

# Diagramme de vraisemblance
x11()
discrimination.plot(glb,pred.log.qm,main="Diagramme de vraisemblance modele complet test \n Regression log avec interaction")

# Calcul des scores de Pierce et de Brier pour GLBp
pss(glbp>seuil_qm_t,data_test[,"GLBo"])
brier(glb,glbp)


### MODELE 2 : Sélection de modèle par sélection de variables


## AIC backward


log.qm.AIC=stepAIC(log.qm,direction="backward")
summary(log.qm.AIC)

# Matrice de confusion de l’échantillon d’apprentissage
table(log.qm.AIC$fitted.values>0.5,data_appr[,"GLBo"])

# Courbe ROC
x11()
ROC=roc.plot(glbr,log.qm.AIC$fitted.values)
i_qm_aic_r=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_qm_aic_r=ROC$plot.data[i_qm_aic_r,1,1]
seuil_qm_aic_r

# Calcul des scores de Pierce et de Brier
pss(log.qm.AIC$fitted.values>seuil_qm_aic_r,data_appr[,"GLBo"])
brier(glbr,log.qm.AIC$fitted.values)

# Diagramme d'attribut
verif=verify(glbr,log.qm.AIC$fitted.values)
x11()
attribute(verif,main="Diagramme d'attribut modele reduit apprentissage \n Regression log avec interaction")

# Diagramme de vraisemblance
x11()
discrimination.plot(glbr,log.qm.AIC$fitted.values,main="Diagramme de vraisemblance modele reduit apprentissage\n Regression log avec interaction")


## PREVISION


pred.qm.AIC=predict(log.qm.AIC,newdata=data_test,type="response")

# Matrice de confusion pour la prévision
table(pred.qm.AIC>0.5,data_test[,"GLBo"])

# Courbe ROC
x11()
ROC=roc.plot(glb,pred.qm.AIC)
i_qm_aic_t=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_qm_aic_t=ROC$plot.data[i_qm_aic_t,1,1]
seuil_qm_aic_t

# Calcul des scores de Pierce et de Brier
pss(pred.qm.AIC>seuil_qm_aic_t,data_test[,"GLBo"])
brier(glb,pred.qm.AIC)

# Diagramme d'attribut
verif=verify(glb,pred.qm.AIC)
x11()
attribute(verif,main="Diagramme d'attribut modele complet test \n Regression log avec interaction")

# Diagramme de vraisemblance
x11()
discrimination.plot(glb,pred.qm.AIC,main="Diagramme de vraisemblance modele complet test \n Regression log avec interaction")

# Calcul des scores de Pierce et de Brier pour GLBp
pss(glbp>seuil_qm_aic_t,data_test[,"GLBo"])
brier(glb,glbp)



#### ANALYSE DISCRIMINANTE SANS INTERACTION



### MODELE

AD_lin=lda(GLBo~.,data=data_appr,CV=T)

# Matrice de confusion de l’échantillon d’apprentissage
table(AD_lin>0.5,data_appr[,"GLBo"])

# Courbe ROC
x11()
ROC=roc.plot(glbr,AD_lin$posterior[,2])
i_lda_r=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_lda_r=ROC$plot.data[i_lda_r,1,1]
seuil_lda_r

# Calcul des scores de Pierce et de Brier
pss(AD_lin$posterior[,2]>seuil_lda_r,data_appr[,"GLBo"])
brier(glbr,AD_lin$posterior[,2])

# Diagramme d'attribut
verif=verify(glbr,AD_lin$posterior[,2])
x11()
attribute(verif,main="Diagramme d'attribut apprentissage \n Analyse discriminante lineaire")

# Diagramme de vraisemblance
x11()
discrimination.plot(glbr,AD_lin$posterior[,2],main="Diagramme de vraisemblance apprentissage\n Analyse discriminante lineaire")


### PREVISIONS


AD_lin=lda(GLBo~.,data=data_appr)
pred.lda=predict(AD_lin,newdata=data_test)$posterior[,2]

# Courbe ROC
x11()
ROC=roc.plot(glb,pred.lda)
i_lda_t=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_lda_t=ROC$plot.data[i_lda_t,1,1]
seuil_lda_t

# Matrice de confusion pour la prévision
table(pred.lda>seuil_lda_t,data_test[,"GLBo"])

# Calcul des scores de Pierce et de Brier
pss(pred.lda>seuil_lda_t,data_test[,"GLBo"])
brier(glb,pred.lda)

# Diagramme d'attribut
verif=verify(glb,pred.lda)
x11()
attribute(verif,main="Diagramme d'attribut test \n Analyse discriminante lineaire")

# Diagramme de vraisemblance
x11()
discrimination.plot(glb,pred.lda,main="Diagramme de vraisemblance test \n Analyse discriminante lineaire")

# Calcul des scores de Pierce et de Brier pour GLBp
pss(glbp>seuil_lda_t,data_test[,"GLBo"])
brier(glb,glbp)



#### ANALYSE DISCRIMINANTE AVEC INTERACTION



### MODELE


AD_q=qda(GLBo~.,data=data_appr,CV=T)

# Matrice de confusion de l’échantillon d’apprentissage
table(AD_q>0.5,data_appr[,"GLBo"])

# Courbe ROC
x11()
ROC=roc.plot(glbr,AD_q$posterior[,2])
i_qda_r=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_qda_r=ROC$plot.data[i_qda_r,1,1]
seuil_qda_r

# Calcul des scores de Pierce et de Brier
pss(AD_q$posterior[,2]>seuil_qda_r,data_appr[,"GLBo"])
brier(glbr,AD_q$posterior[,2])

# Diagramme d'attribut
verif=verify(glbr,AD_q$posterior[,2])
x11()
attribute(verif,main="Diagramme d'attribut apprentissage \n Analyse discriminante lineaire")

# Diagramme de vraisemblance
x11()
discrimination.plot(glbr,AD_q$posterior[,2],main="Diagramme de vraisemblance apprentissage\n Analyse discriminante lineaire")


### PREVISIONS


AD_q=qda(GLBo~.,data=data_appr)
pred.qda=predict(AD_q,newdata=data_test)$posterior[,2]

# Matrice de confusion pour la prévision
table(pred.qda>seuil_qda_t,data_test[,"GLBo"])

# Courbe ROC
x11()
ROC=roc.plot(glb,pred.qda)
i_qda_t=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_qda_t=ROC$plot.data[i_qda_t,1,1]
seuil_qda_t

# Calcul des scores de Pierce et de Brier
pss(pred.qda>seuil_qda_t,data_test[,"GLBo"])
brier(glb,pred.qda)

# Diagramme d'attribut
verif=verify(glb,pred.qda)
x11()
attribute(verif,main="Diagramme d'attribut test \n Analyse discriminante quadratique")

# Diagramme de vraisemblance
x11()
discrimination.plot(glb,pred.qda,main="Diagramme de vraisemblance test \n Analyse discriminante quadratique")

# Calcul des scores de Pierce et de Brier pour GLBp
pss(glbp>seuil_qda_t,data_test[,"GLBo"])
brier(glb,glbp)



#### FORETS ALEATOIRES



foret=randomForest(GLBo~.,data=data_appr,xtest=data_test[,-7],ytest=data_test[,"GLBo"],ntree=500,do.trace=50,importance=TRUE)

# Matrice de confusion de l’échantillon d’apprentissage
table(foret$votes[,2]>0.5,data_appr[,"GLBo"])

# Courbe ROC
x11()
ROC=roc.plot(glbr,foret$votes[,2])
i_foret_r=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_foret_r=ROC$plot.data[i_foret_r,1,1]
seuil_foret_r

# Calcul des scores de Pierce et de Brier
pss(foret$votes[,2]>seuil_foret_r,data_appr[,"GLBo"])
brier(glbr,foret$votes[,2])

# Diagramme d'attribut
verif=verify(glbr,foret$votes[,2])
x11()
attribute(verif,main="Diagramme d'attribut apprentissage \n Foret aleatoire")

# Diagramme de vraisemblance
x11()
discrimination.plot(glbr,foret$votes[,2],main="Diagramme de vraisemblance apprentissage\n Foret aleatoire")


### PREVISIONS


pred.rfq=foret$test$votes[,2]

# Matrice de confusion pour la prévision du dépassement de seuil
table(pred.rfq,data_test$GLBo)

# Courbe ROC
ROC_rfq=roc.plot(glb,pred.rfq)
i_foret_t=which.max(ROC_rfq$plot.data[,2,1]-ROC_rfq$plot.data[,3,1])
seuil_foret_t=ROC_rfq$plot.data[i_foret_t,1,1]
seuil_foret_t

# Matrice de confusion pour la prévision du dépassement de seuil
table(pred.rfq>seuil_foret_t,data_test[,"GLBo"])

# Calcul des scores de Pierce et de Brier
pss(pred.rfq>seuil_foret_t,data_test[,"GLBo"])
brier(glb,pred.rfq)

# Diagramme d'attribut
verif=verify(glb,pred.rfq)
x11()
attribute(verif,main="Diagramme d'attribut test \n Foret aleatoire")

# Diagramme de vraisemblance
x11()
discrimination.plot(glb,pred.rfq,main="Diagramme de vraisemblance test \n Foret aleatoire")

# Calcul des scores de Pierce et de Brier pour GLBp
pss(glbp>seuil_foret_t,data_test[,"GLBo"])
brier(glb,glbp)

# Conclusion : par comparaison des differents scores obtenus sur les fichiers test et apprentissage pour chaque modele, il semblerait que la methode de foret aleatoire soit la plus pertinente (PSS : 0.922 le plus élevé et BS : 0.024 le plus faible). De plus la courbe d'attribut est celle qui représente le meilleur compromis entre fiabilité et résolution. La courbe de discrimination montre bien que le modèle prévient bien l'occurence de non gel et un peu moins bien l'occurence de gel. Tous les modèles sont meilleurs que les donnees de prevision GLBp.


##############################################
#          3. Etude variabilite              #
#                 H,F,PSS                    #
##############################################


pss_appr_lm=0
pss_test_lm=0
pss_appr_qm=0
pss_test_qm=0
pss_appr_aic=0
pss_test_aic=0
pss_appr_lda=0
pss_test_lda=0
pss_appr_qda=0
pss_test_qda=0
pss_appr_foret=0
pss_test_foret=0

for(i in 100:120){
# Extraction d'un échantillon aléatoire
set.seed(i)                   # initialisation du générateur
test.ratio=.2                   # part de l'échantillon test (20%)
npop=nrow(data)                 # nombre de lignes dans les données
nvar=ncol(data)                 # nombre de colonnes
ntest=ceiling(npop*test.ratio)  # taille de l'échantillon test

testi=sample(1:npop,ntest)   # indices (tirés aléatoirement) de l'échantillon test
appri=setdiff(1:npop,testi)  # indices complémentaires de l'échantillon d'apprentissage

# Construction des échantillons test et apprentissage
data_appr = data[appri,] 
data_test = data[testi,]

# donnees converties en numerique

glb = as.numeric(data_test[,"GLBo"])-1
glbr = as.numeric(data_appr[,"GLBo"])-1
glbp = as.numeric(data_test[,"GLBp"])-1
glbpr = as.numeric(data_appr[,"GLBp"])-1

# modele regression logistique sans interaction

log.lm=glm(GLBo~.,data=data_appr,family=binomial)
ROC=roc.plot(glbr,log.lm$fitted.values)
i_lm_r=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_lm_r=ROC$plot.data[i_lm_r,1,1]
seuil_lm_r
pss_appr_lm=rbind(pss_appr_lm,pss(log.lm$fitted.values>seuil_lm_r,data_appr[,"GLBo"]))

pred.log=predict(log.lm,newdata=data_test,type="response")
ROC=roc.plot(glb,pred.log)
i_lm_t=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_lm_t=ROC$plot.data[i_lm_t,1,1]
seuil_lm_t
pss_test_lm=rbind(pss_test_lm,pss(pred.log>seuil_lm_t,data_test[,"GLBo"]))

log.qm=glm(GLBo~(.)^2,data=data_appr,family=binomial)
ROC=roc.plot(glbr,log.qm$fitted.values)
i_qm_r=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_qm_r=ROC$plot.data[i_qm_r,1,1]
seuil_qm_r
pss_appr_qm=rbind(pss_appr_qm,pss(log.qm$fitted.values>seuil_qm_r,data_appr[,"GLBo"]))

pred.log.qm=predict(log.qm,newdata=data_test,type="response")
ROC=roc.plot(glb,pred.log.qm)
i_qm_t=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_qm_t=ROC$plot.data[i_qm_t,1,1]
seuil_qm_t
pss_test_qm=rbind(pss_test_qm,pss(pred.log.qm>seuil_qm_t,data_test[,"GLBo"]))

log.qm.AIC=stepAIC(log.qm,direction="backward")
ROC=roc.plot(glbr,log.qm.AIC$fitted.values)
i_qm_aic_r=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_qm_aic_r=ROC$plot.data[i_qm_aic_r,1,1]
seuil_qm_aic_r
pss_appr_aic=rbind(pss_appr_aic,pss(log.qm.AIC$fitted.values>seuil_qm_aic_r,data_appr[,"GLBo"]))

pred.qm.AIC=predict(log.qm.AIC,newdata=data_test,type="response")
ROC=roc.plot(glb,pred.qm.AIC)
i_qm_aic_t=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_qm_aic_t=ROC$plot.data[i_qm_aic_t,1,1]
seuil_qm_aic_t
pss_test_aic=rbind(pss_test_aic,pss(pred.qm.AIC>seuil_qm_aic_t,data_test[,"GLBo"]))

AD_lin=lda(GLBo~.,data=data_appr,CV=T)
ROC=roc.plot(glbr,AD_lin$posterior[,2])
i_lda_r=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_lda_r=ROC$plot.data[i_lda_r,1,1]
seuil_lda_r
pss_appr_lda=rbind(pss_appr_lda,pss(AD_lin$posterior[,2]>seuil_lda_r,data_appr[,"GLBo"]))

AD_lin=lda(GLBo~.,data=data_appr)
pred.lda=predict(AD_lin,newdata=data_test)$posterior[,2]
ROC=roc.plot(glb,pred.lda)
i_lda_t=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_lda_t=ROC$plot.data[i_lda_t,1,1]
seuil_lda_t
pss_test_lda=rbind(pss_test_lda,pss(pred.lda>seuil_lda_t,data_test[,"GLBo"]))

AD_q=qda(GLBo~.,data=data_appr,CV=T)
ROC=roc.plot(glbr,AD_q$posterior[,2])
i_qda_r=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_qda_r=ROC$plot.data[i_qda_r,1,1]
seuil_qda_r
pss_appr_qda=rbind(pss_appr_qda,pss(AD_q$posterior[,2]>seuil_qda_r,data_appr[,"GLBo"]))

AD_q=qda(GLBo~.,data=data_appr)
pred.qda=predict(AD_q,newdata=data_test)$posterior[,2]
ROC=roc.plot(glb,pred.qda)
i_qda_t=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_qda_t=ROC$plot.data[i_qda_t,1,1]
seuil_qda_t
pss_test_qda=rbind(pss_test_qda,pss(pred.qda>seuil_qda_t,data_test[,"GLBo"]))

foret=randomForest(GLBo~.,data=data_appr,xtest=data_test[,-7],ytest=data_test[,"GLBo"],ntree=500,do.trace=50,importance=TRUE)
ROC=roc.plot(glbr,foret$votes[,2])
i_foret_r=which.max(ROC$plot.data[,2,1]-ROC$plot.data[,3,1])
seuil_foret_r=ROC$plot.data[i_foret_r,1,1]
seuil_foret_r
pss_appr_foret=rbind(pss_appr_foret,pss(foret$votes[,2]>seuil_foret_r,data_appr[,"GLBo"]))

pred.rfq=foret$test$votes[,2]
ROC_rfq=roc.plot(glb,pred.rfq)
i_foret_t=which.max(ROC_rfq$plot.data[,2,1]-ROC_rfq$plot.data[,3,1])
seuil_foret_t=ROC_rfq$plot.data[i_foret_t,1,1]
seuil_foret_t
pss_test_foret=rbind(pss_test_foret,pss(pred.rfq>seuil_foret_t,data_test[,"GLBo"]))
}

pss_appr_lm=pss_appr_lm[-1,]
pss_test_lm=pss_test_lm[-1,]
pss_appr_qm=pss_appr_qm[-1,]
pss_test_qm=pss_test_qm[-1,]
pss_appr_aic=pss_appr_aic[-1,]
pss_test_aic=pss_test_aic[-1,]
pss_appr_lda=pss_appr_lda[-1,]
pss_test_lda=pss_test_lda[-1,]
pss_appr_qda=pss_appr_qda[-1,]
pss_test_qda=pss_test_qda[-1,]
pss_appr_foret=pss_appr_foret[-1,]
pss_test_foret=pss_test_foret[-1,]

list_pss=c("appr_lm","test_lm","appr_qm","test_qm","appr_aic","test_aic","appr_lda","test_lda","appr_qda","test_qda","appr_foret","test_foret")
h_tot=cbind(pss_appr_lm[,1],pss_test_lm[,1],pss_appr_qm[,1],pss_test_qm[,1],pss_appr_aic[,1],pss_test_aic[,1],pss_appr_lda[,1],pss_test_lda[,1],pss_appr_qda[,1],pss_test_qda[,1],pss_appr_foret[,1],pss_test_foret[,1])
f_tot=cbind(pss_appr_lm[,2],pss_test_lm[,2],pss_appr_qm[,2],pss_test_qm[,2],pss_appr_aic[,2],pss_test_aic[,2],pss_appr_lda[,2],pss_test_lda[,2],pss_appr_qda[,2],pss_test_qda[,2],pss_appr_foret[,2],pss_test_foret[,2])
pss_tot=cbind(pss_appr_lm[,4],pss_test_lm[,4],pss_appr_qm[,4],pss_test_qm[,4],pss_appr_aic[,4],pss_test_aic[,4],pss_appr_lda[,4],pss_test_lda[,4],pss_appr_qda[,4],pss_test_qda[,4],pss_appr_foret[,4],pss_test_foret[,4])

x11()
boxplot(t(pss_tot)~list_pss,main="boxplot pss")
x11()
boxplot(t(h_tot)~list_pss,main="boxplot hit rate")
x11()
boxplot(t(f_tot)~list_pss,main="boxplot false alarm")
