#*****************************************************************
#Modèles appliqués au jeu de données diabète après imputation mice
#*****************************************************************


library(data.table)
library(tidyverse)
library(gbm)
library(e1071)
library(glmnet)
library(randomForest)
library(plyr)
library(caret)




#===============================================================================
#SANS INTERACTION
#===============================================================================

nhanes <- read.csv("nhanes_dia_avant_transco.csv",header=TRUE,sep=",",dec=".")
str(nhanes)
nhanes[,1]<-NULL
dim(nhanes)#5221x85
str(nhanes)

nhanes$SEQN<-factor(nhanes$SEQN)
nhanes$RIAGENDR_demo<-factor(nhanes$RIAGENDR_demo)
nhanes$BPQ020_bpq<-as.factor(nhanes$BPQ020_bpq)
nhanes$BPQ080_bpq<-as.factor(nhanes$BPQ080_bpq)
nhanes$MCQ080_mcq<-as.factor(nhanes$MCQ080_mcq)
nhanes$MCQ160A_mcq<-as.factor(nhanes$MCQ160A_mcq)
nhanes$MCQ160B_mcq<-as.factor(nhanes$MCQ160B_mcq)
nhanes$MCQ160C_mcq<-as.factor(nhanes$MCQ160C_mcq)
nhanes$MCQ160D_mcq<-as.factor(nhanes$MCQ160D_mcq)
nhanes$MCQ160E_mcq<-as.factor(nhanes$MCQ160E_mcq)
nhanes$MCQ160F_mcq<-as.factor(nhanes$MCQ160F_mcq)
nhanes$MCQ160G_mcq<-as.factor(nhanes$MCQ160G_mcq)
nhanes$MCQ160M_mcq<-as.factor(nhanes$MCQ160M_mcq)
nhanes$MCQ160N_mcq<-as.factor(nhanes$MCQ160N_mcq)
nhanes$SLQ050_slq<-as.factor(nhanes$SLQ050_slq)
nhanes$HEQ010_heq<-as.factor(nhanes$HEQ010_heq)
nhanes$HEQ030_heq<-as.factor(nhanes$HEQ030_heq)


str(nhanes)

#Regression logistique avec les 84 variables (AIC 3634)
#------------------------------------------------------
set.seed(1234)
reglog<-glm(DIQ010_diq~.,data=nhanes[,-1],family="binomial")
summary(reglog)




#Regression logistique améliorée avec le step
#--------------------------------------------


null=glm(DIQ010_diq~1, data=nhanes[,-1],family=binomial)
null
full=glm(DIQ010_diq~., data=nhanes[,-1],family=binomial)
full
summary(full)

stepreglog=step(null, scope=list(lower=null, upper=full), direction="forward")
summary(stepreglog)
names(stepreglog)
stepreglog$formula




#Mise en facteur du Y pour la forêt
nhanesfac<-nhanes
nhanesfac$DIQ010_diq<-factor(nhanesfac$DIQ010_diq)




#validation croisee avec le dataset complet
#------------------------------------------
#
bloc<-4
set.seed(1234)
ind <- sample(1:nrow(nhanes)%%bloc+1)
RES<- data.frame(Y=nhanes$DIQ010_diq,
                 logcomp=0,
                 logstep=0,
                 logred15=0,
                 ridge=0,
                 lasso=0,
                 elastic=0,
                 foret=0,
                 adaboost=0,
                 logitboost=0,
                 svmlin=0,
                 svmrad=0)

XX<-model.matrix(nhanes$DIQ010_diq~.,data=nhanes[,-1])

LAridge=list()
LAlasso=list()
LAelas=list()

foreach (i = 1:bloc, .packages = c("gbm","glmnet","randomForest","e1071")) %dopar% {
  XXA <- XX[ind!=i,]
  YYA <- as.matrix(nhanes[ind!=i,"DIQ010_diq"])
  print(i)
  #1-logistique complète
  #-------------------
  print("log")
  mod <- glm(DIQ010_diq~.,data=nhanes[ind!=i,-1],family="binomial")
  RES[ind==i,"logcomp"] <- predict(mod,nhanes[ind==i,-1],type="response")
  
  # #2-logistique avec step
  # #--------------------
  print("logstep")
  RES[ind==i,"logstep"] <- predict(stepreglog,nhanes[ind==i,-1],type="response")
  
  
  #3-logistique réduite 15
  #-----------------------
  print("logred15")
  mod <- glm(DIQ010_diq~RIDAGEYR_demo +DR1TSUGR_dr1tot+BPQ080_bpq+MCQ080_mcq+BPQ020_bpq+DR1TALCO_dr1tot+DR1TMOIS_dr1tot+
               RIAGENDR_demo+DR1.320Z_dr1tot+DR1TCHOL_dr1tot+DR1TPROT_dr1tot+INDFMPIR_demo+DR1TIRON_dr1tot+Var_TENSIONDI+DR1TCAFF_dr1tot,
             data=nhanes[ind!=i,-1],family="binomial")
  RES[ind==i,"logred15"] <- predict(mod,nhanes[ind==i,-1],type="response")
  
  
  #4-ridge
  #------
  print("ridge")
  tmp <- cv.glmnet(XXA,YYA,alpha=0,family="binomial")
  LAridge <- c(LAridge,tmp$lambda.min)
  mod <- glmnet(XXA,YYA,alpha=0,lambda=tmp$lambda.min,family="binomial")
  RES[ind==i,"ridge"] <- predict(mod,newx=XX[ind==i,],type="response")
  
  #5-lasso
  #-----
  print("lasso")
  tmp <- cv.glmnet(XXA,YYA,alpha=1,family="binomial")
  LAlasso <- c(LAridge,tmp$lambda.min)
  mod <- glmnet(XXA,YYA,alpha=0,lambda=tmp$lambda.min,family="binomial")
  RES[ind==i,"lasso"] <- predict(mod,newx=XX[ind==i,],type="response")
  
  #6-elastic
  #-------
  print("elastic")
  tmp <- cv.glmnet(XXA,YYA,alpha=0.5,family="binomial")
  LAelas <- c(LAridge,tmp$lambda.min)
  mod <- glmnet(XXA,YYA,alpha=0,lambda=tmp$lambda.min,family="binomial")
  RES[ind==i,"elastic"] <- predict(mod,newx=XX[ind==i,],type="response")
  
  #7-foret
  #-----
  print("foret")
  mod <- randomForest(DIQ010_diq~.,data=nhanesfac[ind!=i,-1])
  RES[ind==i,"foret"] <- predict(mod,nhanesfac[ind==i,-1],type="prob")[,2]
  
  #8-adaboost
  #--------
  print("adaboost")
  tmp <- gbm(DIQ010_diq~.,data = nhanes[ind!=i,-1], distribution = "adaboost", interaction.depth = 2,
             shrinkage = 0.1,n.trees = 500)
  M <- gbm.perf(tmp)[1]
  mod <- gbm(DIQ010_diq~.,data = nhanes[ind!=i,-1], distribution = "adaboost", interaction.depth = 2,
             shrinkage = 0.1,n.trees = M)
  RES[ind==i, "adaboost"] <- predict(mod, newdata=nhanes[ind==i,-1], type = "response", n.trees = M)
  
  #9-logitboost
  #---------
  print("logitboost")
  tmp <- gbm(DIQ010_diq~.,data = nhanes[ind!=i,-1], distribution = "bernoulli", interaction.depth = 2,
             shrinkage = 0.1,n.trees = 500)
  M <- gbm.perf(tmp)[1]
  mod <- gbm(DIQ010_diq~.,data = nhanes[ind!=i,-1], distribution = "bernoulli", interaction.depth = 2,
             shrinkage = 0.1,n.trees = M)
  RES[ind==i, "logitboost"] <- predict(mod, newdata=nhanes[ind==i,-1], type = "response", n.trees = M)
  
  #10-svm linéaire
  #---------------
  print("svmlin")
  mod <- svm(DIQ010_diq~.,data=nhanes[ind!=i,-1], kernel="linear")
  RES[ind==i,"svmlin"] <- predict(mod,newdata = nhanes[ind==i,-1])
  
  #11-svm radial
  #-------------
  print("svmrad")
  mod <- svm(DIQ010_diq~.,data=nhanes[ind!=i,-1], kernel="radial")
  RES[ind==i,"svmrad"] <- predict(mod,newdata = nhanes[ind==i,-1])
  
  }





#METRIQUES
#---------

#courbes ROC et aires sous la courbe
#-----------------------------------
roclogcomp <- roc(RES[,1],RES[,2])
plot(roclogcomp)
roclogstep <- roc(RES[,1],RES[,3])
lines(roclogstep,col="red",lty=3,lwd=3)
roclogred15 <- roc(RES[,1],RES[,4])
lines(roclogred15,col="blue",lty=3,lwd=3)
rocrid<-roc(RES[,1],RES[,5])
lines(rocrid,col="blue")
roclas<-roc(RES[,1],RES[,6])
lines(roclas,col="orange")
rocela<-roc(RES[,1],RES[,7])
lines(rocela,col="brown")
rocfor <- roc(RES[,1],RES[,8])
lines(rocfor,col="green")
rocada<-roc(RES[,1],RES[,9])
lines(rocada,col="purple")
roclogib<-roc(RES[,1],RES[,10])
lines(roclogib,col="dodger blue")
rocsvmlin<-roc(RES[,1],RES[,11])
lines(rocsvmlin,col="grey")
rocsvmrad<-roc(RES[,1],RES[,12])
lines(rocsvmrad,col="darkgrey")

legend("bottomright", legend=c("logcomp","logstep","logred15","ridge","lasso","elastic","foret","adaboost","logitboost","svmlin","svmrad"),
       col=c("black", "red","blue","blue","orange","brown","green","purple","dodger blue","grey","darkgrey"),lty=c(1,3,3,1,1,1,1,1,1,1),lwd=c(1,2,2,1,1,1,1,1,1,2), cex=1)


auc(RES[,1],RES[,2])
auc(RES[,1],RES[,3])
auc(RES[,1],RES[,4])
auc(RES[,1],RES[,5])
auc(RES[,1],RES[,6])
auc(RES[,1],RES[,7])
auc(RES[,1],RES[,8])
auc(RES[,1],RES[,9])
auc(RES[,1],RES[,10])
auc(RES[,1],RES[,11])
auc(RES[,1],RES[,12])


monerreur <- function(X,Y,seuil=0.5){
  table(cut(X,breaks=c(0,seuil,1)),as.factor(Y))
}



.
write.csv(RES,"res_dia.csv")





#=========================================
#Traitement de l'importance des variables
#=========================================

nhanes$SEQN<-NULL
nhanes<-cbind(nhanes[,names(nhanes)!="DIQ010_diq"],DIQ010_diq=nhanes[,c("DIQ010_diq")])

str(nhanes)


XX<- model.matrix(DIQ010_diq~.,data=nhanes)
YY <- as.matrix(nhanes[,c("DIQ010_diq")])
# XX <- as.matrix(model.matrix(~.,nhanes)[,-ncol(model.matrix(~.,nhanes))])
# YY <- as.matrix(model.matrix(~.,nhanes)[,ncol(model.matrix(~.,nhanes))])

# CrÃ©ation de la fonction pour calculer l'importance des variables
# x est le modele; k est le nombre de variable
# valeur de t: (1 pour glm et randomforest, 2 pour glmnet, 3 pour gbm)
variable_imp <- function(x,k=15,t=1,mot=""){
  switch(t,
         x <- varImp(x), #on utilise varImp de Caret pour glm et randomforest
         x <- as.data.frame(as.matrix(x)), # utile pour mise au format dataframe glmnet
         x[,1] <- NULL # utile pour enlever une colonne en trop pour gbm
  )
  tempo <- cbind(row.names(x),x)
  row.names(tempo) <- NULL
  colnames(tempo)[1] <- "variable"
  colnames(tempo)[2] <-  "importance"
  tempo[,1] <- gsub("YES$","",tempo[,1],ignore.case = TRUE)
  tempo[,1] <- gsub("NO$","",tempo[,1],ignore.case = TRUE)
  tempo[,2] <- (tempo[,2]-mean(tempo[,2]))/sqrt(var(tempo[,2])) # on centre rÃ©duit
  tempo <- arrange(tempo, desc(importance))[1:k,]
  colnames(tempo)[2] <-  paste("imp",mot,sep = "_")
  tempo$rank <- seq(1:k)
  colnames(tempo)[3] <- paste("rang",mot,sep = "_")
  return(tempo)
}

#1 Importance variable pour le modele logistique
mod_log <- glm(DIQ010_diq~.,data=nhanes,family="binomial")
varimplog <- variable_imp(x=mod_log,t=1,mot="log")
varimplog

#2Importance variable pour le modele logistique step
mod_log <- glm(DIQ010_diq ~ RIDAGEYR_demo + MCQ080_mcq + BPQ020_bpq + 
                 BPQ080_bpq + DR1TSUGR_dr1tot + BMXWT_bmx + INDFMPIR_demo + 
                 DR1TCHOL_dr1tot + DR1TALCO_dr1tot + RIAGENDR_demo + Var_ACTIVITE + 
                 Var_DENTISTE + Var_TENSIONDI + SLQ050_slq + DR1TACAR_dr1tot + 
                 BMXBMI_bmx + Var_SITUATION + HOD050_hoq + DR1TMOIS_dr1tot + 
                 DR1.320Z_dr1tot + DR1TCAFF_dr1tot + DR1TPROT_dr1tot + DR1TIRON_dr1tot + 
                 MCQ160N_mcq + DR1TRET_dr1tot + MCQ160D_mcq + MCQ160C_mcq + 
                 DR1TSELE_dr1tot + DR1TFOLA_dr1tot,data=nhanes,family="binomial")

varimplogstep <- variable_imp(x=mod_log,t=1,mot="logstep")
varimplogstep

# variable     imp_log rang_log
# 1    RIDAGEYR_demo  2.99402171        1
# 2  DR1TSUGR_dr1tot  2.04803630        2
# 3       BPQ080_bpq  1.76073033        3
# 4       MCQ080_mcq  1.50381126        4
# 5       BPQ020_bpq  1.25641657        5
# 6  DR1TALCO_dr1tot  0.81241564        6
# 7  DR1TMOIS_dr1tot  0.45875503        7
# 8    RIAGENDR_demo  0.12360083        8
# 9  DR1.320Z_dr1tot  0.11885295        9
# 10 DR1TCHOL_dr1tot  0.09586159       10
# 11 DR1TPROT_dr1tot  0.06684281       11
# 12   INDFMPIR_demo  0.04296765       12
# 13 DR1TIRON_dr1tot -0.11450863       13
# 14   Var_TENSIONDI -0.16515044       14
# 15 DR1TCAFF_dr1tot -0.17212308       15


#3)Importance variable pour le modele ridge
tmp <- cv.glmnet(XX,YY,alpha=0,family="binomial")
mod_ridge  <- glmnet(XX,YY,alpha=0,lambda=tmp$lambda.min, family="binomial")
varimpridge <- variable_imp(x=abs(mod_ridge$beta),t=2,mot="ridge")

#4)Importance variable pour le modele lasso (colnames(XX)[mod_lasso$beta@i])
tmp <- cv.glmnet(XX,YY, alpha=1, family="binomial")
mod_lasso <- glmnet(XX,YY,alpha=1, lambda =tmp$lambda.1se,family="binomial" )
varimplasso <- variable_imp(x=abs(mod_lasso$beta),t=2, mot="lasso")

#5)Importance variable pour le modele elastic
tmp <- cv.glmnet(XX,YY, alpha=0.5, family="binomial")
mod_elastic <- glmnet(XX,YY,alpha = 0.5, lambda = tmp$lambda.min, family="binomial")
varimpelastic <- variable_imp(abs(mod_elastic$beta),t=2,mot = "elastic")

#6)Importance variable pour le modele Foret
mod_foret <- randomForest(factor(DIQ010_diq)~., data = nhanes)
varimpforet <- variable_imp(mod_foret,t=1,mot="foret")

#7)Importance variable pour le modele adaboost
tmp <- gbm(as.numeric(DIQ010_diq)~.,data = nhanes, distribution = "adaboost", interaction.depth = 2,
           shrinkage = 0.1,n.trees = 500)
M <- gbm.perf(tmp)[1]
mod_adaboost <- gbm(as.numeric(DIQ010_diq)~.,data = nhanes, distribution = "adaboost", interaction.depth = 2,
                    shrinkage = 0.1,n.trees = M)
varimpada <- variable_imp(summary(mod_adaboost),t=3, mot="adaboost")

#8)Importance variable pour le modele logitboost
tmp <- gbm(as.numeric(DIQ010_diq)~.,data=nhanes, distribution="bernoulli", interaction.depth = 2,
           shrinkage=0.1,n.trees=500)
M <- gbm.perf(tmp)[1]
mod_logitboost <- gbm(as.numeric(DIQ010_diq)~.,data=nhanes, distribution="bernoulli", interaction.depth = 2,
                     shrinkage=0.1,n.trees=M)
varimplogibo <- variable_imp(summary(mod_logitboost),t=3,mot="logitboost")

#9)Importance variable pour le modele svm linéaire
# mod_svmlin <- svm(DIQ010_diq~.,data=nhanes2, kernel="linear",probability=T)
# varimpsvmlin <- variable_imp(summary(mod_svmlin),t=3,mot="svmlin")

#10)Importance variable pour le modele svm radial
# mod_svmrad <- svm(DIQ010_diq~.,data=nhanes2, kernel="radial",probability=T)
# varimpsvmrad <- variable_imp(summary(mod_svmrad),t=3,mot="svmrad")

# Croisement des tables d'importance des variables
choix_var <- varimplog %>%
  full_join(varimplogstep) %>%
  full_join(varimpridge) %>%
  full_join(varimplasso) %>%
  full_join(varimpelastic) %>%
  full_join(varimpforet) %>%
  full_join(varimpada) %>%
  full_join(varimplogibo) 
  # full_join(varimpsvmlin) %>%
  # full_join(varimpsvmrad)

choix_var <- cbind(choix_var[,1],choix_var[,c(which(grepl("^imp",names(choix_var))))])

choix_var <- as.data.frame(choix_var)
names(choix_var)[1] <- "variable"
write.csv2(choix_var,"choix_var_dia.csv",row.names = FALSE)




#tuning de la random forest
#--------------------------

# set.seed(1234)
# mod_complet <- randomForest(DIQ010_diq~.,data=nhanesfac[,-1],ntree=1000)
# plot(mod_complet$err.rate[, 1], type = "l", xlab = "nombre d'arbres", ylab = "erreur OOB")
# legend("topright", legend=c("mod_complet"),col=c("red"),lty=1,cex=1)
# 
# print(mod_complet)
# 
# #? l'issue de cette ?tape je retiens le mod?le r?duit et ntree=200
# 
# 
# set.seed(1234)
# oob<-list()
# for (i in seq(2,10,1)) {
#   mod_complet <- randomForest(DIQ010_diq~.,data=nhanesfac[,-1],ntree=200,mtry=i)  
#   oob[[i]]<-mod_complet$err.rate[, 1][200]
#   }
# 
# oob 
# oob<-unlist(oob)
# plot(x=seq(2,10,1),y=oob,xlab="mtry")
# lines(x=seq(2,10,1),y=oob)
# #mtry=6 donne le meilleur résultat








#Utilisation de Metrics et ModelMetrics 
#--------------------------------------
# #Ã  partir de RES qui contient les probabilitÃ©s, je crÃ©e RESp qui contient les prÃ©dictions sous un seuil donnÃ©
# RESp<-data.table(RES)
# seuil<-0.5
# RESp[logcomp>=seuil,logcompp:=1,]
# RESp[logcomp<seuil,logcompp:=0,]
# RESp[logstep>=seuil,logstepp:=1,]
# RESp[logstep<seuil,logstepp:=0,]
# RESp[logred15>=seuil,logstepp:=1,]
# RESp[logred15<seuil,logstepp:=0,]
# RESp[foret>=seuil,foretp:=1,]
# RESp[foret<seuil,foretp:=0,]
# RESp[ridge>=seuil,ridgep:=1,]
# RESp[ridge<seuil,ridgep:=0,]
# RESp[lasso>=seuil,lassop:=1,]
# RESp[lasso<seuil,lassop:=0,]
# RESp[elastic>=seuil,elasticp:=1,]
# RESp[elastic<seuil,elasticp:=0]
# RESp[adaboost>=seuil,adaboostp:=1,]
# RESp[adaboost<seuil,adaboostp:=0,]
# RESp[logitboost>=seuil,logitboostp:=1,]
# RESp[logitboost<seuil,logitboostp:=0]
# 
# RESp<-data.frame(RESp)
# 
# str(RESp)


# resultat<-data.frame(logcomp=0,logstep=0,foret=0,ridge=0,lasso=0,elastic=0,adaboost=0,logitboost=0)
# for (i in 1:8){
#   resultat[1,i]<-Metrics::auc(RESp[,1],RESp[,i+1])  
#   resultat[2,i]<-Metrics::accuracy(RESp[,1],RESp[,i+9])
#   resultat[3,i]<-Metrics::ce(RESp[,1],RESp[,i+9])
#   resultat[4,i]<-ModelMetrics::sensitivity(RESp[,1],RESp[,i+1],cutoff=0.25)
#   resultat[5,i]<-ModelMetrics::specificity(RESp[,1],RESp[,i+1],cutoff=0.25)
#   resultat[6,i]<-ModelMetrics::logLoss(RESp[,1],RESp[,i+9])
#   }
# rownames(resultat)<-c("auc","accuracy","taux mal classÃ©s","sensibilitÃ©","spÃ©cificitÃ©","logloss")
# resultat
# 
# 
# 
# 
# 
# #1/MATRICES DE CONFUSION
# #prÃ©alablement on calcule les taux de mal class?s triviaux
# nhanes<-data.table(nhanes)
# nhanes[,.N,]
# nhanes[DIQ010_diq=="0",.N,]
# nhanes[DIQ010_diq=="1",.N,]
# taux_positifs<-nhanes[DIQ010_diq=="1",.N,]/nhanes[,.N,]*100
# taux_positifs
# taux_negatifs<-nhanes[DIQ010_diq=="0",.N,]/nhanes[,.N,]*100
# taux_negatifs
# #si tout le monde est class? en n?gatif alors on commet une erreur de 4376/5221=16,2%
# #si tout le monde est class? en positif alors on commet une erreur de 4376/5221=83,8%
# #la logistique am?liore un petit peu le mod?le trivial avec un taux de mal class?s de 14,55 pour le seuil de 0.5
# 
# 
# monerreur <- function(X,Y,seuil){
#   t<-table(cut(X,breaks=c(0,seuil,1)),as.factor(Y))
#   df<-data.frame(table(cut(X,breaks=c(0,seuil,1)),as.factor(Y)))
#   taux_mal_classes<-(df[2,3]+df[3,3])/sum(df$Freq)*100
#   liste<-list(t,taux_mal_classes)
#   return(liste)
#   }
# 
# monerreur(RES[,3],RES[,1],0.5)
# 
# taux_synthese<-data.frame()
# for (i in 2:9) {
#   taux_synthese[1,i-1]<-taux_negatifs
#   taux_synthese[11,i-1]<-taux_positifs
#   k<-2
#   for (j in seq(0.1,0.9,0.1)) {
#     taux_synthese[k,i-1]<-monerreur(RES[,i],RES[,1],j)[[2]]
#     k<-k+1
#   }  }
# colnames(taux_synthese)<-c("logcomp","logstep","foret","ridge","lasso","elastic","adaboost","logitboost")
# rownames(taux_synthese)<-seq(0,1,0.1)
# taux_synthese
# 
# 
# #Graphique des taux de mal classÃ©s
# #---------------------------------
# plot(seq(0,1,0.1),taux_synthese$logcomp,xlab="seuil",ylab="taux de mal classÃ©s (%)")
# lines(seq(0,1,0.1),taux_synthese$logcomp,lwd=2)
# lines(seq(0,1,0.1),taux_synthese$logstep,lty=3,lwd=5,col="red")
# lines(seq(0,1,0.1),taux_synthese$foret,col="green",lwd=2)
# lines(seq(0,1,0.1),taux_synthese$ridge,col="blue",lwd=2)
# lines(seq(0,1,0.1),taux_synthese$lasso,col="orange",lwd=2)
# lines(seq(0,1,0.1),taux_synthese$elastic,col="brown",lwd=2)
# lines(seq(0,1,0.1),taux_synthese$adaboost,col="purple",lwd=2)
# lines(seq(0,1,0.1),taux_synthese$logitboost,col="dodger blue",lwd=2)
# 
# legend("topright", legend=c("logcomp","logstep","foret","ridge","lasso","elastic","adaboost","logitboost"),
#        col=c("black", "red","green","blue","orange","brown","purple","dodger blue"),lwd=2,lty=c(1,3,1,1,1,1,1,1), cex=1)






