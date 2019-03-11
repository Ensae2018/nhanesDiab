#**********************************************************************************************
#Ce programme utilise le dataset nhanes_16012019 avec 8339 individus et 150 variables en entree
#Etude DIABETE pour la presentation du 27 fevrirer 2019
#**********************************************************************************************

#permet de ne pas limiter le nombre des ecritures dans la console lors des executions
options(max.print=1000000)

#declaration des librairies
library(data.table)
library(mice)
library(randomForest)
library(glmnet)
library(pROC)
library(gbm)
library(foreach)
library(ModelMetrics)
library(Metrics)


#lecture des donnees
nhanes<-fread("nhanes_16012019.csv")
nhanes[,V1:=NULL,]
dim(nhanes)#8339x150



#0- Statistiques sur le dataset nhanes 8339x150
#------------------------------------------
#nombre total de data:1 250 850
nrow(nhanes)*ncol(nhanes)
#nombre total de NA:351 779 soit
nhanes[,sum(is.na(.SD)),]
#poids des NA:28%
nhanes[,sum(is.na(.SD)),]/(nrow(nhanes)*ncol(nhanes))*100


#1- Recodage / sélection additionnelle par les NA et les variables au cas par cas
#--------------------------------------------------------------------------------

#je conserve les individus qui ont un DIQ010_diq renseigne ? 1 (Yes)  2 (No) 3 (Borderline)
nhanes1<-nhanes[!(is.na(DIQ010_diq)) & !(DIQ010_diq %in% c("7","9")),,]
dim(nhanes1)#8115x150

#les jeunes de moins de 18 ans generent la moitie du nombre total des NA alors que seulement 25 d'entre eux sont diabetiques
#je decide de les eliminer 
nhanes1[,sum(is.na(.SD)),]
nhanes1[RIDAGEYR_demo<18,sum(is.na(.SD)),]
nhanes1[RIDAGEYR_demo>=18,sum(is.na(.SD)),]
nhanes1[RIDAGEYR_demo<18 & DIQ010_diq %in% c("1","3"),.N,]
nhanes1<-nhanes1[RIDAGEYR_demo>=18,,]
dim(nhanes1)#5264x150

#j'enleve les 43 individus avec un regime low sugar car l'alimentation est modifiee dans ce cas (DRQSDT4_dr1tot=4)
nhanes1[DRQSDT4_dr1tot=="4",.N,]
nhanes1<-nhanes1[is.na(DRQSDT4_dr1tot),,]
dim(nhanes1)#5221x150

#traitement de la table alq (alcool)
apply(nhanes1[,.(ALQ120Q_alq,ALQ120U_alq,ALQ130_alq),],2,function (x) sum(is.na(x)))

nhanes1[ALQ120U_alq=="1",temp:=52,]
nhanes1[ALQ120U_alq=="2",temp:=12,]
nhanes1[ALQ120U_alq=="3",temp:=1,]

nhanes1[!(is.na(ALQ120Q_alq) | ALQ120Q_alq=="777" | ALQ120Q_alq=="999" | is.na(ALQ130_alq) | ALQ130_alq=="777" | ALQ130_alq=="999"),Var_ALQ:=ALQ120Q_alq*temp*ALQ130_alq,]
nhanes1[,c("temp","ALQ120Q_alq","ALQ120U_alq","ALQ130_alq"):=NULL,]
dim(nhanes1)#5221x148


#traitement de la table demo (demographique)
apply(nhanes1[,.(RIDSTATR_demo,RIAGENDR_demo,RIDAGEYR_demo,DMDEDUC3_demo,DMDEDUC2_demo,
                 DMDMARTL_demo,DMDHHSIZ_demo, 
                 DMDFMSIZ_demo,INDHHIN2_demo,INDFMIN2_demo,INDFMPIR_demo),],2,function (x) sum(is.na(x)))

nhanes1[DMDEDUC2_demo %in% c("1","2","3"),Var_EDUCATION:="secondaire",]
nhanes1[DMDEDUC2_demo %in% c("4","5"),Var_EDUCATION:="superieur",]
nhanes1[is.na(DMDEDUC2_demo) | DMDEDUC2_demo %in% c("7","9"),Var_EDUCATION:=NA,]
nhanes1[,c("DMDEDUC2_demo","DMDEDUC3_demo"):=NULL,]
nhanes1$Var_EDUCATION<-factor(nhanes1$Var_EDUCATION)
dim(nhanes1)#5221x147

#traitement de la table bmx (body measures)
#on a deja le BMI quantitatif pour toute la population avec BMXBMI; BMDBMIC ne porte que sur les jeunes
nhanes1[,BMDBMIC_bmx:=NULL,]
dim(nhanes1)#5221x146

#traitement de la table ocq (travail ou recherche d'emploi);je regroupe les 7/9/NA vu leur faible nombre
nhanes1[OCD150_ocq %in% c("1","2"),Var_TRAVAIL:="oui",]
nhanes1[OCD150_ocq %in% c("3", "4"),Var_TRAVAIL:="non",]
nhanes1[OCD150_ocq %in% c("7", "9") | is.na(OCD150_ocq),Var_TRAVAIL:=NA,]
nhanes1[,c("OCD150_ocq","OCQ180_ocq"):=NULL,]
nhanes1$Var_TRAVAIL<-factor(nhanes1$Var_TRAVAIL)
dim(nhanes1)#5221x145

#traitement de la table duq (drogue)
#seules les personnes de 18 a 69 ans sont incluses dans le fichier
nhanes1[DUQ200_duq=="1" | DUQ240_duq=="1" | DUQ370_duq=="1",Var_DROGUE:="oui",]
nhanes1[DUQ200_duq=="2" & DUQ240_duq=="2" & DUQ370_duq=="2",Var_DROGUE:="non",]
nhanes1[!(Var_DROGUE=="oui" | Var_DROGUE=="non"),Var_DROGUE:=NA,]
nhanes1[,c("DUQ200_duq","DUQ240_duq","DUQ370_duq"):=NULL,]
nhanes1$Var_DROGUE<-factor(nhanes1$Var_DROGUE)
table(nhanes1$Var_DROGUE,useNA="always")
dim(nhanes1)#5221x143

#traitement de la table dpq (depressif)
nhanes1[DPQ020_dpq=="0",Var_DEPRESSION:="non",]
nhanes1[DPQ020_dpq %in% c("1","2","3"),Var_DEPRESSION:="oui",]
nhanes1[DPQ020_dpq %in% c("7","9") | is.na(DPQ020_dpq),Var_DEPRESSION:=NA, ]
nhanes1[,c("DPQ020_dpq"):=NULL,]
nhanes1$Var_DEPRESSION<-factor(nhanes1$Var_DEPRESSION)
table(nhanes1$Var_DEPRESSION,useNA="always")
dim(nhanes1)#5221x143

#traitement du statut marital de la table demo (DMDMARTL_demo)
nhanes1[DMDMARTL_demo %in% c("1","6"),Var_SITUATION:="couple",]
nhanes1[DMDMARTL_demo %in% c("2","3","4","5"),Var_SITUATION:="seul",]
nhanes1[DMDMARTL_demo %in% c("77","99") | is.na(DMDMARTL_demo),Var_SITUATION:=NA,]
nhanes1[,c("DMDMARTL_demo"):=NULL,]
nhanes1$Var_SITUATION<-factor(nhanes1$Var_SITUATION)
dim(nhanes1)#5221x143

#traitement du sexe de la table demo
nhanes1$RIAGENDR_demo<-as.factor(nhanes1$RIAGENDR_demo)

#traitement de la table diq (DIQ160_diq); mieux vaut rester sur le diabete pur
nhanes1[,c("DIQ160_diq"):=NULL,]
dim(nhanes1)#5221x142

#la table slq ne posera pas de probleme vu le faible nombre de NA
nhanes1[is.na(SLQ050_slq),.N,]
nhanes1[is.na(SLD012_slq),.N,]

#traitement de la table smq
nhanes1[SMQ040_smq %in% c("1","2"),Var_FUMEUR:="oui",]
nhanes1[SMQ040_smq %in% c("3"),Var_FUMEUR:="non",]
nhanes1[,c("SMQ040_smq"):=NULL,]
nhanes1$Var_FUMEUR<-factor(nhanes1$Var_FUMEUR)
table(nhanes1$Var_FUMEUR,useNA="always")
dim(nhanes1)#5221x142


#traitement de la table smqfam
nhanes1[SMD460_smqfam=="0",Var_COFUMEUR:="non",]
nhanes1[SMD460_smqfam %in% c("1","2","3"),Var_COFUMEUR:="oui",]
nhanes1[,c("SMD460_smqfam"):=NULL,]
nhanes1$Var_COFUMEUR<-factor(nhanes1$Var_COFUMEUR)
table(nhanes1$Var_COFUMEUR,useNA="always")
dim(nhanes1)#5221x142


#exploration des variables du nombre des personnes qui vivent dans la famille ou le household
nhanes1[DMDHHSIZ_demo!=DMDFMSIZ_demo,.(DMDHHSIZ_demo,DMDFMSIZ_demo),]
nhanes1[DMDFMSIZ_demo>DMDHHSIZ_demo,.N,]
dim(nhanes1)#5221x142

#traitement de la consommation habituelle de sel
nhanes1[DBD100_dr1tot %in% c("7","9") | is.na(DBD100_dr1tot),DBD100_dr1tot:=NA,]

#traitement des tensions arterielles
#nhanes1[,Var_TENSIONSY:=mean(BPXSY1_bpx,BPXSY2_bpx,BPXSY3_bpx),] ne marche pas
nhanes1[,Var_TENSIONSY:=(BPXSY1_bpx+BPXSY2_bpx+BPXSY3_bpx)/3,]
nhanes1[,Var_TENSIONDI:=(BPXDI1_bpx+BPXDI2_bpx+BPXDI3_bpx)/3,]
nhanes1[,c("BPXSY1_bpx","BPXSY2_bpx","BPXSY3_bpx","BPXSY4_bpx","BPXDI1_bpx","BPXDI2_bpx","BPXDI3_bpx","BPXDI4_bpx"):=NULL,]
dim(nhanes1)#5221x136

#traitement de l'argent destine a la consommation
nhanes1[,Var_ARGENTALIM:=(CBD071_cbq+CBD111_cbq+CBD121_cbq+CBD131_cbq),]
nhanes1[,c("CBD071_cbq","CBD111_cbq","CBD121_cbq","CBD131_cbq","CBD091_cbq"):=NULL,]
nhanes1[is.na(Var_ARGENTALIM),.N,]
dim(nhanes1)#5221x132

#traitement de la table OHXREF
nhanes1[OHAREC_ohxref=="1",Var_DENTISTE:="Immediatement",]
nhanes1[OHAREC_ohxref=="2",Var_DENTISTE:="Visite sous 2 semaines",]
nhanes1[OHAREC_ohxref=="3",Var_DENTISTE:="Visite sous convenance",]
nhanes1[OHAREC_ohxref=="4",Var_DENTISTE:="Pas de visite prochaine",]
nhanes1[,c("OHAREC_ohxref"):=NULL,]
nhanes1$Var_DENTISTE<-as.factor(nhanes1$Var_DENTISTE)
dim(nhanes1)#5221x132

#traitement de l'assurance
nhanes1[HIQ011_hiq=="1",Var_ASSURE:="oui",]
nhanes1[HIQ011_hiq=="2",Var_ASSURE:="non",]
nhanes1[HIQ011_hiq %in% c("7","9"),Var_ASSURE:=NA,]
nhanes1[,c("HIQ011_hiq"):=NULL,]
nhanes1$Var_ASSURE<-as.factor(nhanes1$Var_ASSURE)
dim(nhanes1)#5221x132

#traitement de statut de propriete ou location
nhanes1[HOQ065_hoq=="1",Var_MAISON:="proprietaire",]
nhanes1[HOQ065_hoq=="2",Var_MAISON:="en location",]
nhanes1[HOQ065_hoq %in% c("7","9"),Var_MAISON:=NA,]
nhanes1[,c("HOQ065_hoq"):=NULL,]
nhanes1$Var_MAISON<-as.factor(nhanes1$Var_MAISON)
dim(nhanes1)#5221x132

#traitement de la table hoq
nhanes1[HOD050_hoq %in% c("777","999"),HOD050_hoq:=NA,]

#traitement de la table paq
nhanes1[PAQ605_paq=="1" | PAQ620_paq=="1" | PAQ635_paq=="1" | PAQ650_paq=="1" | PAQ665_paq=="1",Var_ACTIVITE:="actif",]
nhanes1[!(PAQ605_paq=="1" | PAQ620_paq=="1" | PAQ635_paq=="1" | PAQ650_paq=="1" | PAQ665_paq=="1"),Var_ACTIVITE:="inactif",]
nhanes1[,c("PAQ605_paq","PAQ620_paq","PAQ635_paq","PAQ650_paq","PAQ665_paq"):=NULL,]
nhanes1$Var_ACTIVITE<-as.factor(nhanes1$Var_ACTIVITE)
nhanes1[PAD680_paq=="7777" | PAD680_paq=="9999" | is.na(PAD680_paq),PAD680_paq:=NA,]
summary(nhanes1$PAD680_paq)


#mise en facteurs des variables ad hoc
nhanes1[BPQ020_bpq %in% c("7","9"),BPQ020_bpq:=NA,]
nhanes1[BPQ080_bpq %in% c("7","9"),BPQ080_bpq:=NA,]
nhanes1[MCQ080_mcq %in% c("7","9"),MCQ080_mcq:=NA,]
nhanes1[MCQ160A_mcq %in% c("7","9"),MCQ160A_mcq:=NA,]
nhanes1[MCQ160B_mcq %in% c("7","9"),MCQ160B_mcq:=NA,]
nhanes1[MCQ160C_mcq %in% c("7","9"),MCQ160C_mcq:=NA,]
nhanes1[MCQ160D_mcq %in% c("7","9"),MCQ160D_mcq:=NA,]
nhanes1[MCQ160E_mcq %in% c("7","9"),MCQ160E_mcq:=NA,]
nhanes1[MCQ160F_mcq %in% c("7","9"),MCQ160F_mcq:=NA,]
nhanes1[MCQ160G_mcq %in% c("7","9"),MCQ160G_mcq:=NA,]
nhanes1[MCQ160M_mcq %in% c("7","9"),MCQ160M_mcq:=NA,]
nhanes1[MCQ160N_mcq %in% c("7","9"),MCQ160N_mcq:=NA,]
nhanes1[MCQ160A_mcq %in% c("7","9"),MCQ160A_mcq:=NA,]
nhanes1[SLQ050_slq %in% c("7","9"),SLQ050_slq:=NA,]
nhanes1[HEQ010_heq %in% c("7","9"),HEQ010_heq:=NA,]
nhanes1[HEQ030_heq %in% c("7","9"),HEQ030_heq:=NA,]

nhanes1$BPQ020_bpq<-as.factor(nhanes1$BPQ020_bpq)
nhanes1$BPQ080_bpq<-as.factor(nhanes1$BPQ080_bpq)
nhanes1$MCQ080_mcq<-as.factor(nhanes1$MCQ080_mcq)
nhanes1$MCQ160A_mcq<-as.factor(nhanes1$MCQ160A_mcq)
nhanes1$MCQ160B_mcq<-as.factor(nhanes1$MCQ160B_mcq)
nhanes1$MCQ160C_mcq<-as.factor(nhanes1$MCQ160C_mcq)
nhanes1$MCQ160D_mcq<-as.factor(nhanes1$MCQ160D_mcq)
nhanes1$MCQ160E_mcq<-as.factor(nhanes1$MCQ160E_mcq)
nhanes1$MCQ160F_mcq<-as.factor(nhanes1$MCQ160F_mcq)
nhanes1$MCQ160G_mcq<-as.factor(nhanes1$MCQ160G_mcq)
nhanes1$MCQ160M_mcq<-as.factor(nhanes1$MCQ160M_mcq)
nhanes1$MCQ160N_mcq<-as.factor(nhanes1$MCQ160N_mcq)
nhanes1$SLQ050_slq<-as.factor(nhanes1$SLQ050_slq)
nhanes1$HEQ010_heq<-as.factor(nhanes1$HEQ010_heq)
nhanes1$HEQ030_heq<-as.factor(nhanes1$HEQ030_heq)

#traitement de suppression de variables
#--------------------------------------

#les 16 variables DR1TOT qui sont mal renseignees ou ne nous apportent rien pour l'etude diabete
#les 2  variables de prise de tension n?4 BPXSY4 et BPXDI4 tres mal renseignees (8042NA) (on se basera sur 3 prises de tensions)
#les 3 variable ecq ne concernent que les enfants de 0 a 15 ans donc on ne pourra pas traiter un modele general (adultes+enfants)
#DMDHHSIZ_demo est toujours supeieur ou egal a DMDFMSIZ_demo donc je garde  DMDHHSIZ_demo
#pour les revenus je garde les revenus du household INDHHIN2_demo par coherence; data aussi bien renseignee que la famille

dput(names(nhanes1))
length(dput(names(nhanes1)))
nhanes1<-nhanes1[,c("SEQN",
                    # "RIDSTATR_demo",
                    "RIAGENDR_demo", "RIDAGEYR_demo", "DMDHHSIZ_demo",
                    # "DMDFMSIZ_demo",
                    "INDHHIN2_demo",
                    # "INDFMIN2_demo",
                    "INDFMPIR_demo","BMXWT_bmx", "BMXHT_bmx", "BMXBMI_bmx",
                    "BPQ020_bpq", "BPQ080_bpq",
                    "DIQ010_diq",
                    # "ECD010_ecq", "ECQ020_ecq", "ECD070B_ecq",
                    "HEQ010_heq", "HEQ030_heq",
                    # "HIQ011_hiq",
                    "HOD050_hoq",
                    # "HOQ065_hoq",
                    # "IMQ011_imq", "IMQ020_imq", "MCQ010_mcq", 
  "MCQ080_mcq", "MCQ160A_mcq", "MCQ160N_mcq", "MCQ160B_mcq", "MCQ160C_mcq", 
  "MCQ160D_mcq", "MCQ160E_mcq", "MCQ160F_mcq", "MCQ160G_mcq", "MCQ160M_mcq", 
  # "MCQ160K_mcq", "MCQ160L_mcq", "MCQ220_mcq",
  # "MCQ230A_mcq", "MCQ230B_mcq", "MCQ230C_mcq", "MCQ230D_mcq",
  # "OHAREC_ohxref",
  #"PAQ605_paq", "PAQ620_paq", "PAQ635_paq", "PAQ650_paq", "PAQ665_paq","PAD680_paq", 
  "SLD012_slq", "SLQ050_slq",
  # "SMD030_smq", "SMQ720_smqrtu",
  "LBXBPB_pbcd", "LBDBPBSI_pbcd", "LBDBPBLC_pbcd",
  "DR1TKCAL_dr1tot", 
  "DR1TPROT_dr1tot", "DR1TCARB_dr1tot", "DR1TSUGR_dr1tot", "DR1TFIBE_dr1tot", 
  "DR1TTFAT_dr1tot", "DR1TSFAT_dr1tot", "DR1TMFAT_dr1tot", "DR1TPFAT_dr1tot", 
  "DR1TCHOL_dr1tot", "DR1TATOC_dr1tot", "DR1TATOA_dr1tot", "DR1TRET_dr1tot", 
  "DR1TVARA_dr1tot", "DR1TACAR_dr1tot", "DR1TBCAR_dr1tot", "DR1TCRYP_dr1tot", 
  "DR1TLYCO_dr1tot", "DR1TLZ_dr1tot", "DR1TVB1_dr1tot", "DR1TVB2_dr1tot", 
  "DR1TNIAC_dr1tot", "DR1TVB6_dr1tot", "DR1TFOLA_dr1tot", "DR1TFA_dr1tot", 
  "DR1TFF_dr1tot", "DR1TFDFE_dr1tot", "DR1TCHL_dr1tot", "DR1TVB12_dr1tot", 
  "DR1TB12A_dr1tot", "DR1TVC_dr1tot", "DR1TVD_dr1tot", "DR1TVK_dr1tot", 
  "DR1TCALC_dr1tot", "DR1TPHOS_dr1tot", "DR1TMAGN_dr1tot", "DR1TIRON_dr1tot", 
  "DR1TZINC_dr1tot", "DR1TCOPP_dr1tot", "DR1TSODI_dr1tot", "DR1TPOTA_dr1tot", 
  "DR1TSELE_dr1tot", "DR1TCAFF_dr1tot", "DR1TTHEO_dr1tot", "DR1TALCO_dr1tot", 
  "DR1TMOIS_dr1tot", "DR1.320Z_dr1tot",
  # "DRABF_dr1tot","DRDINT_dr1tot", 
  "DBD100_dr1tot",
  # "DRQSDIET_dr1tot", "DRQSDT1_dr1tot", "DRQSDT2_dr1tot", 
  # "DRQSDT3_dr1tot", "DRQSDT4_dr1tot", "DRQSDT5_dr1tot", "DRQSDT6_dr1tot", 
  # "DRQSDT7_dr1tot", "DRQSDT8_dr1tot", "DRQSDT9_dr1tot", "DRQSDT10_dr1tot", 
  # "DRQSDT11_dr1tot", "DRQSDT12_dr1tot", "DRQSDT91_dr1tot",
  "Var_ALQ", 
  "Var_EDUCATION", "Var_TRAVAIL", "Var_DROGUE", "Var_DEPRESSION", 
  "Var_SITUATION", "Var_FUMEUR", "Var_COFUMEUR", "Var_TENSIONSY", 
  "Var_TENSIONDI", "Var_ARGENTALIM","Var_DENTISTE","Var_ASSURE","Var_ACTIVITE")]


dim(nhanes1)#5221x92

dt<-data.table(var=names(nhanes1),nbna=apply(nhanes1,2,function (x) {sum(is.na(x))}))
dt

nhanes1$RIAGENDR_demo<-factor(nhanes1$RIAGENDR_demo)

#il faut passer le Y en numerique 0/1 pour le glm 
nhanes1[nhanes1$DIQ010_diq =="2",12]<-0
nhanes1[nhanes1$DIQ010_diq %in% c("1","3"),12]<-1

str(nhanes1)

nhanes1<-data.frame(nhanes1)
dim(nhanes1)




#je ne garde que les variables qui ont moins de 10% de NA et performe le mice dessus
#-----------------------------------------------------------------------------------

dt<-data.table(var=names(nhanes1),nbna=apply(nhanes1,2,function (x) {sum(is.na(x))}))
dt
var_ok<-dt[,nbna<=0.1*nrow(nhanes1)]
var_ok
nhanes1<-nhanes1[,var_ok]
dim(nhanes1)#nhanes3 comporte maintenant 5221 individus et 85 variables

write.csv(nhanes1,"nhanes_diab_mice_avant.csv")

dim(nhanes1)#5221x85


class(nhanes1$SLD012_slq)
summary(nhanes1$SLD012_slq)



#2- Imputation par mice
#----------------------


imp<-mice(nhanes1,m=1)

nhanes2<-complete(imp)
dim(nhanes2)#5221x85


str(nhanes1)
dt<-data.table(var=names(nhanes2),nbna=apply(nhanes2,2,function (x) {sum(is.na(x))}))
dt[,sum(is.na(.SD)),]

write.csv(nhanes2,"nhanes_diab_mice_apres.csv")

#vérification de la non déformation de la distribution des variables
str(nhanes2)
barplot(table(nhanes2$Var_ASSURE,useNA="always"))
barplot(table(nhanes1$Var_ASSURE,useNA="always"))
table(nhanes1$Var_TRAVAIL,useNA="always")


nhanes3<-nhanes2
dim(nhanes3)


#=====================================================
#=====================================================

?rownames

#3- Modèles
#----------

nhanes3 <- read.csv("nhanes_diab_mice_apres.csv")
str(nhanes3)
nhanes3[,1]<-NULL
dim(nhanes3)#5221x84
str(nhanes3)
class(summary(nhanes3$SEQN))
summary(nhanes3$DR1TKCAL_dr1tot)
xdia<-nhanes3[,"Var_ACTIVITE"]
d<-data.frame(table(xdia))
d
summary(xdia)
#dernieres vérifications avant lancement des modèles
dim(nhanes3)#5221x84
str(nhanes3)
colnames(summary(xdia))

?ifelse


#Regression logistique avec les 84 variables
#-------------------------------------------
reglog<-glm(DIQ010_diq~.,data=nhanes3[,-1],family="binomial")
summary(reglog)

#les variables qui ressortent de la regression logistique complete sont:
#3*
#age (ridageyr_demo)
#hypertension NON (bpq020_bpq)
#cholesterol NON (bpq080_bpq)
#surpoids NON (mcq080_mcq)
#hydratation (dr1tmois_dr1tot)

#2*
#sexe femme (riagendr_demo)
#indice de pauvret? (indfmpir_demo)
#sucre (dr1tsugr_dr1tot)
#eau (dr1.320z_dr1tot)
#tension diastolique (var_tensiondi)

#1*
#nb de pieces dans la maison (hod050_hoq)
#goutte NON (mcq160n_mcq)
#angine de poitrine NON (mcq160d_mcq)
#troubles du sommeil NON (slq050_slq)
#fer (dr1tiron)
#situation seul (var_situation)
#activite NON (var_activite)


#Regression logistique améliorée avec le step
#--------------------------------------------

null=glm(DIQ010_diq~1, data=nhanes3[,-1],family=binomial)
null
full=glm(DIQ010_diq~., data=nhanes3[,-1],family=binomial)
full
summary(full)

bestmodelfor=step(null, scope=list(lower=null, upper=full), direction="forward")
summary(bestmodelfor)


#3* (12)
#age                              (ridageyr_demo)
#surpoids NON                     (mcq080_mcq)
#hypertension NON                 (bpq020_bpq)
#cholesterol NON                  (bpq080_bpq)
#sucre                            (dr1tsugr_dr1tot)
#taux de pauvret?                 (indfmpir_demo)
#cholesterol                      (dr1tchol_dr1tot)
#alcool                           (dr1talco_dr1tot)
#sexe f?minin                     (riagendr_demo)
#humidite                         (dr1tmois_dr1tot)
#eau                              (dr1.320z_dr1tot)
#proteines                        (dr1tprot_dr1tot)

#2* (5)
#activit? NON                     (var_activite)
#troubles du sommeil NON          (slq050_slq)
#situation seul                   (var_situation)
#cafeine                          (dr1tcaff_dr1tot)
#fer                              (dr1tiron_dr1tot)

#1* (6)
#alpha-carotene                   (dr1tacar_dr1tot)
#IMC                              (bmxbmi_bmx)
#angine de poitrine               (mcq160d_mcq)
#goutte                           (mcq160n_mcq)
#tension diastolique              (var_tension)
#nombre de pieces dans la maison  (hod050_hoq)

#others (6)
#poids                            (bmxwt_bmx)
#dentiste                         (var_dentiste)
#carbohydrate                     (dr1tcarb_dr1tot)
#problemes coronaires             (mcq160c_mcq)
#folates                          (dr1tfola_dr1tot)
#selenium                         (dr1tsele_dr1tot)


#================================================

nhanes3fac<-nhanes3
nhanes3fac$DIQ010_diq<-factor(nhanes3fac$DIQ010_diq)


#tuning de la random forest
#--------------------------

# set.seed(1234)
# mod_complet <- randomForest(DIQ010_diq~.,data=nhanes3fac[,-1],ntree=1000)
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
#   mod_complet <- randomForest(DIQ010_diq~.,data=nhanes3fac[,-1],ntree=200,mtry=i)  
#   oob[[i]]<-mod_complet$err.rate[, 1][200]
#   }
# 
# oob 
# oob<-unlist(oob)
# plot(x=seq(2,10,1),y=oob,xlab="mtry")
# lines(x=seq(2,10,1),y=oob)
# #mtry=6 donne le meilleur résultat


#validation croisee avec le dataset complet
#------------------------------------------
#
bloc<-4
set.seed(1234)
ind <- sample(1:nrow(nhanes3)%%bloc+1)
RES<- data.frame(Y=nhanes3$DIQ010_diq,
                 logcomp=0,
                 logstep=0,
                 logred15=0,
                 foret=0,
                 ridge=0,
                 lasso=0,
                 elastic=0,
                 adaboost=0,
                 logiboost=0)


XX<-model.matrix(nhanes3$DIQ010_diq~.,data=nhanes3[,-1])

LAridge=list()

foreach (i = 1:bloc, .packages = c("gbm","glmnet","randomForest")) %dopar% {
  XXA <- XX[ind!=i,]
  YYA <- as.matrix(nhanes3[ind!=i,"DIQ010_diq"])
  # 
  #1-logistique compl?te
  #-------------------
  mod <- glm(DIQ010_diq~.,data=nhanes3[ind!=i,-1],family="binomial")
  RES[ind==i,"logcomp"] <- predict(mod,nhanes3[ind==i,-1],type="response")
  
  # #2a-logistique avec step
  # #--------------------
  RES[ind==i,"logstep"] <- predict(bestmodelfor,nhanes3[ind==i,-1],type="response")
  # 
  # #2b-logistique réduite 15
  # #-------------------
  mod <- glm(DIQ010_diq~RIDAGEYR_demo +DR1TSUGR_dr1tot+BPQ080_bpq+MCQ080_mcq+BPQ020_bpq+DR1TALCO_dr1tot+DR1TMOIS_dr1tot+
               RIAGENDR_demo+DR1.320Z_dr1tot+DR1TCHOL_dr1tot+DR1TPROT_dr1tot+INDFMPIR_demo+DR1TIRON_dr1tot+Var_TENSIONDI+DR1TCAFF_dr1tot,
             data=nhanes3[ind!=i,-1],family="binomial")
  RES[ind==i,"logred15"] <- predict(mod,nhanes3[ind==i,-1],type="response")
  # 
  # #3-foret
  # #-----
  mod <- randomForest(DIQ010_diq~.,data=nhanes3fac[ind!=i,-1],ntree=200,mtry=6)
  RES[ind==i,"foret"] <- predict(mod,nhanes3fac[ind==i,-1],type="prob")[,2]
  # 
  # #4-ridge
  # #-----
  tmp <- cv.glmnet(XXA,YYA,alpha=0,family="binomial")
  LAridge <- c(LAridge,tmp$lambda.min)
  mod <- glmnet(XXA,YYA,alpha=0,lambda=tmp$lambda.min,family="binomial")
  RES[ind==i,"ridge"] <- predict(mod,newx=XX[ind==i,],type="response")
  # 
  # #5-lasso
  # #-----
  mod <- cv.glmnet(XXA,YYA,alpha=1,family="binomial")
  RES[ind==i,"lasso"] <- predict(mod,newx=XX[ind==i,],lambda=mod$lambda.1se,type="response")
  # 
  # #6-elastic
  # #-------
  mod <- cv.glmnet(XXA,YYA,alpha=0.5,family="binomial")
  mod <- glmnet(XXA,YYA,alpha=.5,lambda=tmp$lambda.min,family="binomial")
  RES[ind==i,"elastic"] <- predict(mod,newx=XX[ind==i,],type="response")
  
  # #7-adaboost
  # #--------
  tmp <- gbm(DIQ010_diq~.,data = nhanes3[ind!=i,-1], distribution = "adaboost", interaction.depth = 2,
         shrinkage = 0.1,n.trees = 500)
  M <- gbm.perf(tmp)[1]
  mod <- gbm(DIQ010_diq~.,data = nhanes3[ind!=i,-1], distribution = "adaboost", interaction.depth = 2,
           shrinkage = 0.1,n.trees = M)
  RES[ind==i, "adaboost"] <- predict(mod, newdata=nhanes3[ind==i,-1], type = "response", n.trees = M)
  
  # #8-logiboost
  # #---------
  tmp <- gbm(DIQ010_diq~.,data=nhanes3[ind!=i,-1], distribution="bernoulli", interaction.depth = 2,
          shrinkage=0.1,n.trees=500)
  M <- gbm.perf(tmp)[1]
  mod <- gbm(DIQ010_diq~.,data=nhanes3[ind!=i,-1], distribution="bernoulli", interaction.depth = 2,
          shrinkage=0.1,n.trees=M)
  RES[ind==i, "logiboost"] <- predict(mod,newdata=nhanes3[ind==i,-1], type= "response", n.trees = M)
  # 
  }

str(RES)


#write.csv(RES,"res_dia.csv")


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
# RESp[logiboost>=seuil,logiboostp:=1,]
# RESp[logiboost<seuil,logiboostp:=0]
# 
# RESp<-data.frame(RESp)
# 
# str(RESp)

#METRIQUES
#---------

#courbes ROC et aires sous la courbe
#-----------------------------------
roclogcomp <- roc(RES[,1],RES[,2])
plot(roclogcomp)
roclogstep <- roc(RES[,1],RES[,3])
lines(roclogstep,lty=3,lwd=5,col="red")
roclogred15 <- roc(RES[,1],RES[,4])
lines(roclogred15,lwd=5,col="blue")
rocfor <- roc(RES[,1],RES[,5])
lines(rocfor,col="green")
rocrid<-roc(RES[,1],RES[,6])
lines(rocrid,col="blue")
roclas<-roc(RES[,1],RES[,7])
lines(roclas,col="orange")
rocela<-roc(RES[,1],RES[,8])
lines(rocela,col="brown")
rocada<-roc(RES[,1],RES[,9])
lines(rocada,col="purple")
roclogib<-roc(RES[,1],RES[,10])
lines(roclogib,col="dodger blue")

legend("bottomright", legend=c("logcomp","logstep","logred15","foret","ridge","lasso","elastic","adaboost","logiboost"),
       col=c("black", "red","blue","green","blue","orange","brown","purple","dodger blue"),lty=c(1,3,3,1,1,1,1,1,1),lwd=1, cex=1)


auc(RES[,1],RES[,2])
auc(RES[,1],RES[,3])
auc(RES[,1],RES[,4])
auc(RES[,1],RES[,5])
auc(RES[,1],RES[,6])
auc(RES[,1],RES[,7])
auc(RES[,1],RES[,8])
auc(RES[,1],RES[,9])
auc(RES[,1],RES[,10])

precision <- function(X,Y,seuil=input$seuilmod){
  Xc <- cut(X,breaks=c(0,seuil,1),labels=c(0,1),include.lowest=TRUE)
  round(sum(as.factor(Y)==Xc)/(sum(as.factor(Y)==Xc)+sum(as.factor(Y)!=Xc)),3)}
RES

monerreur <- function(X,Y,seuil=0.5){
  table(cut(X,breaks=c(0,seuil,1)),as.factor(Y))
}
823/(876+823)
monerreur(RES[,2],RES[,1])
apply(RES,2,monerreur,Y=RES[,1],seuil=0.4)
(218+183)/(218+627)

precision(RES[,5],RES[,1],0.5)
str(RES)
class(RES)
RES[RES$foret=="0",]
Xc <- cut(RES[,5],breaks=c(0,0.5,1),include.lowest=TRUE,labels=c(0,1))
class(Xc)
table(Xc,useNA="always")
?cut
nhanes3fac[1:54,]
# 
# 
# resultat<-data.frame(logcomp=0,logstep=0,foret=0,ridge=0,lasso=0,elastic=0,adaboost=0,logiboost=0)
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
# nhanes3<-data.table(nhanes3)
# nhanes3[,.N,]
# nhanes3[DIQ010_diq=="0",.N,]
# nhanes3[DIQ010_diq=="1",.N,]
# taux_positifs<-nhanes3[DIQ010_diq=="1",.N,]/nhanes3[,.N,]*100
# taux_positifs
# taux_negatifs<-nhanes3[DIQ010_diq=="0",.N,]/nhanes3[,.N,]*100
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
# colnames(taux_synthese)<-c("logcomp","logstep","foret","ridge","lasso","elastic","adaboost","logiboost")
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
# lines(seq(0,1,0.1),taux_synthese$logiboost,col="dodger blue",lwd=2)
# 
# legend("topright", legend=c("logcomp","logstep","foret","ridge","lasso","elastic","adaboost","logiboost"),
#        col=c("black", "red","green","blue","orange","brown","purple","dodger blue"),lwd=2,lty=c(1,3,1,1,1,1,1,1), cex=1)



#=========================================
#Traitement de l'importance des variables
#=========================================


library(data.table)
library(mice)
library(randomForest)
library(glmnet)
library(pROC)
library(gbm)
library(foreach)
library(ModelMetrics)
library(Metrics)
library(plyr)
library(caret)#pour l'utilisation des donnÃ©es Caret
library(dplyr)
library(e1071)

don <- read.csv("nhanes_diab_mice_apres.csv")
don$X <- NULL
don$SEQN<-NULL


XX <- as.matrix(model.matrix(~.,don)[,-ncol(model.matrix(~.,don))])
YY <- as.matrix(model.matrix(~.,don)[,ncol(model.matrix(~.,don))])

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
str(don)
# 1a)Importance variable pour le modele logistique
mod_log <- glm(DIQ010_diq~.,data=don,family="binomial")
varimplog <- variable_imp(x=mod_log,t=1,mot="log")
varimplog

# 1b)Importance variable pour le modele logistique step
mod_log <- glm(DIQ010_diq ~ RIDAGEYR_demo + MCQ080_mcq + BPQ020_bpq + 
                 BPQ080_bpq + DR1TSUGR_dr1tot + BMXWT_bmx + INDFMPIR_demo + 
                 DR1TCHOL_dr1tot + DR1TALCO_dr1tot + RIAGENDR_demo + Var_ACTIVITE + 
                 Var_DENTISTE + Var_TENSIONDI + SLQ050_slq + DR1TACAR_dr1tot + 
                 BMXBMI_bmx + Var_SITUATION + HOD050_hoq + DR1TMOIS_dr1tot + 
                 DR1.320Z_dr1tot + DR1TCAFF_dr1tot + DR1TPROT_dr1tot + DR1TIRON_dr1tot + 
                 MCQ160N_mcq + DR1TRET_dr1tot + MCQ160D_mcq + MCQ160C_mcq + 
                 DR1TSELE_dr1tot + DR1TFOLA_dr1tot,data=don,family="binomial")

varimplogstep <- variable_imp(x=mod_log,t=1,mot="log")
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


# 2)Importance variable pour le modele ridge
tmp <- cv.glmnet(XX,YY,alpha=0,family="binomial")
mod_ridge  <- glmnet(XX,YY,alpha=0,lambda=tmp$lambda.min, family="binomial")
varimpridge <- variable_imp(x=mod_ridge$beta,t=2,mot="ridge")

# 3)Importance variable pour le modele lasso (colnames(XX)[mod_lasso$beta@i])
tmp <- cv.glmnet(XX,YY, alpha=1, family="binomial")
mod_lasso <- glmnet(XX,YY,alpha=1, lambda =tmp$lambda.1se,family="binomial" )
varimplasso <- variable_imp(x=mod_lasso$beta,t=2, mot="lasso")

# 4)Importance variable pour le modele elastic
tmp <- cv.glmnet(XX,YY, alpha=0.5, family="binomial")
mod_elastic <- glmnet(XX,YY,alpha = 0.5, lambda = tmp$lambda.min, family="binomial")
varimpelastic <- variable_imp(mod_elastic$beta,t=2,mot = "elastic")

# 5)Importance variable pour le modele Foret
mod_foret <- randomForest(factor(DIQ010_diq)~., data = don, ntree=200,mtry=6)
varimpforet <- variable_imp(mod_foret,t=1,mot="foret")

# 6)Importance variable pour le modele adaboost
tmp <- gbm(as.numeric(DIQ010_diq)~.,data = don, distribution = "adaboost", interaction.depth = 2,
           shrinkage = 0.1,n.trees = 500)
M <- gbm.perf(tmp)[1]
mod_adaboost <- gbm(as.numeric(DIQ010_diq)~.,data = don, distribution = "adaboost", interaction.depth = 2,
                    shrinkage = 0.1,n.trees = M)
varimpada <- variable_imp(summary(mod_adaboost),t=3, mot="adaboost")

# 7)Importance variable pour le modele logiboost
tmp <- gbm(as.numeric(DIQ010_diq)~.,data=don, distribution="bernoulli", interaction.depth = 2,
           shrinkage=0.1,n.trees=500)
M <- gbm.perf(tmp)[1]
mod_logiboost <- gbm(as.numeric(DIQ010_diq)~.,data=don, distribution="bernoulli", interaction.depth = 2,
                     shrinkage=0.1,n.trees=M)
varimplogibo <- variable_imp(summary(mod_logiboost),t=3,mot="logiboost")

# # ?)Importance variable pour le modele SVM (je ne sais pas appliquer la feature selection)
# mod_svm <- svm(DIQ010_diq~.,data=don, kernel="linear",probability=T)
# tmp <- tune(svm,DIQ010_diq~.,data=don, kernel="linear",probability=T,range=list(cost=c(0.1,1,10)))
# mod <- tmp$best.model

# Croisement des tables d'importance des variables
choix_var <- varimplog %>%
  full_join(varimplogstep) %>%
  full_join(varimpridge) %>%
  full_join(varimplasso) %>%
  full_join(varimpelastic) %>%
  full_join(varimpforet) %>%
  full_join(varimpada) %>%
  full_join(varimplogibo)

choix_var <- cbind(choix_var[,1],choix_var[,c(which(grepl("^imp",names(choix_var))))])

choix_var <- as.data.frame(choix_var)
names(choix_var)[1] <- "variable"
write.csv2(choix_var,"choix_var_dia.csv",row.names = FALSE)






