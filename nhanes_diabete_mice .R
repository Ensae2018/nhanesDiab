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

write.csv(nhanes1,"nhanes_dia_avant_mice.csv")

dim(nhanes1)#5221x85


class(nhanes1$SLD012_slq)
summary(nhanes1$SLD012_slq)



#2- Imputation par mice
#----------------------


imp<-mice(nhanes1,m=1)

nhanes2<-complete(imp)
dim(nhanes2)#5221x85


dt<-data.table(var=names(nhanes2),nbna=apply(nhanes2,2,function (x) {sum(is.na(x))}))
dt[,sum(is.na(.SD)),]

str(nhanes2)

write.csv(nhanes2,"nhanes_dia_avant_transco.csv")

#3- Renomnage des colonnes
#-------------------------

nhanes_transco<-read.csv("nhanes_dia_transcodifie.csv")
str(nhanes_transco)
nhanes_transco$X<-NULL
dim(nhanes_transco)
dput(names(nhanes_transco))


str(nhanes2)

# #vérification de la non déformation de la distribution des variables
# str(nhanes2)
# barplot(table(nhanes2$Var_ASSURE,useNA="always"))
# barplot(table(nhanes1$Var_ASSURE,useNA="always"))
# table(nhanes1$Var_TRAVAIL,useNA="always")

names(nhanes2)<-c("SEQN","Gender", "Age.in.years.at.screening.", "Total.number.of.people.in.the.Household", 
                   "Annual.household.income", "Ratio.of.family.income.to.poverty", 
                   "Weight..kg.", "Standing.Height..cm.", "Body.Mass.Index..kg.m..2.", 
                   "Ever.told.you.had.high.blood.pressure", "high.cholesterol.level", 
                   "Doctor.told.you.have.diabetes", "Ever.told.you.have.Hepatitis.B.", 
                   "Ever.told.you.have.Hepatitis.C.", "Number.of.rooms.in.home", 
                   "Doctor.ever.said.you.were.overweight", "Doctor.ever.said.you.had.arthritis", "Doctor.ever.told.you.that.you.had.gout", 
                   "Ever.told.had.congestive.heart.failure", "Ever.told.you.had.coronary.hart.disease", "Ever.told.you.had.angina.pectoris", "Ever.told.you.had.heart.attack", "Ever.told.you.had.a .stoke", 
                   "Ever.told.you.had.emphysema", "Ever.told.you.had.thyroid.problem", "Sleep.hours.", "Ever.told.doctor.had.trouble.sleeping.", 
                   "Energy..kcal.", "Protein..gm.", "Carbohydrate..gm.", "Total.sugars..gm.", 
                   "Dietary.fiber..gm.", "Total.fat..gm.", "Total.saturated.fatty.acids..gm.", 
                   "Total.monounsaturated.fatty.acids..gm.", "Total.polyunsaturated.fatty.acids..gm.", 
                   "Cholesterol..mg.", "tocopherol..mg.", "tocopherol..Vitamin.E...mg.", 
                   "Retinol..mcg.", "Vitamin.A..RAE..mcg.", "carotene..mcg.", "carotene..mcg..1", 
                   "cryptoxanthin..mcg.", "Lycopene..mcg.", "Lutein...zeaxanthin..mcg.", 
                   "Thiamin..Vitamin.B1...mg.", "Riboflavin..Vitamin.B2...mg.", 
                   "Niacin..mg.", "Vitamin.B6..mg.", "Total.folate..mcg.", "Folic.acid..mcg.", 
                   "Food.folate..mcg.", "Folate..DFE..mcg.", "Total.choline..mg..", 
                   "Vitamin.B12..mcg.", "Added.vitamin.B12..mcg.", "Vitamin.C..mg.", 
                   "Vitamin.D..D2...D3...mcg.", "Vitamin.K..mcg.", "Calcium..mg.", 
                   "Phosphorus..mg.", "Magnesium..mg.", "Iron..mg.", "Zinc..mg.", 
                   "Copper..mg.", "Sodium..mg.", "Potassium..mg.", "Selenium..mcg.", 
                   "Caffeine..mg.", "Theobromine..mg.", "Alcohol..gm.", "Moisture..gm.", 
                   "DR1.320Z_DR1TOT", "VAR_EDUCATION", "VAR_TRAVAIL", "VAR_DEPRESSION", 
                   "VAR_SITUATION", "VAR_COFUMEUR", "VAR_TENSIONSY", "VAR_TENSIONDI", 
                   "VAR_ARGENTALIM", "VAR_DENTISTE", "VAR_ASSURE", "VAR_ACTIVITE")

write.csv(nhanes2,"nhanes_dia_apres_transco.csv")
