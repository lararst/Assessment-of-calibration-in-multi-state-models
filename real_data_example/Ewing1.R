##########################################################
#INSTALL PACKAGES
##########################################################

#install.packages("readxl")
library(readxl)
#install.packages("mstate")
library(mstate)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("colorspace")
library(colorspace)
#install.packages("dplyr")
library(dplyr)
#install.packages("arsenal")
library(arsenal)

#####################################################
#Importing all the files
#####################################################
baseline0 <- read_excel("/Baseline.xlsx", na = "NA")
death0 <- read_excel("/Death.xlsx", na = "NA")
followUp0 <- read_excel("/Follow Up.xlsx", na = "NA")
indChem0 <- read_excel("/InductionChemotherapy.xlsx", na = "NA")
patient0 <- read_excel("/Patient.xlsx", na = "NA")
progression0 <- read_excel("/Progression.xlsx", na = "NA")
radiotherapy0 <- read_excel("/Radiotherapy.xlsx", na = "NA")
surgery0 <- read_excel("/Surgery.xlsx", na = "NA")
chemo0 <- read_excel("/InductionChemotherapy.xlsx", na = "NA")

#save a copy
baseline <- baseline0
death <- death0
followUp <- followUp0
indChem <- indChem0
patient <- patient0
progression <- progression0
radiotherapy <- radiotherapy0
surgery <- surgery0
chemo <- chemo0

###################################################################
#FOLLOW UP AND DEATH
###################################################################
#removing all the lines where PatAlive is NA or the patient is alive but there is not a time value
followUp0 <- followUp0[
  !(is.na(followUp0$PatAlive) |
      (followUp0$PatAlive == "Yes" & is.na(followUp0$Time_to_last_seen))),
]

#I consider just three columns of the set followup
followUp1 <- followUp0[c(1,2,3)]
#now I want to obtain the last observation for each patient in the set
a <- c()
for(i in 1:(nrow(followUp1)-1) ) {
  if(followUp1$Patient_ID[i] == followUp1$Patient_ID[i+1]){
    a <- c(a,i) 
  }
}
followUp1 <- followUp1[-a, ] 
#Comparing the patients labeled as not alive in the set followUp and the ones listed in the set Death
a <- c()
for(j in 1:nrow(followUp1) ) {
  if (!is.na(followUp1$PatAlive[j]) & followUp1$PatAlive[j] == "No"){
    a <- c(a,j)
  }
}
#patients dead in FollowUp
deadFollowUp <- followUp1[a,]
#patients alive in FollowUp
aliveFollowUp <- followUp1[-a,]

#eliminating dead patients from aliveFollowUp
a <- c()
for(i in 1:nrow(aliveFollowUp)){
  if(aliveFollowUp$Patient_ID[i] %in% death0$Patient_ID){
    a <- c(a,i) 
  }
}
aliveFollowUp<-aliveFollowUp[-a,]

deadFollowUp$Patient_ID %in% death0$Patient_ID
#all the patient reported dead in FollowUp are considered in Death

#Total number of deaths:
nrow(death0)

########################################################################
#Some patients only appear at the beginning of the study, but never again in further evaluations
#These patients will not be considered anymore
########################################################################

a <- c()
for (i in 1:nrow(baseline0) ){
  if(baseline0$Patient_ID[i] %in% aliveFollowUp$Patient_ID | baseline0$Patient_ID[i] %in% death0$Patient_ID | baseline0$Patient_ID[i] %in% progression0$Patient_ID | baseline0$Patient_ID[i] %in% radiotherapy0$Patient_ID | baseline0$Patient_ID[i] %in% surgery0$Patient_ID){
    a <- c(a,i)
  }
}
lostpatients<-baseline0[-a,]

a <- c()
for (i in 1:nrow(baseline0) ){
  if(baseline0$Patient_ID[i] %in% lostpatients$Patient_ID){
    a <- c(a,i)
  }
}
baseline0<-baseline0[-a,]
patient0<-patient0[-a,]

######################################################################
#Obtaining categorical baseline variables
#New dataframe at baseline
######################################################################
databaseline<-subset(patient0,select = c(Patient_ID, PAT_Sex, PAT_R1TumourVolume))
for (i in 1:nrow(patient0) ) {
  if(patient0$PAT_Age[i]<14){
    databaseline$Pat_Age[i]="<14"
  }else{
    databaseline$Pat_Age[i]=">=14"
  }
}

a<-c()
for (i in 1:nrow(baseline0) ) {
  if(!is.na(baseline0$TumourBodyPart[i])){
    if(baseline0$TumourBodyPart[i]=="Lower extremity" | baseline0$TumourBodyPart[i]=="Upper
       extremity"){
      databaseline$TumourBodyPart[i]="Limb"
    }else if(baseline0$TumourBodyPart[i]=="Unknown origin"){
      a <- c(a,i)
    }else if(baseline0$TumourBodyPart[i]=="Pelvis"|baseline0$TumourBodyPart[i]=="Spine" | baseline0$TumourBodyPart[i]=="Abdomen" | baseline0$TumourBodyPart[i]=="Head and neck" | baseline0$TumourBodyPart[i]=="Chest"){
      databaseline$TumourBodyPart[i]="Axis"
    }
  }else{
    a <- c(a,i)
  }
}
deletedpat<-databaseline[a,]
deletedpat<-deletedpat$Patient_ID
databaseline<-databaseline[-a,]
baseline0<-baseline0[-a,]
print(table(baseline0[7]))
print(prop.table(table(baseline0[7])))

############################################################
#CHEMOTHERAPY
#brunch A (VIDE) and B (VDC, IE)
##############################################################

chemo1<-distinct(chemo0, Patient_ID, .keep_all = TRUE)
a<-c()
for (i in 1:nrow(chemo1) ) {
  if(chemo$Patient_ID[i]%in%deletedpat){
    a <- c(a,i)
  }
  if(!is.na(chemo1$VIDE[i]) & chemo1$VIDE[i]==1){
    databaseline$ChemoArm[i]="A"
  }
  if((!is.na(chemo1$VDC[i]) & chemo1$VDC[i]==1)|(!is.na(chemo1$IE[i]) & chemo1$IE[i]==1)){
    databaseline$ChemoArm[i]="B"
  }
}
chemo1<-chemo1[-a,]

################################################################
#Setting reference levels in baseline variables
###############################################################
databaseline$PAT_Sex <- as.factor(databaseline$PAT_Sex)
databaseline$PAT_Sex <- relevel(databaseline$PAT_Sex, ref = "Male")
databaseline$PAT_R1TumourVolume <- as.factor(databaseline$PAT_R1TumourVolume)
databaseline$PAT_R1TumourVolume <- relevel(databaseline$PAT_R1TumourVolume, ref = "<200ml")
databaseline$Pat_Age <- as.factor(databaseline$Pat_Age)
databaseline$Pat_Age <- relevel(databaseline$Pat_Age, ref = ">=14")
databaseline$TumourBodyPart <- as.factor(databaseline$TumourBodyPart)
databaseline$TumourBodyPart <- relevel(databaseline$TumourBodyPart, ref = "Axis")
databaseline$ChemoArm <- as.factor(databaseline$ChemoArm)
databaseline$ChemoArm <- relevel(databaseline$ChemoArm, ref = "B")

############################################################################
#CHECKING LAST OBSERVATION OF PATIENTS AND SAVING TIMES OF EVENTS OF INTEREST
#############################################################################
baseline0<-baseline0[!(baseline0$Patient_ID%in%deletedpat),]
alltimes<-baseline0[c(1)]
for (i in 1:nrow(alltimes) ){
  alltimes$survtime[i] = NA
  alltimes$lastfollowup[i] = NA
  for (j in 1:nrow(death0) ){
    if (alltimes$Patient_ID[i] == death0$Patient_ID[j] ){
      alltimes$survtime[i] = (death0$Time_to_death[j])
    }
  }
  for (j in 1:nrow(aliveFollowUp) ){
    if (alltimes$Patient_ID[i] == aliveFollowUp$Patient_ID[j] ){
      alltimes$lastfollowup[i] = (aliveFollowUp$Time_to_last_seen[j])
    }
  }
}

#####################################################################
#PROGRESSION
#####################################################################
progression0<-progression0[(progression0$Patient_ID%in%alltimes$Patient_ID),]
#keeping only the first line of every patient
progression1<-progression0
progression0<-distinct(progression0, Patient_ID, .keep_all = TRUE)

#all the metastasis that are not pulmonary
progrMetast <- subset(progression0, select = c(Bone, BoneMarrow, CNS, Liver, LymphNode, OtherMetastas))

for (i in 1:nrow(progression0) ){
  progression0$LRstat[i]=0
  progression0$DMPstat[i]=0
  progression0$DMOstat[i]=0
  progression0$DEATHstat[i]=0
  
  if(!is.na(progression0$ProgressionType[i]) & progression0$ProgressionType[i]=="Local"){
    progression0$LRstat[i]=1
  }
  if(!is.na(progression0$ProgressionType[i]) & progression0$ProgressionType[i]=="Metastatic"){
    if(!is.na(progression0$LungPleura[i]) & progression0$LungPleura[i]=="Yes" ){
      if("Yes" %in% progrMetast[i,]){
        progression0$DMOstat[i]=1
      }else{
        progression0$DMPstat[i]=1
      }
    }else{
      progression0$DMOstat[i]=1
    }
  }
}
Mod1 <- progression0[,c(1,3,15,16,17,18)]

#######################################################################
#wide format
#######################################################################

wideMod1 <- alltimes
for  (i in 1:nrow(wideMod1) ){
  wideMod1$LRstat[i]=0
  wideMod1$DMPstat[i]=0
  wideMod1$DMOstat[i]=0
  wideMod1$DEATHstat[i]=0
  
  wideMod1$LRtime[i]=NA
  wideMod1$DMPtime[i]=NA
  wideMod1$DMOtime[i]=NA
  
  #If a patient died the status of death is 1
  if(!is.na(wideMod1$survtime[i])){
    wideMod1$DEATHstat[i]=1
  }
  
  #saving other states and times of events, if an event time equals the survtime, there is actually no progression to that event
  for (j in 1:nrow(Mod1) ){
    if(wideMod1$Patient_ID[i]==Mod1$Patient_ID[j]){
      if(Mod1$LRstat[j]==1){
        wideMod1$LRstat[i]=1
        wideMod1$LRtime[i]=Mod1$Time_to_prog_confirmed[j]
        if(!is.na(wideMod1$LRtime[i]) & !is.na(wideMod1$survtime[i]) & (wideMod1$LRtime[i]==wideMod1$survtime[i])){
          wideMod1$LRtime[i] = wideMod1$survtime[i] - 1
        }
      }
      if(Mod1$DMPstat[j]==1){
        wideMod1$DMPstat[i]=1
        wideMod1$DMPtime[i]=Mod1$Time_to_prog_confirmed[j]
        if(!is.na(wideMod1$DMPtime[i]) & !is.na(wideMod1$survtime[i])  & wideMod1$DMPtime[i]==wideMod1$survtime[i]){
          wideMod1$DMPtime[i] = wideMod1$survtime[i] - 1
        }
      }
      if(Mod1$DMOstat[j]==1){
        wideMod1$DMOstat[i]=1
        wideMod1$DMOtime[i]=Mod1$Time_to_prog_confirmed[j]
        if(!is.na(wideMod1$DMOtime[i]) & !is.na(wideMod1$survtime[i])  & wideMod1$DMOtime[i]==wideMod1$survtime[i]){
          wideMod1$DMOtime[i] = wideMod1$survtime[i] - 1
        }
      }
    }
  }
}

#adjustung na times: na values are substituted by the time of death (survtime) 
# or, if they did not die during the study, the time of the last followup or the 
#higher time event recorded.
oldwide<-wideMod1
a<-c()
for  (i in 1:nrow(wideMod1) ){
  if(!is.na(wideMod1$survtime[i])){
    if(is.na(wideMod1$LRtime[i])){
      wideMod1$LRtime[i]=wideMod1$survtime[i]
    }
    if(is.na(wideMod1$DMPtime[i])){
      wideMod1$DMPtime[i]=wideMod1$survtime[i]
    }
    if(is.na(wideMod1$DMOtime[i])){
      wideMod1$DMOtime[i]=wideMod1$survtime[i]
    }
  }else{
    if(is.na(wideMod1$lastfollowup[i]) & is.na(wideMod1$LRtime[i]) & 
       is.na(wideMod1$DMPtime[i]) & is.na( wideMod1$DMOtime[i])){
      a <- c(a,i)
    }else{
      tmax=max(wideMod1$lastfollowup[i], wideMod1$LRtime[i], wideMod1$DMPtime[i], wideMod1$DMOtime[i],  na.rm=T)
      if(is.na(wideMod1$LRtime[i])){
        wideMod1$LRtime[i]=tmax
      }
      if(is.na(wideMod1$DMPtime[i])){
        wideMod1$DMPtime[i]=tmax
      }
      if(is.na(wideMod1$DMOtime[i])){
        wideMod1$DMOtime[i]=tmax
      }
      if(is.na(wideMod1$survtime[i])){
        wideMod1$survtime[i]=tmax
      }
    }
  }
}
a

#eliminate patients that had no events and no follow-up
wideMod1<-wideMod1[-a,]

#removing lastfollowup
wideMod1 <- subset(wideMod1, select = -lastfollowup) 

##################################################################
#Adding baseline variables to the wide format dataset
##################################################################
baselinevar <- databaseline[(databaseline$Patient_ID%in%wideMod1$Patient_ID),]
widedata <- cbind(wideMod1,baselinevar$PAT_Sex,baselinevar$PAT_R1TumourVolume,
                  baselinevar$Pat_Age,baselinevar$TumourBodyPart,baselinevar$ChemoArm) 
colnames(widedata)[14]<-"ChemoArm"
colnames(widedata)[13]<-"TumourBodyPart"
colnames(widedata)[12]<-"Pat_Age"
colnames(widedata)[11]<-"TumourVolume"
colnames(widedata)[10]<-"PAT_Sex"

#time calculated in years
widedata$LRtime<-widedata$LRtime/365
widedata$DMPtime<-widedata$DMPtime/365
widedata$DMOtime<-widedata$DMOtime/365
widedata$survtime<-widedata$survtime/365

#######################################################################
#Function that checks if t-start and t-stop are correct
#######################################################################
check.start.stop <- function(data){
  #tstop cannot be smaller than tstart
  for  (i in 1:nrow(data)){
    if(data$Tstart[i]>=data$Tstop[i]){
      data$Tstop[i]=data$Tstart[i]+0.001
      data$time[i]=data$Tstop[i]-data$Tstart[i]
      print(i)
      #adjusting following observations of the same patient
      for (j in (i+1):nrow(data)){
        if (data$id[j]==data$id[i] & data$from[j]!=data$from[i]){
          if(data$Tstart[j]<data$Tstop[i]){
            data$Tstart[j]=data$Tstop[i]
          }
        }
      }
    }
    #check on infinite values
    if(is.infinite(data$Tstop[i])){
      data$Tstop[i]=data$Tstart[i]+0.001
      data$time[i]=data$Tstop[i]-data$Tstart[i]
      print(i)
      for (j in (i+1):nrow(data)){
        if (data$id[j]==data$id[i] & data$from[j]!=data$from[i]){
          if(data$Tstart[j]<data$Tstop[i]){
            data$Tstart[j]=data$Tstop[i]
          }
        }
      }
    }
  }
  return(data)
}

#################################################################
#FREQUENCY TABLES
#################################################################
for (i in c(10:14) ){
  print(table(widedata[i]))
  print(prop.table(table(widedata[i])))
}

######################################################################
#Frequency tables baseline under different treatment
#####################################################################
chemA<-widedata[widedata$ChemoArm=="A",]
chemB<-widedata[widedata$ChemoArm=="B",]
chemA<-chemA[,-14]
chemB<-chemB[,-14]

for (i in c(10:13) ){
  print(table(chemA[i]))
  print(prop.table(table(chemA[i])))
}
for (i in c(10:13) ){
  print(table(chemB[i]))
  print(prop.table(table(chemB[i])))
}

#SURVIVAL ANALYSIS
#########################################################################
#transition matrix
tmat <- transMat(x = list(c(2,3,4,5), 
                          c(5),
                          c(5),
                          c(5),
                          c()),
                 names = c("Rand", "LR", "DMP", "DMO","Death"))
tmat

#long format
covs <- c("PAT_Sex", "TumourVolume","Pat_Age", "TumourBodyPart","ChemoArm")
longdata <- msprep(time = c(NA, "LRtime", "DMPtime", "DMOtime","survtime"), 
                   status = c(NA, "LRstat", "DMPstat", "DMOstat", "DEATHstat"), 
                   data = widedata, trans = tmat, keep = covs)
longdata = check.start.stop(longdata)

#Cox model with no covariates
CoxNoCov <- coxph(Surv(Tstart, Tstop, status ) ~ strata(trans), data =longdata, method = "breslow")
CoxNoCov
#number of events
events(longdata)

longdataCovs <- expand.covs(longdata, covs, longnames = FALSE)

colnames(longdataCovs)

#MODEL WITH ALL TRANS-SPEC COV
cfullMod <- coxph(Surv(Tstart, Tstop, status ) ~PAT_Sex.1 + PAT_Sex.2 +PAT_Sex.3 
                  + PAT_Sex.4 +PAT_Sex.5 +PAT_Sex.6 +PAT_Sex.7 + TumourVolume.1 
                  +TumourVolume.2+TumourVolume.3+TumourVolume.4+TumourVolume.5
                  +TumourVolume.6 + TumourVolume.7 + Pat_Age.1+ Pat_Age.2+ Pat_Age.3
                  +Pat_Age.4+ Pat_Age.5+ Pat_Age.6+ Pat_Age.7 +TumourBodyPart.1+TumourBodyPart.2
                  +TumourBodyPart.3+TumourBodyPart.4+TumourBodyPart.5+TumourBodyPart.6 
                  +TumourBodyPart.7 + ChemoArm.1+ChemoArm.2+ChemoArm.3+ChemoArm.4
                  +ChemoArm.5 +ChemoArm.6 +ChemoArm.7 +
                    + strata(trans), data =longdataCovs, method = "breslow")
cfullMod


#final model
CoxModelFinal <- coxph(Surv(Tstart, Tstop, status ) ~ PAT_Sex.2 +PAT_Sex.3 +PAT_Sex.5
                       +PAT_Sex.6 +PAT_Sex.7 + TumourVolume.1 +TumourVolume.2+TumourVolume.3
                       +TumourVolume.5 +TumourVolume.6 +TumourVolume.7 + Pat_Age.2+ Pat_Age.3 
                       +Pat_Age.4 +Pat_Age.6 +Pat_Age.7 +TumourBodyPart.1 +TumourBodyPart.2 
                       +TumourBodyPart.3 +TumourBodyPart.5 +TumourBodyPart.6 +TumourBodyPart.7
                       +ChemoArm.1 +ChemoArm.2 +ChemoArm.3 +ChemoArm.6 +ChemoArm.7
                       + strata(trans), data =longdataCovs, method = "breslow")
CoxModelFinal

#plot predicted probabilities
msf0 <- msfit(CoxNoCov, trans = tmat)
pt0  <- probtrans(msf0, predt = 0)

png("/PredProb_plot.png",
    width = 2000, height = 1600, res = 300)
plot(pt0, ord=c(2,3,4,5,1),
     #main=list("Overall population", cex=1.5),
     las=1,
     xlab=list("Years since randomization", cex=1.3),
     xlim=c(0,6),
     type="filled",
     col=c("red2","orange","lightyellow","lightgreen","#00B81F"))
dev.off()



#predictions on patients with different chemo
whA <- which(longdata$PAT_Sex=="Male" & longdata$TumourVolume=="<200ml" &
               longdata$Pat_Age==">=14" & longdata$TumourBodyPart=="Axis" & longdata$ChemoArm=="A")

patA <- longdata[rep(whA[1],7), 9:13]
patA$trans <- 1:7
attr(patA, "trans") <- tmat
patA <- expand.covs(patA, covs, longnames = FALSE)
patA$strata <- patA$trans
msfA <- msfit(CoxModelFinal, patA, trans = tmat)
ptA <- probtrans(msfA, predt = 0)
#msfA0 <- msfit(CoxNoCov, patA, trans = tmat)
#ptA0 <- probtrans(msfA0, predt = 0)
#rowSums(ptA0[[1]][,2:6])

#save graph
png("/PatientA_plot.png",
    width = 2000, height = 1600, res = 300)
plot(ptA, ord=c(2,3,4,5,1), main=list("Patient A", cex=1.5), las=1, 
     xlab=list("Years since randomization", cex=1.5), xlim=c(0,6), type="filled", 
     col= c("red2", "orange", "light yellow", "light green", "#00B81F"))
dev.off()

#patB
whB <- which(longdata$PAT_Sex=="Male" & longdata$TumourVolume=="<200ml" &
               longdata$Pat_Age==">=14" & longdata$TumourBodyPart=="Axis" & longdata$ChemoArm=="B")

patB <- longdata[rep(whB[1],7), 9:13]
patB$trans <- 1:7
attr(patB, "trans") <- tmat
patB <- expand.covs(patB, covs, longnames = FALSE)
patB$strata <- patB$trans
msfB <- msfit(CoxModelFinal, patB, trans = tmat)
ptB <- probtrans(msfB, predt = 0)

#save graph
png("/PatientB_plot.png",
    width = 2000, height = 1600, res = 300)
plot(ptA, ord=c(2,3,4,5,1), main=list("Patient B", cex=1.5), las=1, 
     xlab=list("Years since randomization", cex=1.5), xlim=c(0,6), type="filled", 
     col= c("red2", "orange", "light yellow", "light green", "#00B81F"))
dev.off()


#predictions on patients with different tumor volume
whA <- which(longdata$PAT_Sex=="Male" & longdata$TumourVolume=="<200ml" &
               longdata$Pat_Age==">=14" & longdata$TumourBodyPart=="Axis" & longdata$ChemoArm=="A")

patA <- longdata[rep(whA[1],7), 9:13]
patA$trans <- 1:7
attr(patA, "trans") <- tmat
patA <- expand.covs(patA, covs, longnames = FALSE)
patA$strata <- patA$trans
msfA <- msfit(CoxModelFinal, patA, trans = tmat)
ptA <- probtrans(msfA, predt = 0)

#save graph
png("/PatientAtv_plot.png",
    width = 2000, height = 1600, res = 300)
plot(ptA, ord=c(2,3,4,5,1), main=list("Patient A", cex=1.5), las=1, 
     xlab=list("Years since randomization", cex=1.5), xlim=c(0,6), type="filled", 
     col= c("red2", "orange", "light yellow", "light green", "#00B81F"))
dev.off()

#patB
whB <- which(longdata$PAT_Sex=="Male" & longdata$TumourVolume==">=200ml" &
               longdata$Pat_Age==">=14" & longdata$TumourBodyPart=="Axis" & longdata$ChemoArm=="A")

patB <- longdata[rep(whB[1],7), 9:13]
patB$trans <- 1:7
attr(patB, "trans") <- tmat
patB <- expand.covs(patB, covs, longnames = FALSE)
patB$strata <- patB$trans
msfB <- msfit(CoxModelFinal, patB, trans = tmat)
ptB <- probtrans(msfB, predt = 0)

#save graph
png("/PatientBtv_plot.png",
    width = 2000, height = 1600, res = 300)
plot(ptA, ord=c(2,3,4,5,1), main=list("Patient B", cex=1.5), las=1, 
     xlab=list("Years since randomization", cex=1.5), xlim=c(0,6), type="filled", 
     col= c("red2", "orange", "light yellow", "light green", "#00B81F"))
dev.off()

#predictions on patients with different age
whA <- which(longdata$PAT_Sex=="Male" & longdata$TumourVolume=="<200ml" &
               longdata$Pat_Age==">=14" & longdata$TumourBodyPart=="Axis" & longdata$ChemoArm=="A")

patA <- longdata[rep(whA[1],7), 9:13]
patA$trans <- 1:7
attr(patA, "trans") <- tmat
patA <- expand.covs(patA, covs, longnames = FALSE)
patA$strata <- patA$trans
msfA <- msfit(CoxModelFinal, patA, trans = tmat)
ptA <- probtrans(msfA, predt = 0)

#save graph
png("/PatientAtageplot.png",
    width = 2000, height = 1600, res = 300)
plot(ptA, ord=c(2,3,4,5,1), main=list("Patient A", cex=1.5), las=1, 
     xlab=list("Years since randomization", cex=1.5), xlim=c(0,6), type="filled", 
     col= c("red2", "orange", "light yellow", "light green", "#00B81F"))
dev.off()

#patB
whB <- which(longdata$PAT_Sex=="Male" & longdata$TumourVolume=="<200ml" &
               longdata$Pat_Age=="<14" & longdata$TumourBodyPart=="Axis" & longdata$ChemoArm=="A")

patB <- longdata[rep(whB[1],7), 9:13]
patB$trans <- 1:7
attr(patB, "trans") <- tmat
patB <- expand.covs(patB, covs, longnames = FALSE)
patB$strata <- patB$trans
msfB <- msfit(CoxModelFinal, patB, trans = tmat)
ptB <- probtrans(msfB, predt = 0)

#save graph
png("/PatientBage_plot.png",
    width = 2000, height = 1600, res = 300)
plot(ptA, ord=c(2,3,4,5,1), main=list("Patient B", cex=1.5), las=1, 
     xlab=list("Years since randomization", cex=1.5), xlim=c(0,6), type="filled", 
     col= c("red2", "orange", "light yellow", "light green", "#00B81F"))
dev.off()






