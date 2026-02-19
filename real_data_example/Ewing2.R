#this file checks if there are at least 30 patients in each state at some time points

#check time between progression events and death
datatimes<- widedata[,c(1,3,4,5,6)]
for (i in 1:nrow(widedata)) {
  datatimes$times[i] = NA 
  datatimes$death[i] = NA
  datatimes$PROGR[i] = NA 
  datatimes$DMO[i] = NA 
  if(widedata$DEATHstat[i]==1 & widedata$LRstat[i]==1){
    datatimes$times[i]= widedata$survtime[i] - widedata$LRtime[i] 
    datatimes$death[i] = widedata$survtime[i]
    datatimes$PROGR[i] = widedata$LRtime[i] 
  }
  if(widedata$DEATHstat[i]==1 & widedata$DMPstat[i]==1){
    datatimes$times[i]= widedata$survtime[i] - widedata$DMPtime[i] 
    datatimes$death[i] = widedata$survtime[i]
    datatimes$PROGR[i] = widedata$DMPtime[i] 
  }
  if(widedata$DEATHstat[i]==1 & widedata$DMOstat[i]==1){
    datatimes$times[i]= widedata$survtime[i] - widedata$DMOtime[i] 
    datatimes$death[i] = widedata$survtime[i]
    datatimes$DMO[i] = widedata$DMOtime[i] 
  }
}
summary(datatimes)
timestry<- seq(0, 6, by = 0.1)
for (t in timestry) {
  a=0
  b=0
  a = sum(datatimes$PROGR <= t & datatimes$death > t, na.rm = TRUE)
  b = sum(datatimes$DMO <= t & datatimes$death > t, na.rm = TRUE)
  if(a>=30 & b>=30){
    print(t)
  }
}

datatimesVal <- datatimes[301:630,]

summary(datatimesVal)

for (t in timestry) {
  a=0
  b=0
  a = sum(datatimesVal$PROGR <= t & datatimesVal$death > t, na.rm = TRUE)
  b = sum(datatimesVal$DMO <= t & datatimesVal$death > t, na.rm = TRUE)
  if(a>=30 & b>=30){
    print(t)
  }
}

