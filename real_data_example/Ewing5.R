#this file produces predicted transition probabilities

#########################################################################
#SEPARATING SETS (development and validation)
unique_ids <- unique(longdata$id)
ids.devel <- unique_ids[1:300]
ids.valid <- unique_ids[301:630]
SetDevel <- longdata[longdata$id %in% ids.devel,]
events(SetDevel)
SetValid <- longdata[longdata$id %in% ids.valid,]
events(SetValid)
### Set t.eval
t.eval <- 2.7

Cox01 <- coxph(Surv(Tstart, Tstop, status ) ~ strata(trans), data =SetDevel, method = "breslow")
Cox01

Cox02 <- coxph(Surv(Tstart, Tstop, status ) ~ strata(trans), data =SetValid, method = "breslow")
Cox02

#expand covs
SetDevel <- expand.covs(SetDevel, covs, longnames = FALSE)
SetValid <- expand.covs(SetValid, covs, longnames = FALSE)

#Cox model on development set

CoxDev0 <- coxph(Surv(Tstart, Tstop, status ) ~PAT_Sex.1 + PAT_Sex.2 
                 + TumourVolume.1 +TumourVolume.2 + Pat_Age.1+ Pat_Age.2
                 + TumourBodyPart.1+TumourBodyPart.2 + ChemoArm.1+ChemoArm.2
                 + strata(trans), data =SetDevel, method = "breslow")
CoxDev0


### Add a column called strata (required when making predictors for new patients
SetDevel$strata <- SetDevel$trans
SetValid$strata <- SetValid$trans

###########################################################################
#estimated probabilities in the validation set
ce_calc_probtran <- data.frame(matrix(vector(), 0, 7,
                                      dimnames=list(c(), c("person_id", "pt.1", "pt.2",
                                                           "pt.3", "se.1", "se.2", "se.3"))),
                               stringsAsFactors=FALSE)

for (id in ids.valid) {
  ### Paste progress
  print(paste("id =", id, Sys.time(), sep = " "))
  ### Define person_id
  person_id_temp <- id
  
  ### Create temporary dataset to make predictors with
  ## Get location of individual
  person_id_loc <- which(SetValid$id == person_id_temp)
  
  ## Extract dataset with rows of this individual
  temp.data <- SetValid[rep(person_id_loc[1], max(tmat[!is.na(tmat)])), c("id", covs)]
  ## Add trans
  temp.data$trans <- 1:max(tmat[!is.na(tmat)])
  ## Attribute transition matrix
  attr(temp.data, "trans") <- tmat
  ## Expand the covariates
  temp.data <- expand.covs(temp.data, covs, longnames = FALSE)
  ## Add strata
  temp.data$strata <- temp.data$trans
  print(Sys.time())
  
  ### Calculate transition specific hazards
  msm.fit <- msfit(object = CoxDev0, trans = tmat, 
                   newdata = temp.data)
  print(Sys.time())
  
  ### Calculate probtrans
  probtrans <- probtrans(msm.fit, predt = 0)
  print(Sys.time())
  
  ### Extract prob trans
  pt.se <- probtrans[[1]] %>% 
    ## Extract row corresponding to crrecot time point
    slice(min(which(probtrans[[1]]$time > (t.eval)))) %>% 
    ## Extract all the correct columns
    select(c(paste("pstate",1:3, sep = ""), paste("se", 1:3, sep = "")))
  
  ### Create the row to add by adding patid to pt.se
  pt.se <- data.frame("person_id" = person_id_temp, pt.se)
  colnames(pt.se) <- c("person_id", paste("pt.",1:3, sep = ""), paste("se.", 1:3, sep = ""))
  
  ### Subset the validation cohort
  ce_calc_probtran= rbind(ce_calc_probtran, pt.se[1,])
}


tp.predicted <- ce_calc_probtran
