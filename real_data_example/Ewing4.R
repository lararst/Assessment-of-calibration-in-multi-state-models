source("/z_load_packages.R")
data.mstate <- SetValid

#obtain data.raw
unique_ids <- unique(widedata$Patient_ID)
person_ids.devel <- unique_ids[1:300]
person_ids.valid <- unique_ids[301:630]
data.raw <- widedata[widedata$Patient_ID %in% person_ids.valid, ]

#censoring times and indicator variable
data.raw$dtcens <- data.raw$survtime
data.raw$dtcens_s <- 1 - data.raw$DEATHstat

### Add tmat to data.mstat attributes
attributes(data.mstate)[["trans"]] <- tmat

### Create a mapping from person_id to id (according to data.mstate) and create the id variable in pv.comb and data.raw
tp.predicted$id <- tp.predicted$person_id
data.raw$person_id <- tp.predicted$id
data.raw$id <- tp.predicted$id
#pv.comb$person_id <- data.raw$person_id
data.mstate$person_id <- data.mstate$id

#add estimated probabilities
data.raw <- merge(tp.predicted, data.raw, by.x = "id", by.y = "id")
colnames(data.raw)
colnames(data.raw)[2]<-"person_id"
data.raw<-data.raw[,-21]
data.raw <- data.raw %>%
  rename(pstate1 = pt.1,pstate2 = pt.2,pstate3 = pt.3)


#check time between progression events and death
datatimes<- data.raw[,1]
datatimes$times = NA 
datatimes$death = NA
datatimes$PROGR = NA 
for (i in 1:nrow(data.raw)) {
  if(data.raw$DEATHstat[i]==1 & data.raw$Progrstat[i]==1){
    datatimes$times[i]= data.raw$survtime[i] - data.raw$Progrtime[i] 
    datatimes$death[i] = data.raw$survtime[i]
    datatimes$PROGR[i] = data.raw$Progrtime[i] 
  }
}
summary(datatimes)
timestry<- seq(0, 6, by = 0.1)
for (t in timestry) {
  a=0
  a = sum(datatimes$PROGR <= t & datatimes$death > t, na.rm = TRUE)
  if(a>=30){
    print(t)
  }
}

################################################
### Calc calibration using BLR-IPCW approach ###
################################################
### Estimate mean calibration with BLR-IPCW
calib.mean.blr <- calib_msm(data_ms = SetValid,
                            data_raw = data.raw,
                            j = 1,
                            s = 0,
                            t = t.eval,
                            tp_pred = as.matrix(data.raw[, c("pstate1", "pstate2", "pstate3")]),
                            calib_type = "blr",
                            w_covs = covs,
                            assess_moderate = FALSE,
                            assess_mean = TRUE)[["mean"]] 
calib.mean.blr
### Estimate moderate calibration with BLR-IPCW
calib.mod.blr <- calib_msm(data_ms = SetValid,
                           data_raw = data.raw,
                           j=1,s=0,t=t.eval,
                           tp_pred = as.matrix(data.raw[, c("pstate1", "pstate2", "pstate3")]),
                           calib_type="blr",
                           curve_type="rcs",
                           rcs_nk=3,
                           w_covs=covs,
                           CI = 95,
                           CI_R_boot = 100,
                           CI_type = "bootstrap",
                           assess_moderate=TRUE,
                           assess_mean=FALSE)


################################################
### Calc calibration using MLR-IPCW approach ###
################################################
### Estimate mean calibration with MLR-IPCW
calib.mean.mlr <- calib_msm(data_ms = SetValid,
                            data_raw = data.raw,
                            j = 1,
                            s = 0,
                            t = t.eval,
                            tp_pred = as.matrix(data.raw[, c("pstate1", "pstate2", "pstate3")]),
                            calib_type = "mlr",
                            w_covs = covs,
                            assess_moderate = FALSE,
                            assess_mean = TRUE)[["mean"]] 

calib.mean.mlr
### Estimate moderate calibration with MLR-IPCW
calib.mod.mlr <- calib_msm(data_ms = SetValid,
                           data_raw = data.raw,
                           j=1,s=0,t=t.eval,
                           tp_pred = as.matrix(data.raw[, c("pstate1", "pstate2", "pstate3")]),
                           calib_type="mlr",
                           w_covs=covs,
                           assess_moderate=TRUE,
                           assess_mean=FALSE)


#######################################################
### Calc calibration using Aalen-Johansen estimator ###
#######################################################
### Estimate mean calibration with Aalen-Johansen
calib.aj.mean <- calib_msm(data_ms = SetValid,
                           data_raw = data.raw,
                           j = 1,
                           s = 0,
                           t = t.eval,
                           tp_pred = as.matrix(data.raw[, c("pstate1", "pstate2", "pstate3")]),
                           calib_type = "aj",
                           pv_n_pctls = 2,
                           assess_moderate = FALSE,
                           assess_mean = TRUE)[["mean"]]

calib.aj.mean

####################################################
### Calc calibration using pseudo-value approach ###
####################################################
### Estimate moderate calibration with pseudo-values
calib.mod.pv <- calib_msm(data_ms = SetValid,
                          data_raw = data.raw,
                          j=1,s=0,t=t.eval,
                          tp_pred = as.matrix(data.raw[, c("pstate1", "pstate2", "pstate3")]),
                          calib_type="pv",
                          curve_type="rcs",
                          rcs_nk=3,
                          pv_n_pctls = 2,
                          CI = 95,
                          CI_type = "parametric",
                          assess_moderate=TRUE,
                          assess_mean=FALSE)

################################################################
#CALIBRATION PLOTS
###############################################################

### Create plot objects
plot.pv <- plot(calib.mod.pv, nrow = 1, ncol = 3, size.text = 16, size.line = 0.5,
                marg.density = TRUE,
                axis.titles.x = 0, axis.titles.y = 1, axis.titles.text.y = "Pseudo-value\nObserved risk", 
                legend.include = TRUE, legend.seperate = TRUE, 
                titles = c("Rand", "Progr", "Death"))
plot.blr <- plot(calib.mod.blr, nrow = 1, ncol = 3, size.text = 16, size.line = 0.5,
                 marg.density = TRUE,
                 axis.titles.x = 0, axis.titles.y = 1, axis.titles.text.y = "BLR-IPCW\nObserved risk", 
                 legend.include = FALSE, titles.include = FALSE)
plot.mlr <- plot(calib.mod.mlr, nrow = 1, ncol = 3, transparency.plot = 0.05, size.text = 16, size.point = 0.75,
                 marg.density = TRUE,
                 axis.titles.x = 0, axis.titles.y = 1, axis.titles.text.y = "MLR-IPCW\nObserved risk", 
                 titles.include = FALSE)

png(
  filename = "/calib_pv.png",
  width = 3600,
  height = 1200,
  res = 300
)
grid::grid.draw(plot.pv)
dev.off()

png(
  filename = "/calib_blr.png",
  width = 3600,
  height = 1200,
  res = 300
)
grid::grid.draw(plot.blr)
dev.off()

png(
  filename = "/calib_mlr.png",
  width = 3600,
  height = 1200,
  res = 300
)
grid::grid.draw(plot.mlr)
dev.off()
