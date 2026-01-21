###
### This program combines the results from the small sample simulation for moderate calibration, ready to be plotted.
###
scen <- "C1" #C2, C3
n.cohort  <- 700 #500, 400
n.pctls <- 10

### Print these
print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))

### Create list to store calibration of each estimates
calib.pv.mod.list.comb <- vector("list", 3)
calib.blr.mod.list.comb <- vector("list", 3)
calib.mlr.mod.list.comb <- vector("list", 3)

### Combine them
for (set in 1:50){
  
  ### Load data
  load(paste("/small_sample_moderate_analysis_", scen, "_n", n.cohort, "_npctls", n.pctls, "_set", set, ".RData", sep = ""))
  
  ### Concatenate
  calib.pv.mod.list.comb <- lapply(c(1,2,3), function(x) {c(calib.pv.mod.list.comb[[x]], calib.pv.mod.list[[x]])})
  calib.blr.mod.list.comb <- lapply(c(1,2,3), function(x) {c(calib.blr.mod.list.comb[[x]], calib.blr.mod.list[[x]])})
  calib.mlr.mod.list.comb <- lapply(c(1,2,3), function(x) {c(calib.mlr.mod.list.comb[[x]], calib.mlr.mod.list[[x]])})
  
}

### Save image
rm(list=setdiff(ls(), list("calib.pv.mod.list.comb", "calib.blr.mod.list.comb",  "calib.mlr.mod.list.comb", 
                           "n.cohort", "scen", "n.pctls")))
save.image(paste("/small_sample_moderate_analysis_", scen, "_n", n.cohort, "_npctls", n.pctls, ".RData", sep = ""))
print("IMAGE SAVED")
