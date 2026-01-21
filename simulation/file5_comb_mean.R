###
### This program combines the results from the small sample simulation for mean calibration, ready to be plotted/tabled
###

### Set parameters
scen <- "C1"  #C2/C3
n.pctls <- 10
n.cohort  <- 700  #500/400
### Print these
print(paste("scen = ", scen, sep = ""))
print(paste("n.cohort = ", n.cohort, sep = ""))
print(paste("n.pctls = ", n.pctls, sep = ""))

### Create list to store calibration of each estimates
calib.aj.mean.list.comb <- vector("list", 3)
calib.blr.mean.list.comb <- vector("list", 3)
calib.mlr.mean.list.comb <- vector("list", 3)
calib.true.mean.list.comb <- vector("list", 3)

### Combine them
for (set in 1:50){
  
  ### Load data
  load(paste("/small_sample_mean_analysis_", scen, "_n", n.cohort, "_npctls", n.pctls, "_set", set, ".RData", sep = ""))
  
  ### Concatenate into dataset
  calib.aj.mean.list.comb <- lapply(c(1,2,3), function(x) {rbind(calib.aj.mean.list.comb[[x]], calib.aj.mean.list[[x]])})
  calib.blr.mean.list.comb <- lapply(c(1,2,3), function(x) {rbind(calib.blr.mean.list.comb[[x]], calib.blr.mean.list[[x]])})
  calib.mlr.mean.list.comb <- lapply(c(1,2,3), function(x) {rbind(calib.mlr.mean.list.comb[[x]], calib.mlr.mean.list[[x]])})
  calib.true.mean.list.comb <- lapply(c(1,2,3), function(x) {rbind(calib.true.mean.list.comb[[x]], calib.true.mean.list[[x]])})
  
}

### Save image
rm(list=setdiff(ls(), list("calib.aj.mean.list.comb", "calib.blr.mean.list.comb",  "calib.mlr.mean.list.comb", "calib.true.mean.list.comb", 
                           "n.cohort", "scen", "n.pctls")))
save.image(paste("/small_sample_mean_analysis_", scen, "_n", n.cohort, "_npctls", n.pctls, ".RData", sep = ""))
print("IMAGE SAVED")