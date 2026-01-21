###
### This program runs the assessment of moderate calibration in the small sample simulation
###library(foreach)
library(doParallel)

source("/z_functions.R")
source("/z_functions_true_transition_probs.R")
source("/z_load_packages.R")

### Define time at which risks will be generated and  calibration will be assessed
t.eval <- round(365.25*7)

run_one_sim <- function(set, sim, scen, n.cohort, n.pctls) {
  
  load(paste0("/small_sample_prep_data_", scen, ".RData"))
  
  set.seed(set + sim) 
  
  # choose cohort
  patids.vec <- sample(1:1000000, n.cohort, replace = FALSE)
  
  # subset data
  data.mstate.reduc <- data.mstate[data.mstate$id %in% patids.vec, ]
  data.raw.reduc    <- data.raw[data.raw$id %in% patids.vec, ]
  
  # output lists for i = 1,2,3 models
  blr_list <- vector("list", 3)
  mlr_list <- vector("list", 3)
  pv_list  <- vector("list", 3)
  
  # loop over the 3 prediction models
  for(i in 1:3){
    tp_i <- tp.pred[[i]]
    tp_i$id <- 1:(nrow(tp_i))
    
    tp.pred.reduc <- tp_i[tp_i$id %in% patids.vec, ] |>
      dplyr::select(!id)
    
    if(scen == "C1"){
      # BLR
      blr_list[[i]] <- calib_msm(data_ms = data.mstate.reduc,
                                 data_raw = data.raw.reduc,
                                 j=1,s=0,t=t.eval,
                                 tp_pred = tp.pred.reduc,
                                 calib_type="blr",
                                 curve_type="rcs",
                                 rcs_nk=4,
                                 assess_moderate=TRUE,
                                 assess_mean=FALSE)[["plotdata"]] 
      
      # MLR
      mlr_list[[i]] <- calib_msm(data_ms = data.mstate.reduc,
                                 data_raw = data.raw.reduc,
                                 j=1,s=0,t=t.eval,
                                 tp_pred = tp.pred.reduc,
                                 calib_type="mlr",
                                 assess_moderate=TRUE,
                                 assess_mean=FALSE)[["plotdata"]] 
      
      # PV
      pv_list[[i]] <- calib_msm(data_ms = data.mstate.reduc,
                                data_raw = data.raw.reduc,
                                j=1,s=0,t=t.eval,
                                tp_pred = tp.pred.reduc,
                                calib_type="pv",
                                curve_type="rcs",
                                rcs_nk=4,
                                assess_moderate=TRUE,
                                assess_mean=FALSE)[["plotdata"]] 
      
    }
    
    else if (scen %in% c("C2", "C3")){
      
      # BLR
      blr_list[[i]] <- calib_msm(data_ms = data.mstate.reduc,
                                 data_raw = data.raw.reduc,
                                 j=1,s=0,t=t.eval,
                                 tp_pred = tp.pred.reduc,
                                 calib_type="blr",
                                 curve_type="rcs",
                                 rcs_nk=4,
                                 w_covs=c("x12", "x13", "x15", "x24", "x25", "x34", "x35", "x45"),
                                 assess_moderate=TRUE,
                                 assess_mean=FALSE)[["plotdata"]] 
      
      # MLR
      mlr_list[[i]] <- calib_msm(data_ms = data.mstate.reduc,
                                 data_raw = data.raw.reduc,
                                 j=1,s=0,t=t.eval,
                                 tp_pred = tp.pred.reduc,
                                 calib_type="mlr",
                                 w_covs=c("x12", "x13", "x15", "x24", "x25", "x34", "x35", "x45"),
                                 assess_moderate=TRUE,
                                 assess_mean=FALSE)[["plotdata"]] 
      
      # PV
      pv_list[[i]] <- calib_msm(data_ms = data.mstate.reduc,
                                data_raw = data.raw.reduc,
                                j=1,s=0,t=t.eval,
                                tp_pred = tp.pred.reduc,
                                calib_type="pv",
                                curve_type="rcs",
                                rcs_nk=4,
                                pv_n_pctls = n.pctls,
                                assess_moderate=TRUE,
                                assess_mean=FALSE)[["plotdata"]] 
    }
  }
  
  return(list(
    blr = blr_list,
    mlr = mlr_list,
    pv  = pv_list
  ))
}

#settings
n.sim    <- 200        
n.cohort <- 700  #500, 400
n.pctls  <- 10
scen     <- "C1" #C2,C3
set      <- 1/50        

cl = parallel::makeCluster(10)
doParallel::registerDoParallel(cl)

clusterEvalQ(cl, {
  
  # Source all functions ONCE per worker
  source("/z_functions.R")
  source("/z_functions_true_transition_probs.R")
  source("/z_load_packages.R")
  
  NULL
})
### Export specific objects used inside meanCal
clusterExport(cl, c("run_one_sim","t.eval", "n.cohort", "n.pctls", "scen","set"))

results_set1 <- foreach(sim = 1:n.sim) %dopar% {
  run_one_sim(set, sim, scen, n.cohort, n.pctls)
}

calib.blr.mod.list <- vector("list", n.sim)
calib.mlr.mod.list <- vector("list", n.sim)
calib.pv.mod.list  <- vector("list", n.sim)

for (sim in 1:n.sim) {
  calib.blr.mod.list[[sim]] <- results_set1[[sim]]$blr
  calib.mlr.mod.list[[sim]] <- results_set1[[sim]]$mlr
  calib.pv.mod.list[[sim]]  <- results_set1[[sim]]$pv
} 

save(
  list = c("calib.blr.mod.list", "calib.mlr.mod.list", "calib.pv.mod.list",
           "scen", "n.cohort", "n.sim", "n.pctls", "set"),
  file = paste0("/small_sample_moderate_analysis_", scen, "_n", n.cohort, "_npctls", n.pctls, "_set", set, ".RData", sep = "")
)
