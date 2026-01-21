###
### This program runs the assessment of mean calibration in the small sample simulation
###

library(foreach)
library(doParallel)

source("/z_functions.R")
source("/z_functions_true_transition_probs.R")
source("/z_load_packages.R")

### Define time at which risks will be generated and  calibration will be assessed
t.eval <- round(365.25*7)

meanCal <- function(scen, n.cohort, n.pctls, set, n.sim){
  ### Set seed
  set.seed(set)
  
  ### Create list to store calibration of each estimates
  ### (include a list to store true calibration in each simulation iteration)
  calib.blr.mean.list <- vector("list", 3)
  calib.mlr.mean.list <- vector("list", 3)
  calib.aj.mean.list <- vector("list", 3)
  calib.true.mean.list <- vector("list", 3)
  
  ### Calculate the calibration intecepts and slopes
  for (i in 1:3){
    
    ### BLR-IPCW
    calib.blr.mean.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.blr.mean.list[[i]]) <- paste("state", 1:5, sep = "")
    
    ### MLR-IPCW
    calib.mlr.mean.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.mlr.mean.list[[i]]) <- paste("state", 1:5, sep = "")
    
    ### Pseudo-value
    calib.aj.mean.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.aj.mean.list[[i]]) <- paste("state", 1:5, sep = "")
    
    ### True
    calib.true.mean.list[[i]] <- data.frame(matrix(NA, nrow = n.sim, ncol = 5))
    colnames(calib.true.mean.list[[i]]) <- paste("state", 1:5, sep = "")
    
  }
  
  ###
  ### Run the simulation loop
  ###
  for (sim in 1:n.sim){
    
    print(paste("cohort size = ", n.cohort, ", sim =", sim, Sys.time(), sep = " "))
    
    ### Choose patids at random
    patids.vec <- sample((1:1000000), n.cohort, replace = FALSE)
    # We exclude individual 390375, as this person was removed from M1C2 dataset as they caused an error
    
    ### Extract individuals
    data.mstate.reduc <- data.mstate[data.mstate$id %in% patids.vec, ]
    data.raw.reduc <- data.raw[data.raw$id %in% patids.vec, ]
    
    ###
    ### Calculate calibration for each type of predicted risk
    ###
    for (i in 1:3){
      
      ### Extract predicted risks
      tp.pred[[i]]$id <- 1:nrow(tp.pred[[i]])
      tp.pred.reduc <- tp.pred[[i]][tp.pred[[i]]$id %in% patids.vec, ] |>
        dplyr::select(!id)
      
      if (scen == "C1"){
        
        ###
        ### Estimate mean calibration with BLR-IPCW
        calib.mean <- calib_msm(data_ms = data.mstate.reduc,
                                data_raw = data.raw.reduc,
                                j = 1,
                                s = 0,
                                t = t.eval,
                                tp_pred = tp.pred.reduc,
                                calib_type = "blr",
                                assess_moderate = FALSE,
                                assess_mean = TRUE)[["mean"]] 
        
        ### Assign output
        calib.blr.mean.list[[i]][sim, ] <- calib.mean
        
        ###
        ### Estimate mean calibration with MLR-IPCW
        calib.mean <- calib_msm(data_ms = data.mstate.reduc,
                                data_raw = data.raw.reduc,
                                j = 1,
                                s = 0,
                                t = t.eval,
                                tp_pred = tp.pred.reduc,
                                calib_type = "mlr",
                                assess_moderate = FALSE,
                                assess_mean = TRUE)[["mean"]] 
        
        ### Assign output
        calib.mlr.mean.list[[i]][sim, ] <- calib.mean
        
        ###
        ### Estimate mean calibration with Aalen-Johansen
        calib.mean <- calib_msm(data_ms = data.mstate.reduc,
                                data_raw = data.raw.reduc,
                                j = 1,
                                s = 0,
                                t = t.eval,
                                tp_pred = tp.pred.reduc,
                                calib_type = "aj",
                                assess_moderate = FALSE,
                                assess_mean = TRUE)[["mean"]]
        
        ### Assign output
        calib.aj.mean.list[[i]][sim, ] <- calib.mean
        
      } else if (scen %in% c("C2", "C3")){
        
        ###
        ### Estimate mean calibration with BLR-IPCW
        calib.mean <- calib_msm(data_ms = data.mstate.reduc,
                                data_raw = data.raw.reduc,
                                j = 1,
                                s = 0,
                                t = t.eval,
                                tp_pred = tp.pred.reduc,
                                calib_type = "blr",
                                w_covs = c("x12", "x13", "x15", "x24", "x25", "x34", "x35", "x45"),
                                assess_moderate = FALSE,
                                assess_mean = TRUE)[["mean"]] 
        
        ### Assign output
        calib.blr.mean.list[[i]][sim, ] <- calib.mean
        
        ###
        ### Estimate mean calibration with MLR-IPCW
        calib.mean <- calib_msm(data_ms = data.mstate.reduc,
                                data_raw = data.raw.reduc,
                                j = 1,
                                s = 0,
                                t = t.eval,
                                tp_pred = tp.pred.reduc,
                                calib_type = "mlr",
                                w_covs = c("x12", "x13", "x15", "x24", "x25", "x34", "x35", "x45"),
                                assess_moderate = FALSE,
                                assess_mean = TRUE)[["mean"]] 
        
        ### Assign output
        calib.mlr.mean.list[[i]][sim, ] <- calib.mean
        
        ###
        ### Estimate mean calibration with Aalen-Johansen
        calib.mean <- calib_msm(data_ms = data.mstate.reduc,
                                data_raw = data.raw.reduc,
                                j = 1,
                                s = 0,
                                t = t.eval,
                                tp_pred = tp.pred.reduc,
                                calib_type = "aj",
                                pv_n_pctls = n.pctls,
                                assess_moderate = FALSE,
                                assess_mean = TRUE)[["mean"]]
        
        ### Assign output
        calib.aj.mean.list[[i]][sim, ] <- calib.mean
        
      }
      
      ###
      ### Save true calibration
      calib.true.mean.list[[i]][sim, ] <- colMeans(data.raw.reduc[, paste("p.true", 1:5, sep = "")] - tp.pred.reduc)
      
    }
    
  }
  save(
    list = c("calib.blr.mean.list", "calib.mlr.mean.list", "calib.aj.mean.list",
             "calib.true.mean.list", "scen", "n.cohort", "n.sim", "n.pctls", "set"),
    file = paste0("/small_sample_mean_analysis_",
                  scen, "_n", n.cohort, "_npctls", n.pctls, "_set", set, ".RData")
  )
  
}

# settings
nRepl <- 50
n.cohort <- 700 #500/400
n.pctls  <- 10
n.sim    <- 1000
scen <- "C1" #C2/C3

cl = parallel::makeCluster(6)
doParallel::registerDoParallel(cl)

clusterEvalQ(cl, {
  
  # Source all functions ONCE per worker
  source("/z_functions.R")
  source("/z_functions_true_transition_probs.R")
  source("/z_load_packages.R")
  
  # Load big dataset ONCE per worker
  scen_local <- "C1" #C2/C3
  load(paste0("/small_sample_prep_data_", scen_local, ".RData"))
  
  NULL
})
### Export specific objects used inside meanCal
clusterExport(cl, c("meanCal", "t.eval", "n.cohort", "n.pctls", "n.sim", "scen"))


# 2: switch do %dopar% and export necessary packages
results = foreach(i = 1:nRepl) %dopar% {
  meanCal(scen, n.cohort, n.pctls, i, n.sim)
}

stopCluster(cl)
