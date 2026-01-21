###
### Estimating the true calibration curves
###
source("C:/Users/userr/Desktop/ms_thesis/simulation/original_code/z_functions.R")
source("C:/Users/userr/Desktop/ms_thesis/simulation/original_code/z_functions_true_transition_probs.R")
source("C:/Users/userr/Desktop/ms_thesis/simulation/original_code/z_load_packages.R")

### True calibration is the same irrespective of the censoring that has been applied, 
### so just need to load data for one scenario
load(paste0("/small_sample_prep_data_C1.RData"))

# choose cohort
patids.vec <- sample(1:1000000, 5000, replace = FALSE)

# subset data
data.mstate.reduc <- data.mstate[data.mstate$id %in% patids.vec, ]
data.raw.reduc    <- data.raw[data.raw$id %in% patids.vec, ]

### This will be done by regressing the predicted risks on the true risks, and generating predicted-observed risks using this model
### We can treat the true risks as pseudo-values and insert into calib_msm using pv.precalc to achieve this

true.calib <- vector("list", 3) 

rcs.nk= 4

for(i in 1:3){
  tp_i <- tp.pred[[i]]
  tp_i$id <- 1:(nrow(tp_i))
  
  tp.pred.reduc <- tp_i[tp_i$id %in% patids.vec, ] |>
    dplyr::select(!id)
  if(i==1){tp_corr=tp.pred.reduc}
  
  ### Estimate true calibration plot
  true.calib[[i]] <- calib_msm(data_ms = data.mstate.reduc, 
                               data_raw = data.raw.reduc,
                               j = 1,
                               s = 0,
                               t = t.eval,
                               calib_type = "pv",
                               tp_pred = tp.pred.reduc,
                               pv_precalc = tp_corr,
                               curve_type = "rcs",
                               rcs_nk = rcs.nk)[["plotdata"]]
}

saveRDS(true.calib, paste("/true.calib.lineplot.rcs.nk", rcs.nk, ".rds", sep = ""))

