###
### This program sets the parameters of the simulation and generates
### 1000 datasets containing 1000 sujbects each
###
library(foreach)
library(doParallel)

source("/z_functions.R")
source("/z_functions_true_transition_probs.R")
source("/z_load_packages.R")

### Define key data for simulation 

### Set desired survival probabilities at 7 years
#csh.survival.probs.desired <- c("12" = 0.9, "13" = 0.8, "15" = 0.99, "24" = 0.55, "25" = 0.95, "34" = 0.7, "35" = 0.15, "45" = 0.05)
csh.survival.probs.desired <- c("12" = 0.4, "13" = 0.25, "15" = 0.99, "24" = 0.55, "25" = 0.95, "34" = 0.5, "35" = 0.15, "45" = 0.2)

### Apply function to calculate the desired scales
scales.sim <- sapply(csh.survival.probs.desired, calc.scale)
scales.sim

### Define time at which risks will be generated and  calibration will be assessed
t.eval <- round(365.25*7)

### Defne the max follow up to be after t.eval
max.follow <- t.eval + 2

### Choose cohort size
n.cohort <- 1000

#function the produces the dataset
sim_data <- function(seed, n.cohort, scales.sim, t.eval, max.follow) {
  
  print(paste("seed =", seed))
  print("START DATA GEN")
  
  ## Generate baseline data
  set.seed(seed)
  x.baseline <- data.frame("x12" = rnorm(n.cohort, 0, 1),
                           "x13" = rnorm(n.cohort, 0, 1),
                           "x15" = rnorm(n.cohort, 0, 1),
                           "x24" = rnorm(n.cohort, 0, 1),
                           "x25" = rnorm(n.cohort, 0, 1),
                           "x34" = rnorm(n.cohort, 0, 1),
                           "x35" = rnorm(n.cohort, 0, 1),
                           "x45" = rnorm(n.cohort, 0, 1))
  
  Sys.time()
  
  ## Generate data
  data.out <- gen.dat.DGM1(n = n.cohort, #number of patients to simulate
                           max.follow = max.follow, #max follow up, informative censoring 1 day after 7 years
                           shape12 = 1, scale12 = scales.sim["12"], #shape and scale for weibull baseline hazard for transition 1 -> 2
                           shape13 = 1, scale13 = scales.sim["13"], #shape and scale for weibull baseline hazard for transition 1 -> 3
                           shape15 = 1, scale15 = scales.sim["15"], #shape and scale for weibull baseline hazard for transition 1 -> 5
                           shape24 = 1, scale24 = scales.sim["24"], #shape and scale for weibull baseline hazard for transition 2 -> 4
                           shape25 = 1, scale25 = scales.sim["25"], #shape and scale for weibull baseline hazard for transition 2 -> 5
                           shape34 = 1, scale34 = scales.sim["34"], #shape and scale for weibull baseline hazard for transition 3 -> 4
                           shape35 = 1, scale35 = scales.sim["35"], #shape and scale for weibull baseline hazard for transition 3 -> 5
                           shape45 = 1, scale45 = scales.sim["45"], #shape and scale for weibull baseline hazard for transition 4 -> 5
                           beta.x12 = 0.5, #covariate effects for transiion 12
                           beta.x13 = 0.5, #covariate effects for transiion 13
                           beta.x15 = 0.5, #covariate effects for transiion 15
                           beta.x24 = 0.5, #covariate effects for transiion 24
                           beta.x25 = 0.5, #covariate effects for transiion 25
                           beta.x34 = 0.5, #covariate effects for transiion 34
                           beta.x35 = 0.5, #covariate effects for transiion 35
                           beta.x45 = 0.5, #covariate effects for transiion 45
                           x.in = x.baseline, #baseline predictors, dataframe with two columns (x1 continuous, x2 binary)
                           numsteps = max.follow)
  
  str(data.out)
  
  ###
  ### Calculate true transition probabilities for each individual
  ###
  
  ### Write a function to it for a single individual
  print(paste("CALC TRUE RISKS ", Sys.time(), sep = ""))
  calc.true.tp.ind <- function(row){
    return(calc.true.transition.probs.DGM1(0, t.eval, row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8],
                                           shape12 = 1, scale12 = scales.sim["12"], #shape and scale for weibull baseline hazard for transition 1 -> 2
                                           shape13 = 1, scale13 = scales.sim["13"], #shape and scale for weibull baseline hazard for transition 1 -> 3
                                           shape15 = 1, scale15 = scales.sim["15"], #shape and scale for weibull baseline hazard for transition 1 -> 5
                                           shape24 = 1, scale24 = scales.sim["24"], #shape and scale for weibull baseline hazard for transition 2 -> 4
                                           shape25 = 1, scale25 = scales.sim["25"], #shape and scale for weibull baseline hazard for transition 2 -> 5
                                           shape34 = 1, scale34 = scales.sim["34"], #shape and scale for weibull baseline hazard for transition 3 -> 4
                                           shape35 = 1, scale35 = scales.sim["35"], #shape and scale for weibull baseline hazard for transition 3 -> 5
                                           shape45 = 1, scale45 = scales.sim["45"], #shape and scale for weibull baseline hazard for transition 3 -> 5
                                           beta.x12 = 0.5, #covariate effects for transiion 12
                                           beta.x13 = 0.5, #covariate effects for transiion 13
                                           beta.x15 = 0.5, #covariate effects for transiion 15
                                           beta.x24 = 0.5, #covariate effects for transiion 24
                                           beta.x25 = 0.5, #covariate effects for transiion 25
                                           beta.x34 = 0.5, #covariate effects for transiion 34
                                           beta.x35 = 0.5, #covariate effects for transiion 35
                                           beta.x45 = 0.5 #covariate effects for transiion 45
    ))
  }
  Sys.time()
  
  ### Calc true risk for each individual in this dataset
  p.true <- apply(dplyr::select(data.out[["cohort"]], x12, x13, x15, x24, x25, x34, x35, x45), 1, calc.true.tp.ind)
  p.true <- data.frame(t(p.true))
  colnames(p.true) <- paste("p.true", 1:5, sep = "")
  print(paste("FINISH TRUE RISKS ", Sys.time(), sep = ""))
  
  ### Combine data into one data frame
  data.out <- cbind(data.out[[1]], p.true)
  ### Save data object
  saveRDS(data.out, paste("/File1/dataout", seed, ".rds", sep = ""))
}

#parallelize the iterations
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
clusterExport(cl, c("csh.survival.probs.desired","t.eval", "n.cohort", "sim_data", "max.follow","scales.sim"))

results_set1 <- foreach(seed = 1:1000) %dopar% {
  
  sim_data(seed, n.cohort, scales.sim, t.eval, max.follow)
  NULL
}

