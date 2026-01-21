###
### This program will prep the data ready to calculate calibration for the cohort size and scenario of interest
### Aalen-Johansen estimator of the entire cohort will be calculated
###

#-------------------------------C1----------------------------

scen <- "C1"
print(paste("scen = ", scen, sep = ""))

### Load simulated data
data.pop.mstate <- readRDS(paste("/censC1/pop.mstate.dgm1.", scen, ".RData", sep = ""))
data.pop.raw <- readRDS(paste("/censC1/pop.raw.dgm1.", scen, ".RData", sep = ""))
tmat <- readRDS(paste("/censC1/tmat.dgm1.RData", sep = ""))

### Extract all individuals, which is what we're doing the large sample analysis with
data.mstate <- data.pop.mstate
data.raw <- data.pop.raw
str(data.mstate)
str(data.raw)

### Extract true transition probabilities
tp.true <- data.raw %>% select(paste("p.true", 1:5, sep = ""))

### Now defined the predicted transition probabilities based off these
tp.pred <- vector("list", 3)

### First element is perfectly prediction probabilities
tp.pred[[1]] <- dplyr::mutate(tp.true,
                              tp1 = p.true1,
                              tp2 = p.true2,
                              tp3 = p.true3,
                              tp4 = p.true4,
                              tp5 = p.true5) |>
  dplyr::select(tp1, tp2, tp3, tp4, tp5)

### Second set of predicted transition probabilities we add vector (0.5, 0.25, 0, -0.25, -0.5) to the log-odds, convert back and normalise
tp.pred[[2]] <- dplyr::mutate(tp.true,
                              tp1 = exp(log(p.true1/(1 - p.true1)) + 0.5)/(1 + exp(log(p.true1/(1 - p.true1)) + 0.5)),
                              tp2 = exp(log(p.true2/(1 - p.true2)) + 0.25)/(1 + exp(log(p.true2/(1 - p.true2)) + 0.25)),
                              tp3 = p.true3,
                              tp4 = exp(log(p.true4/(1 - p.true4)) - 0.25)/(1 + exp(log(p.true4/(1 - p.true4)) - 0.25)),
                              tp5 = exp(log(p.true5/(1 - p.true5)) - 0.5)/(1 + exp(log(p.true5/(1 - p.true5)) - 0.5)),
                              tp1 = tp1/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp2 = tp2/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp3 = tp3/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp4 = tp4/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp5 = tp5/(tp1 + tp2 + tp3 + tp4 + tp5)) |>
  dplyr::select(tp1, tp2, tp3, tp4, tp5)

### Third set of predicted transition probabilities  (0.5, -0.5, -0.5, -0.5, 0.5) to the log-odds, convert back and normalise
tp.pred[[3]] <- dplyr::mutate(tp.true,
                              tp1 = exp(log(p.true1/(1 - p.true1)) + 0.5)/(1 + exp(log(p.true1/(1 - p.true1)) + 0.5)),
                              tp2 = exp(log(p.true2/(1 - p.true2)) - 0.5)/(1 + exp(log(p.true2/(1 - p.true2)) - 0.5)),
                              tp3 = exp(log(p.true3/(1 - p.true3)) - 0.5)/(1 + exp(log(p.true3/(1 - p.true3)) - 0.5)),
                              tp4 = exp(log(p.true4/(1 - p.true4)) - 0.5)/(1 + exp(log(p.true4/(1 - p.true4)) - 0.5)),
                              tp5 = exp(log(p.true5/(1 - p.true5)) + 0.5)/(1 + exp(log(p.true5/(1 - p.true5)) + 0.5)),
                              tp1 = tp1/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp2 = tp2/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp3 = tp3/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp4 = tp4/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp5 = tp5/(tp1 + tp2 + tp3 + tp4 + tp5)) |>
  dplyr::select(tp1, tp2, tp3, tp4, tp5)

### Need to create the censoring variable which will be used to estimate th weights when assessing calibration
### dtcens is time until either censored, or entering an absorbing state
### dtcens.s == 1 if individual is censored, dtcens.s == 0 if absorbing state is entered
### Also add an "id" variable which is required for calibmsm functionality
data.raw <- mutate(data.raw,
                   dtcens = pmin(State.5, cens.times, na.rm = TRUE),
                   dtcens_s = case_when(dtcens == cens.times ~ 1,
                                        dtcens == State.5 ~ 0),
                   id = patid)

data.mstate$id <- data.mstate$patid

###
### Save image
rm(data.pop.raw, data.pop.mstate)
save.image(paste("/small_sim_C1/small_sample_prep_data_", scen, ".RData", sep = ""))
print("IMAGE SAVED")

#----------------------------C2-------------------------------------

scen <- "C2"
print(paste("scen = ", scen, sep = ""))

### Load simulated data
data.pop.mstate <- readRDS(paste("/censC2/pop.mstate.dgm1.", scen, ".RData", sep = ""))
data.pop.raw <- readRDS(paste("/censC2/pop.raw.dgm1.", scen, ".RData", sep = ""))
tmat <- readRDS(paste("/censC2/tmat.dgm1.RData", sep = ""))

### Extract all individuals, which is what we're doing the large sample analysis with
data.mstate <- data.pop.mstate
data.raw <- data.pop.raw
str(data.mstate)
str(data.raw)

### Extract true transition probabilities
tp.true <- data.raw %>% select(paste("p.true", 1:5, sep = ""))

### Now defined the predicted transition probabilities based off these
tp.pred <- vector("list", 3)

### First element is perfectly prediction probabilities
tp.pred[[1]] <- dplyr::mutate(tp.true,
                              tp1 = p.true1,
                              tp2 = p.true2,
                              tp3 = p.true3,
                              tp4 = p.true4,
                              tp5 = p.true5) |>
  dplyr::select(tp1, tp2, tp3, tp4, tp5)

### Second set of predicted transition probabilities we add vector (0.5, 0.25, 0, -0.25, -0.5) to the log-odds, convert back and normalise
tp.pred[[2]] <- dplyr::mutate(tp.true,
                              tp1 = exp(log(p.true1/(1 - p.true1)) + 0.5)/(1 + exp(log(p.true1/(1 - p.true1)) + 0.5)),
                              tp2 = exp(log(p.true2/(1 - p.true2)) + 0.25)/(1 + exp(log(p.true2/(1 - p.true2)) + 0.25)),
                              tp3 = p.true3,
                              tp4 = exp(log(p.true4/(1 - p.true4)) - 0.25)/(1 + exp(log(p.true4/(1 - p.true4)) - 0.25)),
                              tp5 = exp(log(p.true5/(1 - p.true5)) - 0.5)/(1 + exp(log(p.true5/(1 - p.true5)) - 0.5)),
                              tp1 = tp1/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp2 = tp2/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp3 = tp3/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp4 = tp4/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp5 = tp5/(tp1 + tp2 + tp3 + tp4 + tp5)) |>
  dplyr::select(tp1, tp2, tp3, tp4, tp5)

### Third set of predicted transition probabilities  (0.5, -0.5, -0.5, -0.5, 0.5) to the log-odds, convert back and normalise
tp.pred[[3]] <- dplyr::mutate(tp.true,
                              tp1 = exp(log(p.true1/(1 - p.true1)) + 0.5)/(1 + exp(log(p.true1/(1 - p.true1)) + 0.5)),
                              tp2 = exp(log(p.true2/(1 - p.true2)) - 0.5)/(1 + exp(log(p.true2/(1 - p.true2)) - 0.5)),
                              tp3 = exp(log(p.true3/(1 - p.true3)) - 0.5)/(1 + exp(log(p.true3/(1 - p.true3)) - 0.5)),
                              tp4 = exp(log(p.true4/(1 - p.true4)) - 0.5)/(1 + exp(log(p.true4/(1 - p.true4)) - 0.5)),
                              tp5 = exp(log(p.true5/(1 - p.true5)) + 0.5)/(1 + exp(log(p.true5/(1 - p.true5)) + 0.5)),
                              tp1 = tp1/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp2 = tp2/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp3 = tp3/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp4 = tp4/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp5 = tp5/(tp1 + tp2 + tp3 + tp4 + tp5)) |>
  dplyr::select(tp1, tp2, tp3, tp4, tp5)

### Need to create the censoring variable which will be used to estimate th weights when assessing calibration
### dtcens is time until either censored, or entering an absorbing state
### dtcens.s == 1 if individual is censored, dtcens.s == 0 if absorbing state is entered
### Also add an "id" variable which is required for calibmsm functionality
data.raw <- mutate(data.raw,
                   dtcens = pmin(State.5, cens.times, na.rm = TRUE),
                   dtcens_s = case_when(dtcens == cens.times ~ 1,
                                        dtcens == State.5 ~ 0),
                   id = patid)

data.mstate$id <- data.mstate$patid

###
### Save image
rm(data.pop.raw, data.pop.mstate)
save.image(paste("/small_sim_C2/small_sample_prep_data_", scen, ".RData", sep = ""))
print("IMAGE SAVED")




#-----------------------------C3---------------------------------------

scen <- "C3"
print(paste("scen = ", scen, sep = ""))

### Load simulated data
data.pop.mstate <- readRDS(paste("/censC3/pop.mstate.dgm1.", scen, ".RData", sep = ""))
data.pop.raw <- readRDS(paste("/censC3/pop.raw.dgm1.", scen, ".RData", sep = ""))
tmat <- readRDS(paste("/censC3/tmat.dgm1.RData", sep = ""))

### Extract all individuals, which is what we're doing the large sample analysis with
data.mstate <- data.pop.mstate
data.raw <- data.pop.raw
str(data.mstate)
str(data.raw)

### Extract true transition probabilities
tp.true <- data.raw %>% select(paste("p.true", 1:5, sep = ""))

### Now defined the predicted transition probabilities based off these
tp.pred <- vector("list", 3)

### First element is perfectly prediction probabilities
tp.pred[[1]] <- dplyr::mutate(tp.true,
                              tp1 = p.true1,
                              tp2 = p.true2,
                              tp3 = p.true3,
                              tp4 = p.true4,
                              tp5 = p.true5) |>
  dplyr::select(tp1, tp2, tp3, tp4, tp5)

### Second set of predicted transition probabilities we add vector (0.5, 0.25, 0, -0.25, -0.5) to the log-odds, convert back and normalise
tp.pred[[2]] <- dplyr::mutate(tp.true,
                              tp1 = exp(log(p.true1/(1 - p.true1)) + 0.5)/(1 + exp(log(p.true1/(1 - p.true1)) + 0.5)),
                              tp2 = exp(log(p.true2/(1 - p.true2)) + 0.25)/(1 + exp(log(p.true2/(1 - p.true2)) + 0.25)),
                              tp3 = p.true3,
                              tp4 = exp(log(p.true4/(1 - p.true4)) - 0.25)/(1 + exp(log(p.true4/(1 - p.true4)) - 0.25)),
                              tp5 = exp(log(p.true5/(1 - p.true5)) - 0.5)/(1 + exp(log(p.true5/(1 - p.true5)) - 0.5)),
                              tp1 = tp1/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp2 = tp2/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp3 = tp3/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp4 = tp4/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp5 = tp5/(tp1 + tp2 + tp3 + tp4 + tp5)) |>
  dplyr::select(tp1, tp2, tp3, tp4, tp5)

### Third set of predicted transition probabilities  (0.5, -0.5, -0.5, -0.5, 0.5) to the log-odds, convert back and normalise
tp.pred[[3]] <- dplyr::mutate(tp.true,
                              tp1 = exp(log(p.true1/(1 - p.true1)) + 0.5)/(1 + exp(log(p.true1/(1 - p.true1)) + 0.5)),
                              tp2 = exp(log(p.true2/(1 - p.true2)) - 0.5)/(1 + exp(log(p.true2/(1 - p.true2)) - 0.5)),
                              tp3 = exp(log(p.true3/(1 - p.true3)) - 0.5)/(1 + exp(log(p.true3/(1 - p.true3)) - 0.5)),
                              tp4 = exp(log(p.true4/(1 - p.true4)) - 0.5)/(1 + exp(log(p.true4/(1 - p.true4)) - 0.5)),
                              tp5 = exp(log(p.true5/(1 - p.true5)) + 0.5)/(1 + exp(log(p.true5/(1 - p.true5)) + 0.5)),
                              tp1 = tp1/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp2 = tp2/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp3 = tp3/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp4 = tp4/(tp1 + tp2 + tp3 + tp4 + tp5),
                              tp5 = tp5/(tp1 + tp2 + tp3 + tp4 + tp5)) |>
  dplyr::select(tp1, tp2, tp3, tp4, tp5)

### Need to create the censoring variable which will be used to estimate th weights when assessing calibration
### dtcens is time until either censored, or entering an absorbing state
### dtcens.s == 1 if individual is censored, dtcens.s == 0 if absorbing state is entered
### Also add an "id" variable which is required for calibmsm functionality
data.raw <- mutate(data.raw,
                   dtcens = pmin(State.5, cens.times, na.rm = TRUE),
                   dtcens_s = case_when(dtcens == cens.times ~ 1,
                                        dtcens == State.5 ~ 0),
                   id = patid)

data.mstate$id <- data.mstate$patid

###
### Save image
rm(data.pop.raw, data.pop.mstate)
save.image(paste("/small_sim_C3/small_sample_prep_data_", scen, ".RData", sep = ""))
print("IMAGE SAVED")

