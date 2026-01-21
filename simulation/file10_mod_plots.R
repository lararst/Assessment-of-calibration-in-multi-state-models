source("/z_functions.R")
source("/z_functions_true_transition_probs.R")
source("/z_load_packages.R")
###
### This program will create plots for the small sample analysis, moderate calibration
###

### Load the true calibration data
true.calib <- readRDS(paste("/true.calib.lineplot.rcs.nk", 4, ".rds", sep = ""))

### This function will create the plots
create_plots_moderate_small_sample <- function(scen, scen.tp, n.cohort, n.pctls, n.plots){
  
  ### Load the calibration plot data
  load(paste("/small_sim_", scen, "/mod_", scen, "_", n.cohort, "/small_sample_moderate_analysis_", scen, "_n", n.cohort, "_npctls", n.pctls, ".RData", sep = ""))
  
  ### Reduce true.calib to scenario of interest, and 5,000 observations
  true.calib <- lapply( 1:5, function(x) true.calib[[scen.tp]][[x]])
  
  ### Reduce to scep.tp data
  calib.pv.mod.list.comb <- calib.pv.mod.list.comb[[scen.tp]]
  calib.blr.mod.list.comb <- calib.blr.mod.list.comb[[scen.tp]]
  calib.mlr.mod.list.comb <- calib.mlr.mod.list.comb[[scen.tp]]
  
  ### Reduce each plot to 500 observations (robust to missing plots)
  calib.pv.mod.list.comb <- lapply(seq_along(calib.pv.mod.list.comb), function(x) {
    lapply(calib.pv.mod.list.comb[[x]], function(y) {
      out <- y[1:min(500, nrow(y)), ]
      out$plot <- x
      out
    })
  })
  
  calib.blr.mod.list.comb <- lapply(seq_along(calib.blr.mod.list.comb), function(x) {
    lapply(calib.blr.mod.list.comb[[x]], function(y) {
      out <- y[1:min(500, nrow(y)), ]
      out$plot <- x
      out
    })
  })
  
  calib.mlr.mod.list.comb <- lapply(seq_along(calib.mlr.mod.list.comb), function(x) {
    lapply(calib.mlr.mod.list.comb[[x]], function(y) {
      out <- y[1:min(500, nrow(y)), ]
      out$plot <- x
      out
    })
  })
  
  ### Now want to concatenate into a single data for each state
  calib.pv.mod.list.comb <- lapply(1:5, function(x) {do.call("rbind", lapply(calib.pv.mod.list.comb, function(y) {y[[x]]}))})
  calib.blr.mod.list.comb <- lapply(1:5, function(x) {do.call("rbind", lapply(calib.blr.mod.list.comb, function(y) {y[[x]]}))})
  calib.mlr.mod.list.comb <- lapply(1:5, function(x) {do.call("rbind", lapply(calib.mlr.mod.list.comb, function(y) {y[[x]]}))})
  
  common.cols <- c("id", "pred", "obs", "plot")
  
  
  calib.pv.mod.list.comb <- lapply(1:5, function(x) {
    
    sim.df <- calib.pv.mod.list.comb[[x]][, common.cols]
    
    true.df <- data.frame(
      true.calib[[x]][, c("id", "pred", "obs")],
      plot = 0
    )
    
    rbind(sim.df, true.df)
  })
  
  calib.blr.mod.list.comb <- lapply(1:5, function(x) {
    
    sim.df <- calib.blr.mod.list.comb[[x]][, common.cols]
    
    true.df <- data.frame(
      true.calib[[x]][, c("id", "pred", "obs")],
      plot = 0
    )
    
    rbind(sim.df, true.df)
  })
  
  calib.mlr.mod.list.comb <- lapply(1:5, function(x) {
    
    sim.df <- calib.mlr.mod.list.comb[[x]][, common.cols]
    
    true.df <- data.frame(
      true.calib[[x]][, c("id", "pred", "obs")],
      plot = 0
    )
    
    rbind(sim.df, true.df)
  })
  
  ### Create calibration plots
  
  ### Pseudo-value
  calib.plots.pv <- 
    plot_calib_msm_small_sample(calib.pv.mod.list.comb, 
                                calib.type = "line",
                                marg.density = TRUE, 
                                marg.density.size = 8, 
                                ncol = 5, 
                                nrow = 1, 
                                inclu.legend = TRUE,
                                legend.seperate = TRUE,
                                legend.title = "Line plots: ",
                                inclu.title = TRUE,
                                axis.titles.x = 0,
                                axis.titles.y = 1, 
                                axis.titles.text.y = "Pseudo-value\nObserved risk",
                                size = 16)
  
  ### BLR-IPCW
  calib.plots.blr <- 
    plot_calib_msm_small_sample(calib.blr.mod.list.comb, 
                                calib.type = "line",
                                marg.density = TRUE, 
                                marg.density.size = 8, 
                                ncol = 5, 
                                nrow = 1, 
                                inclu.legend = FALSE,
                                inclu.title = FALSE,
                                axis.titles.text.x = "Predicted risk",
                                axis.titles.y = 1, 
                                axis.titles.text.y = "BLR-IPCW\nObserved risk",
                                size = 16)
  
  #   ### MLR-IPCW
  calib.plots.mlr <- 
    plot_calib_msm_small_sample(calib.mlr.mod.list.comb, 
                                calib.type = "scatter",
                                marg.density = TRUE, 
                                marg.density.size = 8, 
                                ncol = 5, 
                                nrow = 1, 
                                inclu.legend = TRUE,
                                legend.seperate = TRUE,
                                legend.title = "Scatter plots: ",
                                inclu.title = FALSE,
                                axis.titles.x = 0,
                                axis.titles.y = 1, 
                                axis.titles.text.y = "MLR-IPCW\nObserved risk",
                                size = 16)
  
  ### Combine into a single plot for the Figure
  plots.fig <- gridExtra::arrangeGrob(calib.plots.pv[["plots"]], calib.plots.blr, calib.plots.mlr[["plots"]], 
                                      gridExtra::arrangeGrob(calib.plots.pv[["legend"]], calib.plots.mlr[["legend"]], ncol = 2),
                                      nrow = 4, heights = c(1, 1, 1, 0.05))
  
  
  CairoPNG(paste("/small_sample_moderate_", scen, "_tp", scen.tp, "_N", n.cohort, ".png", sep = ""), width = 20, height = 8, unit = "in", res = 150)
  grid::grid.draw(plots.fig)
  dev.off()
  CairoPDF(paste("/small_sample_moderate_", scen,  "_tp", scen.tp,  "_N", n.cohort, ".pdf", sep = ""), width = 20, height = 8)
  grid::grid.draw(plots.fig)
  dev.off()
  
}

lapply(c("C1","C2","C3"), function(x) {lapply(1:3, function(y) {create_plots_moderate_small_sample(x, y, 500, 10, 200)})})
