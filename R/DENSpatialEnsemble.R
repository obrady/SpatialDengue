#' Run a multiple realisations of the spatial dengue model and ensemble them
#'
#' Runs multiple versions of DEN.spatial and calculates confidence intervals
#' @param weekdates Two element vector of the start and end weeks of the simulation over which the mdoel will be evaluated over
#' @param fitdat Data frame of the locations, numbers and timings (in weeks) of cases to fit the model to, see ?sgdat
#' @param pastdat Data frame of the locations, numbers and timings (in weeks) of all cases in the dataset (is used to generate the starting immunity profile), see ?sgdat
#' @param unipix Universal pixel lookup table, see ?make.unipix
#' @param pixdistmat A patch distance matrix, see example
#' @param steprun integer, number of days for which the simulation should run, excluding burn in period
#' @param nruns Integer, number of stochastic model runs over which to average over
#' @param qtiles vector of quantiles over which results should be summarised (one or more values between 0 and 1)
#' @param ncores number of cpu cores to run simulations over (in parallel). Uses snowfall for parallelisation. Specify cores > 1 (default = 1) for parallel processing
#' @param paramsList Optional parameter list. If not supplied returns to defaults, see tutorial for full parameter list, see ?model.run for full list and explanation of parameters
#' @param seasonal Logical, whether to use the defaul DEN.spatial function or the seasonal DEN.spatial.seasonal version. If this equals TRUE please supply a seasonal vector see ?DEN.spatial.seasonal
#' @param seasonal_vector A vector of seasonal transmission intensity. See ?DEN.spatial.seasonal
#' @details Returns a list with each component summarising the results across multiple model runs at a given quantile, e.g. 0.95, 0.5, 0.05, etc.
#' @keywords model simulation ensemble
#' @export
#' @examples
#' data(sgdat)
#' data(sgpop)
#' sgpop <- pop.process(sgpop, agg = 10)
#' unipix <- make.unipix(sgpop)
#' pixdistmat <- distm(cbind(unipix$x, unipix$y))
#' sgdat <- data.frame(sgdat, patchID = apply(cbind(sgdat[, 3:2]), 1, pix.id.find, unipix))
#' weekdates <- c(40, 92)
#' # run model ensemble with default parameters
#' denmod_sim <- DEN.spatial.ensemble(weekdates, sgdat, sgdat, unipix, pixdistmat, 365, 10)

DEN.spatial.ensemble <- function(weekdates,
                                 fitdat,
                                 pastdat,
                                 unipix,
                                 pixdistmat,
                                 steprun,
                                 nruns,
                                 qtiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                 ncores = 1,
                                 paramsList = NULL,
                                 seasonal = FALSE,
                                 seasonal_vector = NA){
  
  # warning message
  if(seasonal & all(is.na(seasonal_vector))){print("please supply a seasonal_vector")}
  
  Smat <- matrix(NA, nrow = steprun, ncol = nruns)
  Imat <- matrix(NA, nrow = steprun, ncol = nruns)
  Rmat <- matrix(NA, nrow = steprun, ncol = nruns)
  Rtmat <- matrix(NA, nrow = steprun, ncol = nruns)
  Dmat <- matrix(NA, nrow = steprun, ncol = nruns)
  newDmat <- matrix(NA, nrow = steprun, ncol = nruns)
  mosEmat <- matrix(NA, nrow = steprun, ncol = nruns)
  mosImat <- matrix(NA, nrow = steprun, ncol = nruns)
  enmodlist <- list()
  # running the models
  # run in parallel if ncores > 1
  if(ncores > 1){
    sfInit(cpus = ncores, parallel = TRUE)
    data("finalWeights")
    sfLibrary(geosphere)
    sfLibrary(raster)
    sfLibrary(automap)
    sfLibrary(mgcv)
    sfLibrary(SpatialDengue)
    sfExportAll()
    
    # if using default parameters
    if(missing(paramsList)){
      # if seasonal
      if(seasonal){
        tg <- function(x, weekdates, fitdat, pastdat, unipix, pixdistmat, steprun, seasonal_vector){
          print(x)
          return(DEN.spatial.seasonal(weekdates, fitdat, pastdat, sgpop, unipix, pixdistmat, steprun,
                                      seasonal_vector = seasonal_vector,
                                      seasonal_start = 7 * weekdates[1]))
        }
        
        enmodlist <- sfLapply(as.list(1:nruns), tg,
                              weekdates = weekdates,
                              fitdat = fitdat,
                              pastdat = pastdat,
                              unipix = unipix,
                              pixdistmat = pixdistmat,
                              steprun = steprun,
                              seasonal_vector = seasonal_vector)
        sfStop()
        
      }else{
        tg <- function(x, weekdates, fitdat, pastdat, unipix, pixdistmat, steprun){
          print(x)
          return(DEN.spatial(weekdates, fitdat, pastdat, unipix, pixdistmat, steprun))
        }
        
        enmodlist <- sfLapply(as.list(1:nruns), tg,
                              weekdates = weekdates,
                              fitdat = fitdat,
                              pastdat = pastdat,
                              unipix = unipix,
                              pixdistmat = pixdistmat,
                              steprun = steprun)
        sfStop()
      }
      
      
# if using custim parameters
    }else{
      
      # if seasonal 
      if(seasonal){
        tg <- function(x, weekdates, fitdat, pastdat, unipix, pixdistmat, steprun, seasonal_vector, paramsList){
          print(x)
          return(DEN.spatial.seasonal(weekdates, fitdat, pastdat, sgpop, unipix, pixdistmat, steprun,
                                      seasonal_vector = seasonal_vector,
                                      seasonal_start = 7 * weekdates[1],
                                      paramsList = paramsList))
        }
        
        enmodlist <- sfLapply(as.list(1:nruns), tg,
                              weekdates = weekdates,
                              fitdat = fitdat,
                              pastdat = pastdat,
                              unipix = unipix,
                              pixdistmat = pixdistmat,
                              steprun = steprun,
                              seasonal_vector = seasonal_vector,
                              paramsList = paramsList)
        sfStop()
      }else{
        tg <- function(x, weekdates, fitdat, pastdat, unipix, pixdistmat, steprun, paramsList){
          print(x)
          return(DEN.spatial(weekdates, fitdat, pastdat, unipix, pixdistmat, steprun, paramsList))
        }
        enmodlist <- sfLapply(as.list(1:nruns), tg,
                              weekdates = weekdates,
                              fitdat = fitdat,
                              pastdat = pastdat,
                              unipix = unipix,
                              pixdistmat = pixdistmat,
                              steprun = steprun,
                              paramsList = paramsList)
        sfStop()
      }
      
    }

  }else{
    # if not parallel:
    for(i in 1:nruns){
      if(missing(paramsList)){
        if(seasonal){
          enmodlist[[i]] <- DEN.spatial.seasonal(weekdates, fitdat, pastdat, unipix, pixdistmat, steprun, seasonal_vector)
        }else{
          enmodlist[[i]] <- DEN.spatial(weekdates, fitdat, pastdat, unipix, pixdistmat, steprun)
        }
        
      }else{
        if(seasonal){
          enmodlist[[i]] <- DEN.spatial.seasonal(weekdates, fitdat, pastdat, unipix, pixdistmat, steprun, seasonal_vector, paramsList)
        }else{
          enmodlist[[i]] <- DEN.spatial(weekdates, fitdat, pastdat, unipix, pixdistmat, steprun, paramsList)
        }
      }
    }
  }


  for(i in 1:nruns){

    # population totals
    Smat[, i] = enmodlist[[i]][[1]][, 1]
    Imat[, i] = enmodlist[[i]][[1]][, 2]
    Rmat[, i] = enmodlist[[i]][[1]][, 3]
    Rtmat[, i] = enmodlist[[i]][[1]][, 4]
    Dmat[, i] = enmodlist[[i]][[1]][, 5]
    newDmat[, i] = enmodlist[[i]][[1]][, 6]
    mosEmat[, i] = enmodlist[[i]][[1]][, 7]
    mosImat[, i] = enmodlist[[i]][[1]][, 8]
  }
  # summarising each compartment
  Sqs <- t(apply(Smat, 1, quantile, probs = qtiles, na.rm = T))
  Iqs <- t(apply(Imat, 1, quantile, probs = qtiles, na.rm = T))
  Rqs <- t(apply(Rmat, 1, quantile, probs = qtiles, na.rm = T))
  Rtqs <- t(apply(Rtmat, 1, quantile, probs = qtiles, na.rm = T))
  Dqs <- t(apply(Dmat, 1, quantile, probs = qtiles, na.rm = T))
  newDqs <- t(apply(newDmat, 1, quantile, probs = qtiles, na.rm = T))
  mosEqs <- t(apply(mosEmat, 1, quantile, probs = qtiles, na.rm = T))
  mosIqs <- t(apply(mosImat, 1, quantile, probs = qtiles, na.rm = T))

  # listing the outputs
  rtnlist <- as.list(rep(NA, length(qtiles)))
  for(i in 1:length(qtiles)){
    rtnlist[[i]] = data.frame(S = Sqs[, i],
                              I = Iqs[, i],
                              R = Rqs[, i],
                              Rt = Rtqs[, i],
                              D = Dqs[, i],
                              newD = newDqs[, i],
                              mosE = mosEqs[, i],
                              mosI = mosIqs[, i])
  }
  names(rtnlist) = sapply(qtiles, function(x) paste("Q_", x, sep = ""))
  # plotting
  plot(rtnlist$Q_0.5$newD, type = "l", ylim = c(0, max(rtnlist$Q_0.975$newD)),
       xlab = "Days", ylab = "Cases per day")
  polygon(c(1:steprun, steprun:1), c(rtnlist$Q_0.975$newD, rev(rtnlist$Q_0.025$newD)),
          col = rgb(0,0,1,0.5))
  polygon(c(1:steprun, steprun:1), c(rtnlist$Q_0.75$newD, rev(rtnlist$Q_0.25$newD)),
          col = rgb(0,0,1,0.8))
  lines(rtnlist$Q_0.5$newD, lwd = 2)

  return(rtnlist)
}
