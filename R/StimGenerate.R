#' Generating a baseline landscape of susceptible individuals
#'
#' Generates a continuous landscape of susceptible individuals given key parameters on the median and range of dengue seroprevalence and the past distribution of recent cases
#' @param pastdat Data frame of the locations, numbers and timings (in weeks) of recent cases, see ?sgdat
#' @param sero a three element vector giving the minimum, median and maximum observed seroprevalence in the landscape
#' @param startweek the week from which the model simulation will begin
#' @param unipix Universal pixel lookup table, see ?make.unipix
#' @details This approach assumes that true seroprevalence will be higher in areas with recent dengue infection. It therefore, takes past data from before the model simulation, uses krigeing (see ?autoKrig) to generate a continuosu surface of predicted past cases, then scales these to match the range and median of the seroprevalence values provided.
#' @keywords susceptible krige landscape
#' @export
#' @examples
#' data(sgdat)
#' data(sgpop)
#' sgpop <- pop.process(sgpop, agg = 10)
#' unipix <- make.unipix(sgpop)
#' sgdat <- data.frame(sgdat, patchID = apply(cbind(sgdat[, 3:2]), 1, pix.id.find, unipix))
#' sero <- c(0.577, 0.659, 0.814)
#' startweek <- 40
#' stim <- stim.generate(sgdat, sero, startweek, unipix)
#' plot.state(stim, sgpop, unipix)

stim.generate <- function(pastdat, sero, startweek, unipix){
  # trim case data to just before the startweek of the model simulation
  sgpast <- pastdat[pastdat$Week < startweek, ]
  # aggregate all past cases by each patch
  sgpastpatch <- aggregate(sgpast$Number_of_cases, by = list(sgpast$patchID), sum)
  
  # add latitude and longitude to the past data
  sgpastpatch <- data.frame(sgpastpatch, Longitude = rep(NA, nrow(sgpastpatch)),
                            Latitude = rep(NA, nrow(sgpastpatch)))
  for(i in 1:length(sgpastpatch[, 1])){
    sgpastpatch[i, 2] = sgpastpatch[i, 2] / unipix[match(sgpastpatch[i, 1], unipix$patchID), "pop"]
    sgpastpatch[i, 3] = sgpast[match(sgpastpatch[i, 1], sgpast$patchID), "Longitude"]
    sgpastpatch[i, 4] = sgpast[match(sgpastpatch[i, 1], sgpast$patchID), "Latitude"]
  }
  names(sgpastpatch)[2] = "z"
  newsgdat <- as.data.frame(coordinates(sgpop))
  coordinates(newsgdat) = ~x+y
  coordinates(sgpastpatch) = ~ Longitude+Latitude
  
  # use Krigeing to interpolate past cases in the landscape
  suppressWarnings(z <- capture.output(predIs <- autoKrige(z~1, sgpastpatch, newsgdat)))
  # only make predictions where people live
  stim = predIs$krige_output$var1.pred[(!is.na(as.vector(sgpop)) & as.vector(sgpop) > 0)]
  # make sure the interpolation doesn't predict negative seroprevalence
  stim[stim < 0] = 0
  
  # scale by median
  stim = stim/(median(stim) / sero[2])
  # scale upper values to be no higher than serohigh
  stimhigh <- stim[stim > sero[2]]
  stimhigh = (sero[3] - sero[2]) * ((stimhigh - min(stimhigh)) / max(stimhigh - min(stimhigh))) + sero[2]
  stim[stim > sero[2]] = stimhigh
  # scale lower values to be no lower than serolow
  stimlow <- stim[stim < sero[2]]
  stimlow = sero[1] + (sero[2] - sero[1]) * ((stimlow - min(stimlow)) / max(stimlow - min(stimlow)))
  stim[stim < sero[2]] = stimlow
  
  return(stim)
}