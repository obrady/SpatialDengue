#' Creating an evaluation dataset for fitting
#'
#' Processes the case data to generate a list of metrics required for different aspects of evaluating model fit
#' @param fitdat Casedata, see ?sgdat for format
#' @param weekdates Two element vector of the start and end weeks of the simulation over which the mdoel will be evaluated over
#' @param unipix Universal pixel lookup table, see ?make.unipix
#' @details returns a list with three elements: i) weekly detialed number of cases per patch,
#' ii) Moran's weights for each pixel (for calculating Morans I)
#' iii) Weekly cumulative total cases
#' @keywords evaluation
#' @export
#' @examples
#' data(sgdat)
#' data(sgpop)
#' sgpop <- pop.process(sgpop, agg = 10)
#' unipix <- make.unipix(sgpop)
#' sgdat <- data.frame(sgdat, patchID = apply(cbind(sgdat[, 3:2]), 1, pix.id.find, unipix))
#' EvalSet <- make.eval.set(sgdat, weekdates = c(40, 92), unipix)
#' summary(EvalSet)

make.eval.set <- function(fitdat, weekdates, unipix){
  
  # 01 Weekly_pixel_detail
  eval_list = list()
  weekpoints <- (weekdates[1] + 1):weekdates[2]
  counter = 1
  for(i in weekpoints){
    # current week- for spatial autocorrelation eval
    sgweek = fitdat[fitdat$Week == i, ]
    if(nrow(sgweek) > 0){
      cases_patch <- aggregate(sgweek$Number_of_cases, by = list(sgweek$patchID), sum)
      cases_patch2 = rep(0, nrow(unipix))
      for(k in 1:nrow(cases_patch)){
        cases_patch2[cases_patch[k, 1]] = cases_patch[k, 2]
      }
      # corresponding week with readjusted start time
      cor_week = (i - (weekdates[1])) * 7 +1
      
      # reference evaluation table
      eval_tab = data.frame(cases = cases_patch2, day = cor_week, unipix[, c("x", "y")])
      
    } else{eval_tab = NA}
    eval_list[[counter]] = eval_tab
    counter = counter + 1
  }
  
  # 02 Moran's weighting per pixel
  coords = eval_list[[which.max(!is.na(eval_list))]][, c("x", "y")]
  radius = 0.04504505 # only interested in spatial correlation out to 5km
  pixdists <- as.matrix(dist(coords))
  pixdists.inv <- 1 / pixdists
  diag(pixdists.inv) <- 0
  pixdists.bin <- (pixdists > 0 & pixdists <= radius)
  pixdists.inv = pixdists.inv * pixdists.bin
  
  # 03 Cumulative case count by week
  casesWeek = rep(NA, length(weekdates[1]:weekdates[2]))
  wdvec = weekdates[1]:weekdates[2]
  for(i in 1:length(weekdates[1]:weekdates[2])){
    casesWeek[i] = sum(fitdat[fitdat$Week == wdvec[i], "Number_of_cases"])
  }
  casesWeek = cumsum(casesWeek)
  
  # geenrate return object
  rtnlist <- list(eval_list, pixdists.inv, casesWeek)
  
  names(rtnlist) <- c("Weekly_pixel_detail", "Morans_weights", "Week_Cum_cases")
  return(rtnlist)
}