#' Updating sample size calculations for trial design
#'
#' 
#' @param modresult a Spatial dengue model object returned from DEN.spatial
#' @param radi radius around the index patch that infections / cases are observed within
#' @param pixdistmat A patch distance matrix, see example
#' @param adultsonly logical, will the trial only include adults (= TRUE), or all ages (= FALSE)
#' @param foi single values range 0-1, long term average force of infection in the area concerned, usually estimated from age-stratified seroprevalence surveys, optional, 
#' defaults to 0.1 ~ an endemic setting, only required if adultsonly = TRUE
#' @param agedist optional vector of counts of individuals in 1 year age bands in years 0-100, defaults to an exponential age distribution otherwise
#' @param dly single value, Days between observign symptomatic case and visiting index case's home, defaults to 0 (day of becoming symptomatic)
#' @details returns four plots showing: i) total people within the specified radius of each index cases, ii) total seronegative people within the specified radius of each index case,
#' iii) infection incidence within the radius of each index case in the 28 days following the visit, iv) symptomatic disease incidence in the radius of the index case in the 28 days following the visit
#' Each of these measures are plotted over time over the course of the outbreak with red splines showing the estimated median value.
#' A data.frame is also returned summarising these distributions over the course of the experiment
#' @keywords sample size
#' @export
#' @examples
#' data(sgdat)
#' data(sgpop)
#' sgpop <- pop.process(sgpop, agg = 10)
#' unipix <- make.unipix(sgpop)
#' pixdistmat <- distm(cbind(unipix$x, unipix$y))
#' sgdat <- data.frame(sgdat, patchID = apply(cbind(sgdat[, 3:2]), 1, pix.id.find, unipix))
#' weekdates <- c(40, 92)
#' # model run with default parameters
#' denmod_sim <- DEN.spatial(weekdates, sgdat, unipix, pixdistmat)
#' set.seed(123)
#' Sample.size.calc(denmod_sim, radi = 1000, pixdistmat, adultsonly = FALSE)
#' # now just adults only
#' set.seed(123)
#' Sample.size.calc(denmod_sim, radi = 1000, pixdistmat, adultsonly = TRUE, foi = 0.12)

# updating sample size calculations:
# dly = days between observign symptomatic case and visiting index case's home
# radi = radius around index patch that are observed
# adults only (logical) do you just enroll adults?
# foi = force of infection (only needed if only looking at adults)
# agedist = adult age distribution (only needed if only looking at adults)
# just observe people in the same patch immediately (all ages)





Sample.size.calc <- function(modresult, radi, pixdistmat, adultsonly = FALSE, foi = 0.1, agedist = NA, dly = 0){
  daymod <- modresult[[2]]
  # add cumulative infections to the individual patch totals
  ctot = rep(0, nrow(daymod[[1]]))
  for(i in 1:length(daymod)){
    ctot = ctot + daymod[[i]]$newD
    daymod[[i]]$cumnewD = ctot
  }
  
  # assign index case data frame
  indcase <- matrix(0, nrow= 0, ncol = 6)
  colnames(indcase) = c("Totsuscep", "Adultsuscep", "NextInf", "NextCase", "TimeVisit", "Tpop")
  
  for(i in 1:(length(daymod) - dly - 28)){
    # bring up index case patches
    APindex = find.Dspots(newsingD = matrix(daymod[[i]]$newD, nrow = 1), pixdistmat, radius = radi)
    
    if(!is.na(APindex[1])){
      indcaseT = matrix(0, nrow = length(APindex), ncol = 6)
      # how many susceptibles will there be in i + dly time?
      indcaseT[, 1] = daymod[[(i + dly)]][APindex, 1]
      # how about susceptible adults? !!! TO DO !!!
      
      # how many infections will occur in these index patches
      # between days (i + dly) and (i + dly+ 28)?
      # N.B. infetion assumed to be detected by point estimate IgG-i.e. recovered individuals
      # Now also assuming IgM and PCR to remove those infected but not yet recovered
      #indcaseT[, 3] = daymod[[(i + dly + 28)]][APindex, 3] - daymod[[(i + dly)]][APindex, 3]
      indcaseT[, 3] = daymod[[(i + dly + 28)]][APindex, 3] - daymod[[(i + dly)]][APindex, 3] - daymod[[(i + dly)]][APindex, 2]
      
      # and how many cases (passive surveillance) over the same time period?
      #indcaseT[, 4] = daymod[[(i + dly + 28)]][APindex, 9] - daymod[[(i + dly)]][APindex, 9]
      # disease per infection ratio:
      Disratio = (daymod[[(i + dly + 28)]][APindex, 9] - daymod[[(i + dly)]][APindex, 9]) /
        (daymod[[(i + dly + 28)]][APindex, 3] + daymod[[(i + dly + 28)]][APindex, 2] - daymod[[(i + dly)]][APindex, 3])
      indcaseT[, 4] = round(indcaseT[, 3] * Disratio, 0)
      
      # asign the time of visit
      indcaseT[, 5] = i + dly
      
      # asign total patch population
      indcaseT[, 6] = as.vector(rowSums(daymod[[(i + dly + 28)]][APindex, 1:3]))
      
      #
      
      # add to main collector
      indcase = rbind(indcase, indcaseT)
    }
  }
  
  # adjustment if just enrolling adults
  if(adultsonly){
    # unless age distribution supplied, assume a rectangle population pyramid
    if(is.na(agedist[1])){
      agedist = 1/exp((0:100) / 100)}
      a = 0:100
    
    # probability of priamry and secondary infection
    pinf <- exp(-foi*a)
    p2inf <- (1 - exp(-foi*a)) * (exp((0.75 * -foi)*a))
    
    # proportion of all susceptibles that are children
    child = sum((agedist * pinf)[1:17]) / sum((agedist * pinf))
    
    # proportion of all disease caes that are in children
    child_D = sum((agedist * p2inf)[1:17]) / sum((agedist * p2inf))
    
    # adjust indcase
    indcase[, 1] = indcase[, 1] * (1 - child)
    indcase[, 3] = indcase[, 3] * (1 - child)
    indcase[, 4] = indcase[, 4] * (1 - child_D)
    indcase[, 6] = indcase[, 6] * sum(agedist[18:81]) 
  }
  
  
  par(mfrow = c(2,2))
  
  # 01 plot people per index patch
  # how many people do you have to test serostatus for per index case?
  plot(indcase[, 5], indcase[, 6], pch = 19, cex = 0.5, xlab = "Day", 
       main = "People per index case", ylab = "n")
  linetab = data.frame(y = indcase[, 6], x = indcase[, 5])
  sline = gam(y ~ s(x), data = linetab)
  spred = predict(sline, newdata = data.frame(x = 0:365), type = "response",
                  se.fit = TRUE)
  lines(0:365, spred$fit, col = "red", lwd = 2)
  lines(0:365, spred$fit + 2*spred$se.fit, col = "red", lwd = 2, lty = 2)
  lines(0:365, spred$fit - 2*spred$se.fit, col = "red", lwd = 2, lty = 2)
  
  # 02 plot total susceptibles in index case patch over time:
  # what is your max enrollment number?
  plot(indcase[, 5], indcase[, 1], pch = 19, cex = 0.5, xlab = "Day", 
       main = "Seronegatives per index case", ylab = "n")
  linetab = data.frame(y = indcase[, 1], x = indcase[, 5])
  sline = gam(y ~ s(x), data = linetab)
  spred = predict(sline, newdata = data.frame(x = 0:365), type = "response",
                  se.fit = TRUE)
  lines(0:365, spred$fit, col = "red", lwd = 2)
  lines(0:365, spred$fit + 2*spred$se.fit, col = "red", lwd = 2, lty = 2)
  lines(0:365, spred$fit - 2*spred$se.fit, col = "red", lwd = 2, lty = 2)
  
  
  # 03 how many infections will occur in these patches over the next 28 days?
  plot(indcase[, 5], indcase[, 3] / indcase[, 1], pch = 19, cex = 0.5, xlab = "Day", 
       main = "Infection incidence", ylab = "p")
  linetab = data.frame(y = indcase[, 3] / indcase[, 1], x = indcase[, 5])
  sline = gam(y ~ s(x), data = linetab)
  spred = predict(sline, newdata = data.frame(x = 0:365), type = "response",
                  se.fit = TRUE)
  lines(0:365, spred$fit, col = "red", lwd = 2)
  lines(0:365, spred$fit + 2*spred$se.fit, col = "red", lwd = 2, lty = 2)
  lines(0:365, spred$fit - 2*spred$se.fit, col = "red", lwd = 2, lty = 2)
  
  # 04 how many passively detected cases will occur in these patches over the next 28 days?
  dstore = cbind(indcase[indcase[, 1] > 0, 5], indcase[indcase[, 1] > 0, 4] / indcase[indcase[, 1] > 0, 1])
  dstore = dstore[!is.na(dstore[, 2]), ]
  dstore = dstore[dstore[, 2] <= 1, ]
  
  plot(dstore, pch = 19, cex = 0.5, xlab = "Day", 
       main = "Disease incidence", ylab = "p")
  linetab = data.frame(y = dstore[, 2], x = dstore[, 1])
  sline = gam(y ~ s(x), data = linetab)
  spred = predict(sline, newdata = data.frame(x = 0:365), type = "response",
                  se.fit = TRUE)
  lines(0:365, spred$fit, col = "red", lwd = 2)
  lines(0:365, spred$fit + 2*spred$se.fit, col = "red", lwd = 2, lty = 2)
  lines(0:365, spred$fit - 2*spred$se.fit, col = "red", lwd = 2, lty = 2)
  dstore = as.vector(summary(dstore[, 2]))
  
  # return some numerical results
  
  
  rtndf <- rbind(rep(nrow(indcase), 6),
                 as.vector(summary(indcase[, 6])),
                 as.vector(summary(indcase[, 1])),
                 as.vector(summary(indcase[indcase[, 1] > 0, 3] / indcase[indcase[, 1] > 0, 1])),
                 dstore)
  colnames(rtndf) = names(summary(indcase[, 1]))
  rownames(rtndf) = c("Total index cases",
                      "People per index case",
                      "Seronegatives per index case",
                      "Infection incidence",
                      "Disease incidence")
  return(rtndf)
}