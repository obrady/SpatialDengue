#' Generating a baseline landscape of vector-human and human-vector contact rates
#'
#' Generates a continuous landscape of human <-> vector contact rates given parameters of mean and variance of contact rate and the degree of spatial correlation with susceptibility
#' @param bmean the mean human <-> vector contact rate across the landscape
#' @param bcorrelation the degree of correlation between human susceptibility and human <-> vector contact rate, must be in interval -1 to 1
#' @param stim a vector of the starting immunity levels in each pixel
#' @details Between the three supplied parameters this controls the intensity of transmission, its spatial variability and its spatial structure (in particularlly whether transmission intensity is most intense in areas that have experience recent transmission, or areas that have not experience recent transmission), i.e. does the outbreak happen in a new area or an existign area of high intensity transmission? The implementation is stochastic and an approximation.
#' @keywords contact rate mosquito human
#' @export
#' @examples
#' data(sgdat)
#' data(sgpop)
#' sgpop <- pop.process(sgpop, agg = 10)
#' unipix <- make.unipix(sgpop)
#' sgdat <- data.frame(sgdat, patchID = apply(cbind(sgdat[, 3:2]), 1, pix.id.find, unipix))
#' sero <- c(0.577, 0.659, 0.814)
#' startweek <- 40
#' stim <- stim.generate(sgdat, sero, startweek, sgpop, unipix)
#' par(mfrow = c(2, 1))
#' plot.state(stim, sgpop, unipix)
#' betaEnv <- betaEnv.generate(bmean = 8, bcorrelation = -1, stim)
#' plot.state(betaEnv, sgpop, unipix)


betaEnv.generate <- function(bmean, bcorrelation, stim){
  # assume sd of stim = sd of betaEnv
  # function for generating vector of correlated values
  correlatedValue = function(x, r){
    r2 = r**2
    ve = 1-r2
    SD = sqrt(ve)
    e  = rnorm(length(x), mean=0, sd=SD)
    y  = r*x + e
    return(y)
  }
  
 if (sd(stim)==0){
    
     betaEnv =  rep(bmean, length(stim))

  } else{ 
  
   # Z scores of starting immunity
  stim_Zscore = (stim - mean(stim)) / sd(stim)
  
  # calculate correlated betaEnv Z scores
  betaEnv_Zscore = correlatedValue(x=stim_Zscore, r=bcorrelation)
  
  # and back to original scale
  betaEnv = betaEnv_Zscore * sd(stim) + bmean }
  
  # scaling for below 0 values
  if(any(betaEnv < 0)){
    betaEnv2 = betaEnv
    betaEnv2[betaEnv2 < 0] = 0
    # now scale non 0 values to match the original mean
    betaEnv2 = betaEnv2 * mean(betaEnv) / mean(betaEnv2)
    betaEnv = betaEnv2
  }
  return(betaEnv)
}
