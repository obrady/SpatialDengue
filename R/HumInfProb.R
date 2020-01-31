#' Calculate probability of human infection
#'
#' Calculates the pixel-specific probability of human infection taking into account their time allocation in different pixels with or without infectious mosquitoes
#' @param mos_I vector, Infectious mosquitoes in each pixel
#' @param betaEnv Starting human <-> mosquito contact rate landscape, see ?betaEnv.generate
#' @param unipix Universal pixel lookup table, see ?make.unipix
#' @param mm_list2 List of human movement matrices adjusted by the hometime parameter, see ?move.matrix.gene and ?DEN.spatial
#' @param mm_list_cum Matrix of patch effective population size (including both resident population and those who move there during the day)
#' @details Calculates proabbility of infection of mosquitoes resident in each patch (mosquitoes don't move around) when infectious humans are moving around.
#' Infectiousness of humans is also dependent on their time since infection, peaking when symptomatic.
#' Returns a vector with length equal to the number of patches giving the proabbility of infection of each susceptible mosquito
#' @keywords mosquito infection movement
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
#' betaEnv <- betaEnv.generate(bmean = 8, blogvariance = log(2), bcorrelation = -1, stim, unipix)
#' mos_I <- runif(nrow(unipix), min = 0, max = 20)
#' mlist <- move.matrix.gene(unipix, hometime = 0.2)
#' HIP <- hum.inf.prob(mos_I, betaEnv, unipix, mlist[[1]], mlist[[2]])
#' plot.state(HIP, sgpop, unipix)






hum.inf.prob <- function(mos_I, betaEnv, unipix, mm_list2, mm_list_cum){
  # collate info for patches with infectious mosquitoes
  infmos_I_i = data.frame(unipix[mos_I > 0, ], mosnum = mos_I[mos_I > 0])
  # effective number of people in each pixel (only those with I mosquitoes)- taking into account travel
  popEff <- mm_list_cum[infmos_I_i[, 1]]
  # infectious mosquito bites in each patch (only those with I mosquitoes)
  mosbite  = infmos_I_i$mosnum * betaEnv[infmos_I_i[, 1]]
  # infectious mosquito bites are distributed over popEff number of people in each patch
  # what is the prob of infection in each patch takign into account all infection sources?
  humInfProb = rep(0, nrow(unipix))
  for(k in 1:nrow(infmos_I_i)){
    humInfProb = humInfProb + mm_list2[[infmos_I_i$patchID[k]]] * mosbite[k] / popEff[k]
  }
  humInfProb[humInfProb > 1] = 1
  return(humInfProb)
}