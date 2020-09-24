#' Establishing model initial conditions
#'
#' Generates the starting values for each of the model states given reported cases at time = 0
#' @param casedat Data frame of the locations, numbers and timings (in weeks) of recent cases, see ?sgdat
#' @param startweek Time (in weeks) at which the model simulation is to begin
#' @param unipix Universal pixel lookup table, see ?make.unipix
#' @param stim Starting susceptibiltiy landscape, see ?stim.generate
#' @param betaEnv Starting human <-> mosquito contact rate landscape, see ?betaEnv.generate
#' @param mospdeath single value between 0 and 1, the daily probability of mosquito mortality
#' @param hometime Parameter that controls spatial concentration of movement, defaults to 0.2
#' @details Uses supplied relationships to work out how many Infectious, exposed, etc are needed over time to generate the supplied reported detected human cases at t= 0.
#' These are then used for model burin.
#' @keywords burnin generation time
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
#' betaEnv <- betaEnv.generate(bmean = 8, blogvariance = log(2), bcorrelation = -1, stim)
#' mospdeath <- 0.1
#' startTable <- case.backtrack.time(sgdat, startweek, unipix, stim, betaEnv, mospdeath, hometime = 0.5)


case.backtrack.time <- function(casedat, startweek, unipix, stim, betaEnv, mospdeath, hometime = 0.2){
  # generating start table of reported case locations at t = 0
  sg_start_inf <- casedat[casedat$Week == startweek, ]

  # extract the relevant pixel IDs for infected start points
  star_inf_ID <- apply(cbind(sg_start_inf$Longitude, sg_start_inf$Latitude), 1, pix.id.find, unipix = unipix)
  star_inf_ID <- data.frame(ID = star_inf_ID, Cases = sg_start_inf$Number_of_cases)
  star_table <- aggregate(star_inf_ID$Cases, by = list(star_inf_ID$ID), sum)
  colnames(star_table) = c("ID", "Count")

  # at which timestep were they infected?, assume infected in current location?
  IIP_samples <- round(rlnorm(round(sum(star_table[, 2]) / 7, 0), meanlog = log(5.9), sdlog = 0.045), 0)
  star_table$infhumtime = -median(IIP_samples)
  # assume only infectious mosquitoes in locale of cases
  star_table$sing_S = round(((1 - stim) * unipix$pop)[star_table$ID], 0)
  star_table$sing_pop = unipix$pop[star_table$ID]

  # how many infectious mosquitoes were there in each of these patches at the start of the burnin period?
  star_table$infmos = round(star_table$Count / (7 * betaEnv[star_table$ID] * hometime * star_table$sing_S / star_table$sing_pop), 0)
  star_table$infmos[!is.finite(star_table$infmos)] = max(star_table$infmos[is.finite(star_table$infmos)])

  # probability a mosquito survives EIP
  deadm_vec <- exp(-rep(mospdeath, sum(star_table$infmos)) * 1)
  EIPsample <- rlnorm(sum(star_table$infmos), meanlog = log(7), sdlog = 0.21)
  post <- deadm_vec^EIPsample

  # add amount of exposed mosquitoes and when they were relative to a reported case
  star_table$expmos = round(star_table$infmos / median(post), 0)
  star_table$expmostime = star_table$infhumtime - median(round(EIPsample, 0))
  return(star_table)
}
