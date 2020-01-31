#' Run a model simulation (workhorse function)
#'
#' Directly implements the daily simulation of the dengue spatial model, usually used nested within a wrapper function, see ?DENSpatial
#' @param sing_S vector, Susceptible human individuals in each pixel
#' @param sing_I matrix, Infectious human individuals in each pixel, rows = day of infection
#' @param sing_R vector, Recovered human individuals in each pixel
#' @param sing_Rt matrix, Recovered (temporary) human individuals in each pixel, rows = day of protection
#' @param mos_S vector, Susceptible mosquitoes in each pixel
#' @param mos_E matrix, Exposed but not yet infectious mosquitoes in each pixel, rows = day since infection
#' @param mos_I vector, Infectious mosquitoes in each pixel
#' @param mos_S vector, Susceptible mosquitoes in each pixel
#' @param betaEnv Starting human <-> mosquito contact rate landscape, see ?betaEnv.generate
#' @param mm_list2 List of human movement matrices adjusted by the hometime parameter, see ?move.matrix.gene and ?DEN.spatial
#' @param mm_list_cum Matrix of patch effective population size (including both resident population and those who move there during the day)
#' @param pdetec Function for sampling probability of detecting a denge ifected individual (by passive surveillance), beta distributied probability
#' @param func_IIP Function for probability of completing human Intrinsic Incubation Period (time between infection and becoming symptomatic), a funciton of time since infection
#' @param func_EIP Function for probability of completing mosquito Extrinsic Incubation Period (time between mosquito becoming infected and becoming infectious), a function of time since infection
#' @param mospdeath Function for daily probability of death of a mosquito, beta distributed hazard, exponential distributed survival
#' @param drugtreat_Eff Effective coverage of prophylactic drug
#' @param drugtreat_Dur Duration of prophylactic drug effectiveness
#' @param drugtreat_Rad Radius of deployment around an index case of prophylactic drug
#' @param drugtreat_Del Delay in days between detection of a dengue case and treatment of target population with drugs
#' @param vectreat_Eff Effective coverage of vector control
#' @param vectreat_Dur Not used
#' @param vectreat_Rad Radius of deployment around an index case of vector control
#' @param vectreat_Del Delay in days between detection of a dengue case and treatment of target population with vector control
#' @param vectreat three element vector of details governing routine reactive vector control, see tutorial
#' @param steps Number of days for which the simulation shoudl run, defaults to 365
#' @param unipix Universal pixel lookup table, see ?make.unipix
#' @param pixdistmat A patch distance matrix, see example in ?DENspatial
#' @param StopCriteria Single numeric value, if the model predicts reported cases exceed this value at any point during the simulation, simulation is terminated and "NA" is returned
#' @details Implements a stochastic daily timestep model predictign the spatial spread of dengue across Singapore. As long as "StopCriteria" are not exceeded, will return a list of final values of each of the 
#' main stages of the model and three summaries: i) "totals" - population wide totals for each model state and ii) detailed_totals- pixel-specific totals for each model state, iii)
#' treatlog - log of length "steps" detailing which patches were treated with drugs or vector control on each day
#' @keywords model workhorse simulation
#' @export
#' @examples
#' See tutorial
model.run <- function(sing_S,
                      sing_I,
                      sing_R,
                      sing_Rt,
                      mos_S,
                      mos_E,
                      mos_I,
                      betaEnv,
                      mm_list2,
                      mm_list_cum,
                      pdetec,
                      func_IIP,
                      func_EIP,
                      mospdeath,
                      drugtreat_Eff,
                      drugtreat_Dur,
                      drugtreat_Rad,
                      vectreat_Eff,
                      vectreat_Dur,
                      vectreat_Rad,
                      steps = 365,
                      unipix,
                      pixdistmat,
                      StopCriteria){
  
  # set up vectors to hold results
  totals = data.frame(S = rep(NA, steps),
                      I = rep(NA, steps),
                      R = rep(NA, steps),
                      Rt = rep(NA, steps),
                      D = rep(NA, steps),
                      newD = rep(NA, steps),
                      mosE = rep(NA, steps),
                      mosI = rep(NA, steps))
  hold = data.frame(S = rep(NA, nrow(unipix)),
                    I = rep(NA, nrow(unipix)),
                    R = rep(NA, nrow(unipix)),
                    Rt = rep(NA, nrow(unipix)),
                    D = rep(NA, nrow(unipix)),
                    newD = rep(NA, nrow(unipix)),
                    mosE = rep(NA, nrow(unipix)),
                    mosI = rep(NA, nrow(unipix)))
  detailed_totals <- list()
  for(i in 1:steps){detailed_totals[[i]] = hold}
  
  # tracks which patches have been treated at each timestep
  treatlog_drug <- list()
  treatlog_vec <- list()
  
  # add detected state
  sing_D <- matrix(0, nrow = nrow(sing_I), ncol = ncol(sing_I))
  # total population size
  tpop = sum(sing_S, sing_I, sing_R, sing_Rt)
  
  # main timestep loop
  for(i in 1:steps){
    # SIR-SEI formulation
    # Part 1 calculates stage transition probabilities
    # Part 2 updates states
    
    # Part 1 stage transition probabilities
    
    # 1A- P(mosquito infected)
    mosInfProb <- mos.inf.prob(sing_I, betaEnv, unipix, mm_list2, mm_list_cum)
    # 1B- P(mosquito completes EIP)
    # func_EIP(days)
    # 1C- P(mosqutio dies naturally)
    mosDeathProb <- mospdeath(nrow(unipix))
    # 1D- calculate drug targetted patches and P(mosquto dies from RVC) and P(human treated with drugs)
    # don't bother doing if Effective coverage = 0
    if(sum(vectreat_Eff, drugtreat_Eff) != 0){
      # identify detected locations
      if(i == 1){newsingDs = matrix(0, nrow = nrow(sing_D), ncol = ncol(sing_D))}
      Dspots1 <- find.Dspots(newsingDs, pixdistmat, radius = drugtreat_Rad)
      Dspots2 <- find.Dspots(newsingDs, pixdistmat, radius = vectreat_Rad)
      mosDeathProbLocal <- rep(0, nrow(unipix))
      humDrugTreatProbLocal <- rep(0, nrow(unipix))
      # Drugs
      if(sum(!is.na(Dspots1)) > 0){
        # assign to treatment log
        treatlog_drug[[i]] = Dspots1
        
        # assign drugs (+/- delay between patient detection and local area visit)
        scheduleday = i - drugtreat_Del
        if(!all(is.na(treatlog_drug[[scheduleday]]))){
          humDrugTreatProbLocal[treatlog_drug[[scheduleday]]] <- drugtreat_Eff
        }
      }else{treatlog_drug[[i]] = NA}
      
      # Vector control
      if(sum(!is.na(Dspots2)) > 0){
        # assign to treatment log
        treatlog_vec[[i]] = Dspots2
        
        # assign drugs (+/- delay between patient detection and local area visit)
        scheduleday = i - vectreat_Del
        if(!all(is.na(treatlog_vec[[scheduleday]]))){
          mosDeathProbLocal[treatlog_vec[[scheduleday]]] <- vectreat_Eff
        }
      }else{treatlog_vec[[i]] = NA}
    }else{
      Dspots1 = NA
      Dspots2 = NA
      mosDeathProbLocal <- rep(0, nrow(unipix))
      humDrugTreatProbLocal <- rep(0, nrow(unipix))
      treatlog_drug[[i]] = NA
      treatlog_vec[[i]] = NA
    }
    # 1F- P(human infected)
    if(sum(mos_I) > 0){humInfProb <- hum.inf.prob(mos_I, betaEnv, unipix, mm_list2, mm_list_cum)} else{humInfProb = 0}
    
    # 1G- P(human completes IIP)
    #func_IIP(days)
    # 1H- P(symptomatic human detected)
    #pdetec(1)
    
    # Part 2 state transitions
    # MOSQUITO 
    #(S never changes so not here)
    # mos E
    mos_E_IO <- multinom(list(list(mos_S, mosInfProb),
                              list(mos_E, func_EIP),
                              list(mos_E,  mospdeath),
                              list(mos_E, mosDeathProbLocal)),
                         type = "mos_E")
    #mos_E = mos_E + mos_E_IO[[1]] - mos_E_IO[[2]] - mos_E_IO[[3]] - mos_E_IO[[4]]
    # extras includes those that reach the end of the holding vector but havent completes EIP or dies 
    # (so we automatically advance them to the next state)
    extraMosE_I = mos_E[nrow(mos_E), ] - mos_E_IO[[2]][nrow(mos_E), ] - mos_E_IO[[3]][nrow(mos_E), ] - mos_E_IO[[4]][nrow(mos_E), ]
    mos_E = transition(base = mos_E,  plus = mos_E_IO[[1]], minus = (mos_E_IO[[2]] + mos_E_IO[[3]] + mos_E_IO[[4]]))
    
    # mos_I
    mos_I_O <- multinom(list(list(mos_I,  mospdeath),
                             list(mos_I, mosDeathProbLocal)),
                        type = "mos_I")
    mos_I = mos_I + colSums(mos_E_IO[[2]]) - mos_I_O[[1]] - mos_I_O[[2]] + extraMosE_I
    
    # HUMAN
    # sing_S
    sing_S_IO <- multinom(list(list(sing_S, humInfProb),
                               list(sing_S, humDrugTreatProbLocal),
                               list(sing_Rt,  function(times) pnorm(times, drugtreat_Dur, 0))),
                          type = "sing_S")
    if(any((sing_S - sing_S_IO[[2]] - sing_S_IO[[3]]) < 0)){break}
    sing_S = sing_S - sing_S_IO[[2]] - sing_S_IO[[3]] + colSums(sing_S_IO[[1]])
    
    # sing I (only drop out at 8 days after symptom onset, so no multinomials needed)
    sing_I_O <- Hum.Drug.Treat.MN(sing_I, humDrugTreatProbLocal)
    extraHumI_R = sing_I[nrow(sing_I), ] - sing_I_O[nrow(sing_I), ]
    sing_I = transition(base = sing_I, plus = sing_S_IO[[2]], minus = sing_I_O)
    
    # sing_R
    sing_R = sing_R + colSums(sing_I_O) + extraHumI_R
    
    # sing_Rt
    extraHumR_S = sing_Rt[nrow(sing_Rt), ] - sing_S_IO[[1]][nrow(sing_Rt), ]
    sing_S = sing_S + extraHumR_S
    sing_Rt = transition(base = sing_Rt, plus = sing_S_IO[[3]], minus = sing_S_IO[[1]])
    
    # sing_D
    # who has completed their EIP and gets detected
    # first transition Sing_D so it matches sing_I
    sing_D = transition(sing_D, plus = 0, minus = 0)
    # only look at infecteds that have not already been detected
    sing_I2 = sing_I - sing_D
    #newsignDopts <- cbind(as.vector(sing_I2), rep(func_IIP(1:nrow(sing_I2)), ncol(sing_I2)) * pdetec(length(sing_I2)))
    newsignDopts <- cbind(as.vector(sing_I2), rep(func_IIP(1:nrow(sing_I2)), ncol(sing_I2)) * pdetec(1))
    poscells <- newsignDopts[, 1] > 0
    newsingDs = rep(0, length(newsignDopts[, 1]))
    if(sum(poscells) > 0){
      if(sum(poscells) == 1){newsingDs[poscells] = rbinom(1, newsignDopts[poscells, 1], newsignDopts[poscells, 2])}else{
        newsingDs[poscells] = apply(newsignDopts[poscells, ], 1, function(u) rbinom(1, u[1], u[2]))
      }
    }
    newsingDs = matrix(newsingDs, nrow = nrow(sing_D), ncol = ncol(sing_D))
    sing_D = sing_D + newsingDs
    
    
    
    
    # store results 
    # broad summary
    totals[i, ] = c(sum(sing_S), sum(sing_I), sum(sing_R), sum(sing_Rt), sum(sing_D), sum(newsingDs),
                    sum(mos_E), sum(mos_I))
    # detailed summary
    detailed_totals[[i]][, 1:8] = cbind(sing_S, colSums(sing_I), sing_R, colSums(sing_Rt), colSums(sing_D), colSums(newsingDs),
                                        colSums(mos_E), mos_I)
    
    # breaking if StopCriteria are breached (efficiency measure)
    if(sum(totals$newD[1:i]) > StopCriteria){
      return(NA)
    }
    
    ## mF5I- model checks
    # check no negative cases
    if(any(c(sing_S, sing_I, sing_R, sing_Rt, sing_D) < 0)){
      print(paste("failed at timestep", i, "due to negative people"))
      break}
    # check for NAs
    if(any(is.na(sing_S), is.na(sing_I), is.na(sing_R), is.na(sing_Rt), is.na(sing_D))){
      print(paste("failed at timestep", i, "due to NAs in population vectors"))
      break}   
    # check dimensions of dataframes are consistent
    if(any(c(nrow(sing_I), nrow(mos_E)) != 15)){
      print(paste("failed at timestep", i, "due to changing matrix size"))
      print(c(nrow(sing_I), nrow(mos_E)))
      break}
    # check human population size is consistent
    if(sum(sing_S, sing_I, sing_R, sing_Rt) != tpop){
      print(paste("failed at timestep", i, "due to changing human population size"))
      break}
  }
  # collate final results
  treatlog = list(treatlog_drug, treatlog_vec)
  
  templist <- list(sing_S,
                   sing_I,
                   sing_R,
                   sing_Rt,
                   sing_D,
                   newsingDs,
                   mos_S,
                   mos_E,
                   mos_I,
                   totals,
                   detailed_totals,
                   treatlog)
  names(templist) = c("sing_S",
                      "sing_I",
                      "sing_R",
                      "sing_Rt",
                      "sing_D",
                      "newsingDs",
                      "mos_S",
                      "mos_E",
                      "mos_I",
                      "totals",
                      "detailed_totals",
                      "treatlog")
  return(templist)
}
