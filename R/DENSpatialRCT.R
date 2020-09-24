#' Simulate a randomised control trial of dengue prophylactic drugs
#'
#' Provide basic information on the RCT and the number of stochastic runs over which model runs should be averaged
#' @param weekdates Two element vector of the start and end weeks of the simulation over which the mdoel will be evaluated over
#' @param fitdat Data frame of the locations, numbers and timings (in weeks) of cases to fit the model to, see ?sgdat
#' @param pastdat Data frame of the locations, numbers and timings (in weeks) of all cases in the dataset (is used to generate the starting immunity profile), see ?sgdat
#' @param unipix Universal pixel lookup table, see ?make.unipix
#' @param pixdistmat A patch distance matrix, see example
#' @param RCTinfo A data frame including details on RCT start and end time, control and treatment patches and effective coverage reached in treatment patch
#' @param nruns integer, number of stochastic runs of the model over which results should be averaged
#' @param paramsList Optional parameter list. If not supplied returns to defaults, see tutorial for full parameter list, see ?model.run for full list and explanation of parameters
#' @details runs "nruns" number of stochastic simulations and treats people in the trial treatment clusters at the specified time. Returns
#' a named list of six matrices detailing the number of seroconversions and symptomatic cases in treatment and control patches with each columns of each matrix representing one realisation from the stochastic model. 
#' The final two elements in the list detail the total person-days of observation and total number of people treated with drugs in the trial.
#' @keywords model simulation RCT
#' @export
#' @examples
#' data(sgdat)
#' data(sgpop)
#' sgpop <- pop.process(sgpop, agg = 10)
#' unipix <- make.unipix(sgpop)
#' pixdistmat <- distm(cbind(unipix$x, unipix$y))
#' sgdat <- data.frame(sgdat, patchID = apply(cbind(sgdat[, 3:2]), 1, pix.id.find, unipix))
#' weekdates <- c(40, 92)
#' potential = (1:nrow(unipix))[unipix$pop > 1000]
#' potentialpops = unipix$pop[unipix$pop > 1000]
#' potential = cbind(potential, potentialpops)
#' randSam = sample(potential[, 1], 100)
#' TreatLocs = randSam[51:100]
#' TreatLocspop = potential[potential[, 1] %in% randSam[51:100], 2]
#' ContLocs = randSam[1:50]
#' ContLocspop = potential[potential[, 1] %in% randSam[1:50], 2]
#' 
#' RCTinfo <- list(Tstart = 50,
#'                 Tend = 54,
#'                 TreatLocs = TreatLocs,
#'                 ContLocs = ContLocs,
#'                 EffectiveCoverage = 0.9)
#' 
#' denmod_RCT = DEN.spatial.RCT(weekdates, sgdat, sgdat, unipix, pixdistmat, RCTinfo, nruns = 10)
#' 
#' # inspect results
#' nruns = 10
#' DEfigs <- rep(NA, nruns)
#' for(i in 1:nruns){
#'      # attack rate untreated
#'      ARU <- (sum(denmod_RCT$Control_infections[, i]) / sum(ContLocspop))
#'      # attack rate treated
#'      ART <- (sum(denmod_RCT$Treat_infections[, i]) / sum(TreatLocspop))
#'      # calculate drug efficacy
#'      DE = 100 * (ARU - ART) / ARU
#'      DEfigs[i] = DE
#' }
#' # ! uncertainty = model uncertainty
#' # other areas of uncertainty to consider = 
#' # site randomization, cluster size and ratio, drug effective coverage, timing of trial 
#' boxplot(DEfigs, ylim = c(0, 100), main = "Drug Efficacy")
#' summary(DEfigs)

DEN.spatial.RCT <- function(weekdates,
                            fitdat,
                            pastdat,
                            unipix,
                            pixdistmat,
                            RCTinfo,
                            nruns,
                            paramsList = NULL){
  multitreatinfs = matrix(NA, nrow = length(RCTinfo$TreatLocs), ncol = nruns)
  multitreatcases = matrix(NA, nrow = length(RCTinfo$TreatLocs), ncol = nruns)
  multicontrolinfs = matrix(NA, nrow = length(RCTinfo$ContLocs), ncol = nruns)
  multicontrolcases = matrix(NA, nrow = length(RCTinfo$ContLocs), ncol = nruns)
  Pdobs_treat = matrix(NA, nrow = length(RCTinfo$TreatLocs), ncol = nruns)
  Pdobs_cont = matrix(NA, nrow = length(RCTinfo$ContLocs), ncol = nruns)
  Pdrugs = matrix(NA, nrow = length(RCTinfo$ContLocs), ncol = nruns)
  
  
  
  # nruns loop
  for(u in 1:nruns){
    ### Part 1- supply default parameters of not supplied
    #if(!exists("paramsList")){paramsList = list(NA)}
    if(missing(paramsList)){paramsList = list(NA)}
    
    # load fitted parameters
    data(finalWeights)
    
    if(!("hometime" %in% names(paramsList))){
      paramsList = c(hometime = sample(finalWeights$hometime[, 1], 1,
                                       prob = finalWeights$hometime[, 2]), 
                     paramsList)
    }
    hometime = paramsList$hometime
    
    if(!("pdetec" %in% names(paramsList))){
      #paramsList = c(pdetec = pdetec <- function(x) rbeta(x, shape1 = 2.322835, shape2 = 28.2446), 
      #               paramsList)
      paramsList = c(pdetec = sample(finalWeights$pdetec[, 1], 1,
                                     prob = finalWeights$pdetec[, 2]), 
                     paramsList)
    }
    pdetec = paramsList$pdetec
    
    if(!("mospdeath" %in% names(paramsList))){
      #paramsList = c(mospdeath = mospdeath <- function(x) rbeta(x, shape1 = 26.84987, shape2 = 107.3995),
      #               paramsList)
      paramsList = c(mospdeath = sample(finalWeights$mospdeath[, 1], 1,
                                        prob = finalWeights$mospdeath[, 2]),
                     paramsList)
    }
    mospdeath = paramsList$mospdeath
    
    if(!("stim" %in% names(paramsList))){
      paramsList = c(stim = list(stim.generate(pastdat, c(0.577, 0.659, 0.814), weekdates[1], unipix)), 
                     paramsList)
    }
    stim = paramsList$stim
    
    if(!("betaEnv" %in% names(paramsList))){
      bmean = sample(finalWeights$BetaEnv_mean[, 1], 1, prob = finalWeights$BetaEnv_mean[, 2])
      bcorrelation = sample(finalWeights$BetaEnv_cor[, 1], 1, prob = finalWeights$BetaEnv_cor[, 2])
      paramsList = c(betaEnv = list(betaEnv.generate(bmean = bmean, bcorrelation = bcorrelation, 
                                                     stim.generate(pastdat, c(0.577, 0.659, 0.814)))), 
                     paramsList)
    }
    betaEnv = paramsList$betaEnv
    
    if(!("drugtreat" %in% names(paramsList))){
      paramsList = c(drugtreat = list(data.frame(EffCoverage = 0,
                                                 Duration = 0,
                                                 Radius = 0,
                                                 Delay = 0)),
                     paramsList)
    }
    drugtreat_Eff = as.numeric(paramsList$drugtreat[1])
    drugtreat_Dur = as.numeric(paramsList$drugtreat[2])
    drugtreat_Rad = as.numeric(paramsList$drugtreat[3])
    drugtreat_Del = as.numeric(paramsList$drugtreat[4])
    
    if(!("vectreat" %in% names(paramsList))){
      paramsList = c(vectreat = list(data.frame(EffCoverage = 0,
                                                Duration = 0,
                                                Radius = 0,
                                                Delay = 0)),
                     paramsList)
    }
    vectreat_Eff = as.numeric(paramsList$vectreat[1])
    vectreat_Dur = as.numeric(paramsList$vectreat[2])
    vectreat_Rad = as.numeric(paramsList$vectreat[3])
    vectreat_Del = as.numeric(paramsList$vectreat[4])
    
    if(!("StopCriteria" %in% names(paramsList))){
      paramsList = c(StopCriteria = sum(unipix$pop),paramsList)
    }
    StopCriteria = paramsList$StopCriteria
    
    if(!("startTable" %in% names(paramsList))){
      paramsList = c(startTable = list(case.backtrack.time(fitdat, weekdates[1], unipix, stim, betaEnv, mospdeath, hometime)),
                     paramsList)
    }
    startTable = paramsList$startTable
    
    if(!("func_IIP" %in% names(paramsList))){
      paramsList = c(func_IIP = func_IIP <- function(times){plnorm(times, meanlog = log(5.9), sdlog = 0.045)},
                     paramsList)
    }
    func_IIP = paramsList$func_IIP
    
    if(!("func_EIP" %in% names(paramsList))){
      paramsList = c(func_EIP = func_EIP <- function(times){plnorm(times, meanlog = log(7), sdlog = 0.21)},
                     paramsList)
    }
    func_EIP = paramsList$func_EIP
    
    
    ### Part 2- model set up
    # pre calculate movement matrices from each of these unique patches
    mm_list <- move.matrix.gene(unipix, hometime)
    mm_list2 <- mm_list[[1]]
    mm_list_cum <- mm_list[[2]]
    
    ## Setting up human model compartments
    # humans - susceptible, infectious, Recovered, Recovered temporary (drugs)
    sing_S = unipix$pop
    sing_I = matrix(0, ncol = nrow(unipix), nrow = 15)
    sing_R = rep(0, nrow(unipix))
    sing_Rt = sing_Rt = matrix(0, nrow = 40, ncol = nrow(unipix)) # 40 day maximum effect of prophylactics
    
    # Setting up mosquito model compartments
    # mosquito compartments - susceptible, exposed, infectious
    # susceptible mosquito populations modelled as inexhaustible
    mos_S = rep(1000, nrow(unipix))
    mos_E = matrix(0, ncol = nrow(unipix), nrow = 15)
    mos_I = rep(0, nrow(unipix))
    
    # Accounting for prior immunity
    pims = round(sing_S * stim, 0)
    sing_S = sing_S - pims
    sing_R = sing_R + pims
    
    # how many exposed mosquitoes in which locations and how many days ago?
    mos_E[1, startTable$ID] = startTable$expmos
    burnin = -(min(startTable$expmostime))
    
    ### Part 2: main model run ###
    steps = burnin + 7 + 365
    trialsteps <- c(burnin + 7 + (RCTinfo$Tstart - weekdates[1])*7,
                    burnin + 7 + (RCTinfo$Tend - weekdates[1])*7)
    
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
    treatlog <- list()
    
    # add detected state
    sing_D <- matrix(0, nrow = nrow(sing_I), ncol = ncol(sing_I))
    # total population size
    tpop = sum(sing_S, sing_I, sing_R, sing_Rt)
    
    
    
    
    # main timestep loop - modifies state directly when reaching trial startpoint
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
      mosDeathProb <- rep(mospdeath, nrow(unipix))
      
      # skip reactive control
      Dspots1 = NA
      Dspots2 = NA
      mosDeathProbLocal <- rep(0, nrow(unipix))
      humDrugTreatProbLocal <- rep(0, nrow(unipix))
      
      # log treatment
      treatspots <- c(Dspots1, Dspots2)
      if(length(treatspots[!is.na(treatspots)]) == 0){treatlog[[i]] = NA} else{
        treatlog[[i]] = treatspots[!is.na(treatspots)]
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
      newsignDopts <- cbind(as.vector(sing_I2), rep(func_IIP(1:nrow(sing_I2)), ncol(sing_I2)) * pdetec)
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
      
      # breakpoint changes for start and end of trial
      if(i == (trialsteps[1] - 1)){
        newRtTrial <- rep(0, length(sing_S))
        newRtTrial[RCTinfo$TreatLocs] <- round(sing_S[RCTinfo$TreatLocs] * RCTinfo$EffectiveCoverage, 0)
        sing_S = sing_S - newRtTrial
        sing_Rt[1, ] = sing_Rt[1, ] + newRtTrial
        drugtreat_Dur = 30
      }
      
      if(i > trialsteps[2]){break}
    }
    
    # trim results to just look at the trial period
    totals = totals[trialsteps[1]:trialsteps[2], ]
    detailed_totals = detailed_totals[trialsteps[1]:trialsteps[2]]
    
    # infections and cases in treatment vs non patches over th
    start = 1
    end = trialsteps[2] - trialsteps[1]
    
    # treated patches (seroconversions)
    Tpatchstart = detailed_totals[[1]][RCTinfo$TreatLocs, 1] + detailed_totals[[1]][RCTinfo$TreatLocs, 4]
    Tpatchend = detailed_totals[[end]][RCTinfo$TreatLocs, 1] + detailed_totals[[end]][RCTinfo$TreatLocs, 4]
    
    Treatseroconv = Tpatchstart - Tpatchend
    
    # control patches (seroconversions)
    Cpatchstart = detailed_totals[[1]][RCTinfo$ContLocs, 1]
    Cpatchend = detailed_totals[[end]][RCTinfo$ContLocs, 1]
    
    Controlseroconv = Cpatchstart - Cpatchend
    
    # detected new cases
    Tpatchcases = rep(0, length(detailed_totals[[1]][RCTinfo$TreatLocs, 1]))
    Cpatchcases = rep(0, length(detailed_totals[[1]][RCTinfo$ContLocs, 1]))
    for(l in 1:(end - start)){
      Tpatchcases = Tpatchcases + detailed_totals[[l]][RCTinfo$TreatLocs, 6]
      Cpatchcases = Cpatchcases + detailed_totals[[l]][RCTinfo$ContLocs, 6]
    }
    
    # return results (first level)
    
    multitreatinfs[, u] = Treatseroconv
    multitreatcases[, u] = Tpatchcases
    multicontrolinfs[, u] = Controlseroconv
    multicontrolcases[, u] = Cpatchcases
    # only observing IgG seronegatives (and those given drugs)
    Pdobs_cont[, u] = sum(detailed_totals[[1]][RCTinfo$ContLocs, c(1, 2, 4)]) * (RCTinfo$Tend - RCTinfo$Tstart) * 7 
    Pdobs_treat[, u] = sum(detailed_totals[[1]][RCTinfo$TreatLocs, c(1, 2, 4)]) * (RCTinfo$Tend - RCTinfo$Tstart) * 7
    Pdrugs[, u] = max(totals$Rt)
      
  }
  
  # return results
  return(list(Treat_infections = multitreatinfs,
              Treat_cases = multitreatcases,
              Control_infections = multicontrolinfs,
              Control_cases = multicontrolcases,
              Person_days_observation_treat = Pdobs_treat,
              Person_days_observation_cont = Pdobs_cont,
              People_treated_drugs = Pdrugs))
}

