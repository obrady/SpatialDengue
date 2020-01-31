#' Simulate an urban cluster trial of dengue prophylactic drugs
#'
#' Provide basic information on the UCT and the number of stochastic runs over which model runs should be averaged
#' @param weekdates Two element vector of the start and end weeks of the simulation over which the mdoel will be evaluated over
#' @param fitdat Data frame of the locations, numbers and timings (in weeks) of cases to fit the model to, see ?sgdat
#' @param pastdat Data frame of the locations, numbers and timings (in weeks) of all cases in the dataset (is used to generate the starting immunity profile), see ?sgdat
#' @param unipix Universal pixel lookup table, see ?make.unipix
#' @param pixdistmat A patch distance matrix, see example
#' @param UCTinfo A data frame including details on Urban Cluster Trial tart time, delay between enrollment and followup, number of treatment and control clustes and drug effective coverage
#' @param nruns integer, number of stochastic runs of the model over which results should be averaged
#' @param paramsList Optional parameter list. If not supplied returns to defaults, see tutorial for full parameter list, see ?model.run for full list and explanation of parameters
#' @details runs "nruns" number of stochastic simulations and enrolls treatment or control clusters once the trail has begun then stops enrolling once the target
#' number of clusters has been reached. Returns the normal DEN.spatial sumamry results but with "treatlog" and "contlog" which detail which patches were enrolled into treatment or control arms at which point
#' @keywords model simulation UCT
#' @export
#' @examples
#' data(sgdat)
#' data(sgpop)
#' sgpop <- pop.process(sgpop, agg = 10)
#' unipix <- make.unipix(sgpop)
#' pixdistmat <- distm(cbind(unipix$x, unipix$y))
#' sgdat <- data.frame(sgdat, patchID = apply(cbind(sgdat[, 3:2]), 1, pix.id.find, unipix))
#' weekdates <- c(40, 92)
#' 
#' UCTinfo <- list(Tstart = 50,
#'                Tfollowup = 28,
#'                ntreat = 30,
#'                ncontrol = 30,
#'                DrugEfficacy = 0.9)
#' 
#' denmod_UCT = DEN.spatial.UCT(weekdates, sgdat, sgdat, unipix, pixdistmat, UCTinfo, nruns = 10)
#' effectiveness <- UCT.Eff.calc(denmod_UCT, unipix)





DEN.spatial.UCT<- function(weekdates,
                           fitdat,
                           pastdat,
                           unipix,
                           pixdistmat,
                           UCTinfo,
                           nruns,
                           paramsList = NULL){
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
                                                   stim.generate(pastdat, c(0.577, 0.659, 0.814), weekdates[1], unipix), 
                                                   unipix)), 
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
  
  
  # collector
  model_outlist2 = list()
  
  for(y in 1:nruns){
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
    steps = burnin + 7 + 365
    # normal model run but with some adjustments for using UCTinfo
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
    treatlog <- as.list(rep(NA, steps))
    contlog <- as.list(rep(NA, steps))
    
    # add detected state
    sing_D <- matrix(0, nrow = nrow(sing_I), ncol = ncol(sing_I))
    # total population size
    tpop = sum(sing_S, sing_I, sing_R, sing_Rt)
    
    # main timestep loop
    for(i in 1:steps){
      # SIR-SEI formulation
      # Part 1 calculates stage transition probabilities
      # Part 2 updates states
      
      # add counter for patches enrolled in Urban Cluster Trial
      UCTpatchnum = 0
      
      # Part 1 stage transition probabilities
      
      # 1A- P(mosquito infected)
      mosInfProb <- mos.inf.prob(sing_I, betaEnv, unipix, mm_list2, mm_list_cum)
      # 1B- P(mosquito completes EIP)
      # func_EIP(days)
      # 1C- P(mosqutio dies naturally)
      mosDeathProb <- rep(mospdeath, nrow(unipix))
      # 1D- calculate drug targetted patches and P(mosquto dies from RVC) and P(human treated with drugs)
      # don't bother doing if Effective coverage = 0
      if((sum(vectreat_Eff, UCTinfo$DrugEfficacy) != 0)  & 
         (i >= ((UCTinfo$Tstart - weekdates[1]) * 7 + burnin + 7)) & 
         (UCTpatchnum < sum(UCTinfo$ntreat, UCTinfo$ncontrol))){
        # identify detected locations
        if(i == 1){newsingDs = matrix(0, nrow = nrow(sing_D), ncol = ncol(sing_D))}
        Dspots1 <- find.Dspots(newsingDs, pixdistmat, radius = 0)
        Dspots2 <- find.Dspots(newsingDs, pixdistmat, radius = vectreat_Rad)
        
        # distill to just positive spots
        if(!is.na(Dspots1[1])){Dspots1 = (1:nrow(unipix))[apply(as.matrix(Dspots1), 2, sum) > 0]}
        if(!is.na(Dspots2[1])){Dspots2 = (1:nrow(unipix))[apply(as.matrix(Dspots2), 2, sum) > 0]}
        
        # remove them if they are already part of the trial
        if(!is.na(Dspots1[1])){Dspots1 = Dspots1[!(Dspots1 %in% c(unlist(treatlog), unlist(contlog)))]}
        
        if((length(Dspots1) > 0) & !is.na(Dspots1[1])){
          # asign to treatment or control clusters
          Treatprob <- UCTinfo$ntreat / (UCTinfo$ntreat + UCTinfo$ncontrol)
          clusD <- sample(c("Treat", "Control"), length(Dspots1), 
                          prob = c(Treatprob, 1 - Treatprob), replace = T)
          # balance to ensure total cluster number for treatment and control is met
          treatexisting <- unique(unlist(treatlog))
          contexisting <- unique(unlist(contlog))
          treatstillfree <- UCTinfo$ntreat - length(treatexisting[!is.na(treatexisting)])
          contstillfree <- UCTinfo$ncontrol - length(contexisting[!is.na(contexisting)])
          # trim clusD to balance numbers
          if(sum(clusD == "Treat") > treatstillfree){DspotsT1 = sample(Dspots1[clusD == "Treat"], treatstillfree)}else{DspotsT1 = Dspots1[clusD == "Treat"]}
          if(sum(clusD == "Control") > contstillfree){DspotsT2 = sample(Dspots1[clusD == "Control"], contstillfree)}else{DspotsT2 = Dspots1[clusD == "Control"]}
          Dspots1 = c(DspotsT1, DspotsT2)
          
          # initialise
          humDrugTreatProbLocal <- rep(0, nrow(unipix))
          mosDeathProbLocal <- rep(0, nrow(unipix))
          
          # treatment
          treatlog[[i]] = c(treatlog[[i]], DspotsT1)
          humDrugTreatProbLocal[DspotsT1] <- UCTinfo$DrugEfficacy
          
          # control
          contlog[[i]] = c(contlog[[i]], DspotsT2)
          
          # updating totals
          UCTpatchnum = UCTpatchnum + sum(DspotsT1, DspotsT2)
        }
      }else{
        Dspots1 = NA
        Dspots2 = NA
        mosDeathProbLocal <- rep(0, nrow(unipix))
        humDrugTreatProbLocal <- rep(0, nrow(unipix))
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
                                 list(sing_Rt,  function(times) pnorm(times, UCTinfo$Tfollowup, 0))),
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
    
    # trim results to remove the burn-in period
    totals = totals[(burnin + 8):nrow(totals), ]
    detailed_totals = detailed_totals[(burnin + 8):length(detailed_totals)]
    treatlog = treatlog[(burnin + 8):length(treatlog)]
    contlog = contlog[(burnin + 8):length(contlog)]
    
    # collate final results
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
                     treatlog,
                     contlog)
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
                        "treatlog",
                        "contlog")
    mainmodrun <- templist
    
    # Return results
    model_outlist <- list()
    model_outlist[[1]] <- mainmodrun$totals
    model_outlist[[2]] <- mainmodrun$detailed_totals
    model_outlist[[3]] = mainmodrun$treatlog
    model_outlist[[4]] = mainmodrun$contlog
    names(model_outlist) = c("Poptotals",
                             "Patchtotals",
                             "treatlog",
                             "contlog")
    
    model_outlist2[[y]] = model_outlist
  }
 return(model_outlist2) 
}