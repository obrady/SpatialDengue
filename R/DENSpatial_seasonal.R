#' Run a single realisation of the spatial dengue model
#'
#' Main wrapper function for running a single realisation of the dengue spatial model, may take a while to run depending on specification
#' @param weekdates Two element vector of the start and end weeks of the simulation over which the model will be evaluated over
#' @param fitdat Data frame of the locations, numbers and timings (in weeks) of cases to fit the model to, see ?sgdat
#' @param pastdat Data frame of the locations, numbers and timings (in weeks) of all cases in the dataset (is used to generate the starting immunity profile), see ?sgdat
#' @param unipix Universal pixel lookup table, see ?make.unipix
#' @param pixdistmat A patch distance matrix, see example
#' @param steprun integer, number of days for which the simulation should run, excluding burn in period
#' @param paramsList Optional parameter list. If not supplied returns to defaults, see tutorial for full parameter list, see ?model.run for full list and explanation of parameters
#' @details This function undertakes three main processes: i) fills in parameters with default options if not supplied, ii) generates a human movement matrix between patches,
#' ii) runs the model simulation.
#' within "paramsList()". Missing parameters will return to their default values. This function returns a outputs in a three element list giving:
#' i) daily counts of each model state (S, I, R, etc) sumed across the whole landscape, ii) daily counts of each model state for every patch,
#' iii) A treatment log (of length = number of steps) detailing which patches were treated each day
#' @keywords model simulation
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
#' denmod_sim <- DEN.spatial(weekdates, sgdat, sgdat, unipix, pixdistmat, 365)
#' # model run with reactive drug deployment
#' drugtreat = data.frame(EffCoverage = 0.9, Duration = 30, Radius = 1000, Delay = 0)
#' denmod_sim_drugs <- DEN.spatial(weekdates, sgdat, unipix, pixdistmat, 365, paramsList = list(drugtreat = drugtreat))



DEN.spatial.season <- function(weekdates,
                        fitdat,
                        pastdat,
                        sgpop,
                        unipix,
                        pixdistmat,
                        steprun,
                        seasonal_vector,
                        seasonal_start,
                        paramsList = NULL){
  ### Part 1- supply default parameters of not supplied
  #if(!exists("paramsList")){paramsList = list(NA)}
  if(missing(paramsList)){paramsList = list(NA)}

  # load fitted parameters
  data(finalWeights)
  # presample parameters
  # if some parameters available, then adjust weighting to preferentially select complementary parameter values
  pamnames = c("hometime", "pdetec", "mospdeath", "betaEnv", "Mov_model_type")
  apams <- pamnames %in% names(paramsList)

  if(any(apams) & (sum(apams) != length(apams))){
    Aweights = matrix(0, nrow = nrow(finalWeights), ncol = length(apams))
    for(i in 1:length(apams)){
      if(apams[i]){
        fname = pamnames[i]
        # weighting proportional to absolute deviation
        if(fname == "betaEnv"){
          val1 = sqrt((finalWeights$betaEnv_mean - mean(paramsList$betaEnv))^2)
          val1 = val1 / max(val1)
          #val2 = sqrt((finalWeights$betaEnv_var - var(paramsList$betaEnv))^2)
          #val2 = val2 / max(val2)

          Aweights = cbind(Aweights, val1)
          #Aweights = cbind(Aweights, val2)
        }else{
          Aweights[, i] = unlist(1 - sqrt((finalWeights[names(finalWeights) == fname] - paramsList[names(paramsList) == pamnames[i]])^2))
        }
      }
    }
    Aweights = Aweights[, colSums(Aweights) != 0]
    if(is.vector(Aweights)){Aweights = matrix(Aweights, ncol = 1)}
    Aweights = apply(Aweights, 1, mean)
    # final max
    fmax = max(finalWeights$totaldevs)
    # scale A weights to be on the same scale as totaldevs
    Aweights = max(finalWeights$totaldevs) * (Aweights / max(Aweights))
    # add finalweights and aweights together then scale back to original minimum and maximum
    sample_probs = finalWeights$totaldevs + Aweights
    sample_probs = sample_probs - min(sample_probs)
    sample_probs = fmax * (sample_probs / max(sample_probs))

    sampams = finalWeights[sample(1:nrow(finalWeights), 1, prob = sample_probs), ]
  }else{
    sampams <- finalWeights[sample(1:nrow(finalWeights), 1, prob = finalWeights$totaldevs), ]
  }


  # now assign parameters back to the list
  if(!("hometime" %in% names(paramsList))){
    paramsList = c(hometime = sampams$hometime, paramsList)
    }
  hometime = paramsList$hometime

  if(!("pdetec" %in% names(paramsList))){
    #paramsList = c(pdetec = pdetec <- function(x) rbeta(x, shape1 = 2.322835, shape2 = 28.2446),
    #               paramsList)
    paramsList = c(pdetec = sampams$pdetec, paramsList)
    }
  pdetec = paramsList$pdetec

  if(!("mospdeath" %in% names(paramsList))){
    #paramsList = c(mospdeath = mospdeath <- function(x) rbeta(x, shape1 = 26.84987, shape2 = 107.3995),
    #               paramsList)
    paramsList = c(mospdeath = sampams$mospdeath, paramsList)
  }
  mospdeath = paramsList$mospdeath

  if(!("stim" %in% names(paramsList))){
    paramsList = c(stim = list(stim.generate(pastdat, c(0.577, 0.659, 0.814), weekdates[1], sgpop, unipix)),
                   paramsList)
  }
  stim = paramsList$stim

  if(!("betaEnv" %in% names(paramsList))){
    bmean_1 = sampams$betaEnv_mean
    bmean = bmean_1 * seasonal_vector[seasonal_start]
    bcorrelation = sampams$betaEnv_cor
    paramsList = c(betaEnv = list(betaEnv.generate(bmean = bmean, bcorrelation = bcorrelation,
                                                   stim = paramsList$stim,
                                                   unipix)),
                   paramsList)
  }
  betaEnv = paramsList$betaEnv

  if(!("drugtreat" %in% names(paramsList))){
    paramsList = c(drugtreat = list(data.frame(EffCoverage = 0,
                                               Duration = 0,
                                               Radius = 0,
                                               Delay = NA)),
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
                                              Delay = NA)),
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

  # !!!no longer used
  #if(!("startTable" %in% names(paramsList))){
  #  paramsList = c(startTable = list(case.backtrack.time(casedat, weekdates[1], unipix, stim, betaEnv, mospdeath, hometime)),
  #                 paramsList)
  #}
  #startTable = paramsList$startTable

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

  if(!("Mov_model_type" %in% names(paramsList))){
    paramsList = c(Mov_model_type = c("exponential", "gravity", "radiation")[sampams$Mov_model_type],
                   paramsList)
  }
  Mov_model_type = paramsList$Mov_model_type


  ### Part 2- model set up
  # pre calculate movement matrices from each of these unique patches
  mm_list <- move.matrix.gene(unipix, hometime, Mov_model_type)
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
  # susceptible mosquito populations modelled as inexhaustible and directly proportional to human population
  mos_S = unipix$pop
  mos_E = matrix(0, ncol = nrow(unipix), nrow = 15)
  mos_I = rep(0, nrow(unipix))

  # Accounting for prior immunity
  pims = round(sing_S * stim, 0)
  sing_S = sing_S - pims
  sing_R = sing_R + pims

  ## Initial conditions:
  # 01 calculate the total numbers of various different states at the begining
  start_sing = data.frame(Number = fitdat$Number_of_cases[fitdat$Week == weekdates[1]],
                            Location = fitdat$patchID[fitdat$Week == weekdates[1]])
  start_sing = aggregate(start_sing$Number, list(start_sing$Location), sum)
  start_sing = data.frame(Number_I = start_sing$x,
                          Number_trueI = round(start_sing$x / (8 * pdetec), 0) - start_sing$x, # need to times pdetec by 8 to go from daily value to value over all symptomatic days
                          Location = start_sing$Group.1,
                          Timing = sample(1:15, nrow(start_sing), replace = T),
                          Pop = unipix$pop[match(start_sing$Group.1, unipix$patchID)])
  # if probability of detection (over course of illness) is greater than 1 (i.e. pdetec > 1/8) set extra true cases to 0
  if(nrow(start_sing) > 0){
    if(any(start_sing$Number_trueI < 0)){start_sing$Number_trueI[start_sing$Number_trueI < 0] = 0}
  }else{
    if(start_sing$Number_trueI < 0){start_sing$Number_trueI = 0}
  }

  # add mosquito compartments
  # how many mosquitoes would need to have been infected to give this number of cases?
  # 1) calculate Rt
  VC <- (1 * (5.457631 / 8) * mean(betaEnv) ^ 2 *  (1 - mospdeath) ^ 7) / (-log((1 - mospdeath)))
  Rt <- 8 * VC * (1 - mean(stim))
  # 2) calculate the number of humans who must have been infected 1 transmission generation ago
  Inf_1g <- (start_sing$Number_I + start_sing$Number_trueI) / Rt

  # number of exposed mosquitoes over their 5.45 effective infectious days
  start_sing$Number_mosE <- round(5.457631 * Inf_1g * mean(betaEnv), 0)
  # and how many of these have survived EIP?
  start_sing$Number_mosI <- round(start_sing$Number_mosE * (1 - mospdeath) ^ 7, 0)

  # 02 allocate these in space and time
  for(i in 1:nrow(start_sing)){
    # frist allocate the reported cases
    fmove = mm_list2[[start_sing$Location[i]]]
    # allocate reported Infectious cases first
    sing_I[cbind(start_sing$Timing[i], start_sing$Location[i])] = sing_I[cbind(start_sing$Timing[i], start_sing$Location[i])] + start_sing$Number_I[i]
    # then unreported infectious cases
    locsam <- sample(1:nrow(unipix), sum(start_sing$Number_trueI[i]), prob = fmove, replace = T)
    timsam <- sample(1:15, sum(start_sing$Number_trueI[i]), replace = T)
    for(k in 1:length(locsam)){
      sing_I[cbind(timsam[k], locsam[k])] = sing_I[cbind(timsam[k], locsam[k])] + 1
    }

    # then infectious mosquitoes
    locsam <- sample(1:nrow(unipix), sum(start_sing$Number_mosI[i]), prob = fmove, replace = T)
    for(k in 1:length(locsam)){
      mos_I[locsam[k]] = mos_I[locsam[k]] + 1
    }

    # then exposed mosquitoes
    locsam <- sample(1:nrow(unipix), sum(start_sing$Number_mosE[i]), prob = fmove, replace = T)
    timsam <- sample(1:15, sum(start_sing$Number_mosE[i]), replace = T)
    for(k in 1:length(locsam)){
      mos_E[cbind(timsam[k], locsam[k])] = mos_E[cbind(timsam[k], locsam[k])] + 1
    }
  }



  ### Part 2: main model run ###
  steps = steprun

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
  treatlog_drug2 <- list()
  treatlog_vec2 <- list()
  

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
    

    #calculate new beta at each stage
    
    if ((i + seasonal_start) <= 365){
      t = i
    bmean = bmean_1 * seasonal_vector[seasonal_start + t]
   betaEnv = betaEnv.generate(bmean = bmean, bcorrelation = bcorrelation,
                                                   stim = stim,
                                                   unipix)
    
    }
    
    if ((i + seasonal_start) > 365) {
      t3 = i + seasonal_start
      t2 = trunc(t3 /365, 0)
      t = t3 - (t2*365)
      bmean = bmean_1 * seasonal_vector[t]
      betaEnv = betaEnv.generate(bmean = bmean, bcorrelation = bcorrelation,
                                 stim = stim,
                                 unipix = unipix)}
  
    
    
    
    # 1A- P(mosquito infected)
    mosInfProb <- mos.inf.prob(sing_I, betaEnv, unipix, mm_list2, mm_list_cum)
    # 1B- P(mosquito completes EIP)
    # func_EIP(days)
    # 1C- P(mosqutio dies naturally)
    #mosDeathProb <- mospdeath(nrow(unipix))
    mosDeathProb <- rep(mospdeath, nrow(unipix))
    # 1D- calculate drug targetted patches and P(mosquto dies from RVC) and P(human treated with drugs)
    # don't bother doing if Effective coverage = 0
    if(sum(vectreat_Eff, drugtreat_Eff) != 0){
      # identify detected locations
      if(i == 1){newsingDs = matrix(0, nrow = nrow(sing_D), ncol = ncol(sing_D))}
      Dspots1 <- find.Dspots(newsingDs, pixdistmat, radius = drugtreat_Rad)
      Dspots2 <- find.Dspots(newsingDs, pixdistmat, radius = vectreat_Rad)
      mosDeathProbLocal <- rep(0, nrow(unipix))
      humDrugTreatProbLocal <- rep(0, nrow(unipix))
      treatlog_drug2[[i]] = NA
      treatlog_vec2[[i]] = NA
      
      # Drugs
      if(sum(!is.na(Dspots1)) > 0){
        # assign to treatment log
        treatlog_drug[[i]] = Dspots1
       
        

        # assign drugs (+/- delay between patient detection and local area visit)
        if((sum(!is.na(drugtreat_Del))> 0)  & i > drugtreat_Del){ # can only start if drugtreat_Del- days in
          scheduleday = i - drugtreat_Del
          if(!all(is.na(treatlog_drug[[scheduleday]]))){
            humDrugTreatProbLocal[treatlog_drug[[scheduleday]]] <- drugtreat_Eff
            treatlog_drug2[[i]] = treatlog_drug[[scheduleday]]
          }
        }
       
      }else{treatlog_drug[[i]] = NA
            }

      # Vector control
      if(sum(!is.na(Dspots2)) > 0){
        # assign to treatment log
        treatlog_vec[[i]] = Dspots2
        
        

        # assign drugs (+/- delay between patient detection and local area visit)
        if( (sum(!is.na(vectreat_Del))> 0) & i > vectreat_Del){ # can only start if vectreat_Del- days in
          scheduleday = i - vectreat_Del
          if(!all(is.na(treatlog_vec[[scheduleday]]))){
            mosDeathProbLocal[treatlog_vec[[scheduleday]]] <- vectreat_Eff
            treatlog_vec2[[i]] = treatlog_vec[[scheduleday]]
          }
        }
        
      }else{treatlog_vec[[i]] = NA
      }
    }else{
      Dspots1 = NA
      Dspots2 = NA
      mosDeathProbLocal <- rep(0, nrow(unipix))
      humDrugTreatProbLocal <- rep(0, nrow(unipix))
      treatlog_drug[[i]] = NA
      treatlog_vec[[i]] = NA
      treatlog_drug2[[i]] = NA
      treatlog_vec2[[i]] = NA
      
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
    #newsignDopts <- cbind(as.vector(sing_I2), rep(func_IIP(1:nrow(sing_I2)), ncol(sing_I2)) * pdetec(1))
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
  }

  # compile treatlog
  treatlog = list(drug1 = treatlog_drug,
                  drug2 = treatlog_drug2,
                  vec1 = treatlog_vec,
                  vec2 = treatlog_vec2)

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
  mainmodrun <- templist


  # if StopCriteria are exceeded, model.run will return NA (useful for fitting),
  # otherwise, Part 4:results:

  if(length(mainmodrun) > 1){

    # Return results
    model_outlist <- list()
    model_outlist[[1]] <- mainmodrun$totals
    model_outlist[[2]] <- mainmodrun$detailed_totals
    model_outlist[[3]] = mainmodrun$treatlog
    names(model_outlist) = c("Poptotals",
                             "Patchtotals",
                             "treatlog")

    return(model_outlist)
  }else{return(NA)}
}
