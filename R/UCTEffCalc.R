#' Summarises the results from multiple simulations of a UCT
#'
#' Calculates Effectiveness and ttoal people treated from the outputs of multiple simulated trials
#' @param denmod_UCT A list of multiple SpatialDengue model objects
#' @param unipix Universal pixel lookup table, see ?make.unipix
#' @keywords summary
#' @details plots a boxplot of Effectiveness against dengue infection and dengue cases and returns a data frame of the Effectiveness and 
#' total person-days of observation in the trial and total people treated with drugs for each realisation of the trial (model run)
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
#' denmod_UCT = DEN.spatial.UCT(weekdates, sgdat, unipix, pixdistmat, UCTinfo, nruns = 10)
#' effectiveness <- UCT.Eff.calc(denmod_UCT, unipix)



UCT.Eff.calc <- function(denmod_UCT, unipix){
  nruns = length(denmod_UCT)
  Efftab = data.frame(Eff_infs = rep(NA, nruns),
                      Eff_cases = rep(NA, nruns))
  trialdeets <- data.frame(Patches_treated = rep(NA, nruns),
                           Patches_control = rep(NA, nruns),
                           People_days_observed = rep(NA, nruns),
                           People_treated = rep(NA, nruns))
  
  for(i in 1:nruns){
    # table out treated and control patches and when they were recruited
    Ptreat = data.frame(daystart = unlist(apply(cbind(1:365, unlist(lapply(denmod_UCT[[i]]$treatlog, length))),
                                                1, function(x) rep(x[1], x[2]))),
                        patch = unlist(denmod_UCT[[i]]$treatlog))
    Ptreat = Ptreat[!is.na(Ptreat[, 2]), ]
    
    Pcont = data.frame(daystart = unlist(apply(cbind(1:365, unlist(lapply(denmod_UCT[[i]]$contlog, length))),
                                               1, function(x) rep(x[1], x[2]))),
                       patch = unlist(denmod_UCT[[i]]$contlog))
    Pcont = Pcont[!is.na(Pcont[, 2]), ]
    
    # if no trial because outbreak was over, just throw NAs- trial never would ahve started
    Ptreatinfs = rep(NA, nrow(Ptreat))
    PtreatNinfs = rep(NA, nrow(Ptreat))
    Ptreatcases = rep(NA, nrow(Ptreat))
    Peopletreated = rep(NA, nrow(Ptreat))
    PTreatdObs = rep(NA, nrow(Ptreat))
    if(nrow(Ptreat) > 0){
      # Q1 efficacy for infections and cases
      # first cases and infections in treatment patches
      for(k in 1:nrow(Ptreat)){
        enroll = Ptreat[k, 1]
        fup = Ptreat[k, 1] + UCTinfo$Tfollowup
        
        # extract goo
        newcases = 0
        for(h in enroll:fup){
          newcases = newcases + denmod_UCT[[i]]$Patchtotals[[h]][, 6][Ptreat[k, 2]]
        }
        Ptreatcases[k] = newcases
        
        enroll_C <- denmod_UCT[[i]]$Patchtotals[[enroll]][, c(1, 2, 4, 6)][Ptreat[k, 2], ]
        fup_C <- denmod_UCT[[i]]$Patchtotals[[fup]][, c(1, 2, 4, 6)][Ptreat[k, 2], ]
        
        # seroconverted, not seroconverted
        Ptreatinfs[k] <- as.numeric(sum(enroll_C[1:3]) - sum(fup_C[1:3]))
        PtreatNinfs[k] <- as.numeric(fup_C[1] + fup_C[3])
        
        # people treated
        Peopletreated[k] = as.numeric(enroll_C[3])
        
        # people-days under observation
        PTreatdObs[k] = (fup - enroll) * sum(enroll_C[1:3])
      }
    }
    
    
    
    
    # then cases and infections in control patches
    Pcontinfs = rep(NA, nrow(Pcont))
    PcontNinfs = rep(NA, nrow(Pcont))
    Pcontcases = rep(NA, nrow(Pcont))
    PcontdObs = rep(NA, nrow(Ptreat))
    if(nrow(Pcont) > 0){
      for(k in 1:nrow(Pcont)){
        enroll = Pcont[k, 1]
        fup = Pcont[k, 1] + UCTinfo$Tfollowup
        
        # extract goo
        newcases = 0
        for(h in enroll:fup){
          newcases = newcases + denmod_UCT[[i]]$Patchtotals[[h]][, 6][Pcont[k, 2]]
        }
        Pcontcases[k] = newcases
        
        enroll_C <- denmod_UCT[[i]]$Patchtotals[[enroll]][, c(1, 2, 4, 6)][Pcont[k, 2], ]
        fup_C <- denmod_UCT[[i]]$Patchtotals[[fup]][, c(1, 2, 4, 6)][Pcont[k, 2], ]
        
        # seroconverted, not seroconverted
        Pcontinfs[k] <- as.numeric(sum(enroll_C[1:3]) - sum(fup_C[1:3]))
        PcontNinfs[k] <- as.numeric(fup_C[1] + fup_C[3])
        
        # people-days under observation
        PcontdObs[k] = (fup - enroll) * sum(enroll_C[1:3])
      }
    }
    
    # now calculate effectiveness
    
    # first gather population data
    Ptreat$pop = unipix[match(Ptreat$patch, unipix[, 1]), 3]
    Pcont$pop = unipix[match(Pcont$patch, unipix[, 1]), 3]
    
    # attack rate untreated clusters
    ARU_inf = sum(Pcontinfs) / sum(Pcont$pop)
    ARU_cases = sum(Pcontcases) / sum(Pcont$pop)
    
    # attack rate treated clusters
    ART_inf = sum(Ptreatinfs) / sum(Ptreat$pop)
    ART_cases = sum(Ptreatcases) / sum(Ptreat$pop)
    
    # efficacy
    Efftab[i, 1] = 100 * (ARU_inf - ART_inf) / ARU_inf
    Efftab[i, 2] = 100 * (ARU_cases - ART_cases) / ARU_cases
    
    # trial details
    trialdeets$Patches_treated[i] = nrow(Ptreat)
    trialdeets$Patches_control[i] = nrow(Pcont)
    trialdeets$People_treated[i] = sum(Peopletreated)
    trialdeets$People_days_observed[i] = sum(PTreatdObs, PcontdObs)
  }
  boxplot(Efftab, ylim = c(0, 100), main = "Drug Efficacy", names = c("Against infection", "Against cases"))
  #print(summary(Efftab))
  # how many people given Px?
  #print(summary(trialdeets$People_treated))
  
  # return dataframe
  rtndf <- data.frame(Inf_eff = Efftab[, 1],
                      Case_eff = Efftab[, 2],
                      People_days_observed = trialdeets$People_days_observed,
                      People_treated = trialdeets$People_treated)
  return(rtndf)
}
