#' Multinomial state transition sampler
#'
#' Ensures arrivals to or departures from a model state have only one source / destination
#' @param x a multi-dimentional list with level 1 giving the input and output processes and each element of level 1 containign a two element list with the starting state, e.g. sing_S and the probability of leaving that state
#' @param type a character indicating for which state ins and outs should be calculated for, e.g. "mos_E"
#' @details Returns a list of new inputs and outputs from the specified state, maintaining the structure of each state.
#' @export
#' @examples
#' See tutorial




multinom <- function(x, type){
  # mos mos_E in and out
  if(type == "mos_E"){
    # Infected mosquitoes
    newmosE = apply(cbind(x[[1]][[1]], x[[1]][[2]]), 1, function(u) rbinom(1, u[1], u[2]))
                    
    for (i in 1:length(newmosE)){
      newmosE[[i]] = ifelse(is.na(newmosE[[i]]), 0, newmosE[[i]]) }           
                    
    # order randomization of out compartments
    ord <- sample(1:3, 3)
    ord2 <- data.frame(name = c("EIP", "natdeath", "cdeath"), orignum = 1:3, newnum = ord)
    
    mos_Es = x[[2]][[1]]
    for(k in 1:3){
      if(ord2[k, 3] == 1){newmosI <- EIP.MN.func(mos_Es, x[[2]][[2]]); mos_Es = mos_Es - newmosI}
      if(ord2[k, 3] == 2){newmosD <- natdeath.MN.func(mos_Es, x[[3]][[2]]); mos_Es = mos_Es - newmosD}
      if(ord2[k, 3] == 3){newmosD2 <- Cdeath.MN.func(mos_Es, x[[4]][[2]]); mos_Es = mos_Es - newmosD2}
    }
    rtnlist  = list(newmosE, newmosI, newmosD, newmosD2)
    names(rtnlist) = c("newmosE", "newmosI", "newmosD", "newmosD2")
    return(rtnlist)
  }
  
  # for mos_I out
  if(type == "mos_I"){
    # order randomization of out compartments
    ord <- sample(1:2, 2)
    ord2 <- data.frame(name = c("natdeath", "cdeath"), orignum = 1:2, newnum = ord)
    
    mos_Is = matrix(x[[1]][[1]], nrow = 1)
    for(k in 1:2){
      if(ord2[k, 3] == 1){newmosD <- natdeath.MN.func(mos_Is, x[[1]][[2]]); mos_Is = mos_Is - newmosD}
      if(ord2[k, 3] == 2){newmosD2 <- Cdeath.MN.func(mos_Is, x[[2]][[2]]); mos_Is = mos_Is - newmosD2}
    }
    rtnlist  = list(as.vector(newmosD), as.vector(newmosD2))
    names(rtnlist) = c("newmosD", "newmosD2")
    return(rtnlist)
  }
  
  # for sing_S in and out
  if(type == "sing_S"){
    # waning prophylactic drugs (not randomised)
    newhumS = matrix(0, ncol = ncol(x[[3]][[1]]), nrow = nrow(x[[3]][[1]]))
    if(sum(x[[3]][[1]]) > 0){
      for(k in 1:nrow(x[[3]][[1]])){
        humageA = cbind(x[[3]][[1]][k, ], x[[3]][[2]](k))
        humageAind = humageA[, 1] > 0
        if(sum(humageAind) > 0){
          # required when just one patch matched
          if(length(humageA[humageAind, ]) == 2){mathumage = matrix(humageA[humageAind, ], nrow = 1)}else{(mathumage = humageA[humageAind, ])}
          newhumS[k, ][ humageAind] = apply(mathumage, 1, function(u) rbinom(1, u[1], u[2]))
        }
      }
    }
    
    # order randomization of outcomponents
    ord <- sample(1:2, 2)
    ord2 <- data.frame(name = c("humInf", "humDrugtreat"), orignum = 1:2, newnum = ord)
    sing_Ss = x[[1]][[1]]
    # if drugs are not in play, omit the sampling- quicker
    if(sum(x[[2]][[2]]) > 0){
      for(k in 1:2){
        if(ord2[k, 3] == 1){newhumI <- apply(cbind(sing_Ss, x[[1]][[2]]), 1, function(u) rbinom(1, u[1], u[2])); 
           for (i in 1:length(newhumI)){
          newhumI[[i]] = ifelse(is.na(newhumI[[i]]), 0, newhumI[[i]])};
            sing_Ss = sing_Ss - newhumI}
        if(ord2[k, 3] == 2){newhumRt <- apply(cbind(sing_Ss, x[[2]][[2]]), 1, function(u) rbinom(1, u[1], u[2])); 
           for (i in 1:length(newhumI)){
          newhumI[[i]] = ifelse(is.na(newhumI[[i]]), 0, newhumI[[i]])};
                                              sing_Ss = sing_Ss - newhumRt}
      }
    }else{newhumI <- apply(cbind(x[[1]][[1]], x[[1]][[2]]), 1, function(u) rbinom(1, u[1], u[2]));
           for (i in 1:length(newhumI)){
          newhumI[[i]] = ifelse(is.na(newhumI[[i]]), 0, newhumI[[i]])};
                           newhumRt = rep(0, length(x[[1]][[2]]))}
    
    rtnlist  = list(newhumS, newhumI, newhumRt)
    names(rtnlist) = c("newhumS", "newhumI", "newhumRt")
    return(rtnlist)
  }
}
