#' Creating a movement matrix
#'
#' Applys a radiation model to a population surface to calcualte relative flux of human movement between patches
#' @param unipix Universal pixel lookup table, see ?make.unipix
#' @param hometime Paramter that controls the spatial concentration of movement
#' @param type Human movement model type. Choose from  "radiation", "gravity" or "exponential"
#' @details Returns a 2 element list: i) a list of length equal to the number of patches. Each item in this list is a movement matrix (in vector form) with fluxes between the home patch and all other patches,
#' ii) a movement matrix (in vector form) that sums total time allocation in each pixel (from residents and visitors)
#' @keywords movement
#' @export
#' @examples
#' data(sgpop)
#' sgpop <- pop.process(sgpop, agg = 10)
#' unipix <- make.unipix(sgpop)
#' mlist <- move.matrix.gene(unipix, hometime= 0.1, type = "radiation")

move.matrix.gene <- function(unipix, hometime, type = "radiation"){
  # create a distance matrix
  distances = distm(unipix[, c("x", "y")])
  masterlist = list()
  for(i in 1:nrow(distances)){
    pophome = unipix[i, 3]
    distab = cbind(unipix[, 3], distances[i, ], 1:nrow(unipix))

    if(type == "exponential"){
      #runvec = 1 / exp(distab[, 2] / 1000)
      runvec = exp(-hometime * distab[, 2] / 100)
      masterlist[[i]] = runvec
    }

    if(type == "gravity"){
      #runvec = apply(distab, 1, function(x, y) x[1] * y / x[2], y = pophome)
      runvec = apply(distab, 1, function(x, y) x[1] * y / (exp(hometime * x[2] / 100)), y = pophome)
      # fix intra patch movements to the maximum
      runvec[!is.finite(runvec)] = max(runvec[is.finite(runvec)])
      masterlist[[i]] = runvec
    }

    if(type == "radiation"){
      distab2 = distab[order(distab[, 2]), ]
      distab2 = cbind(distab2, cumsum(distab2[, 1]))
      runvec = rep(NA, nrow(distab))
      for(k in 1:nrow(distab)){
        radiuspop = distab2[k, 4] - pophome - distab[k, 1]
        if(radiuspop < 0){radiuspop = 0}
        #runvec[k] = rad.model(pophome, distab2[k, 1], radiuspop)
        runvec[k] = rad.model(pophome, distab2[k, 1], radiuspop * hometime)
      }
      # reordering
      runvec2 = cbind(runvec, distab2[, 3])
      runvec2 = runvec2[order(runvec2[, 2]), ]
      masterlist[[i]] = runvec2[, 1]
    }
  }

  ## movement matrix modification for hometime parameter
  # - now hometiem is just an exponent, so just normalizes the movement matrix
  mm_list2 = masterlist
  for(i in 1:length(mm_list2)){
    # normalise matrix
    mm_list2[[i]] = mm_list2[[i]] / sum(mm_list2[[i]])
    #awaytime <- sum(mm_list2[[i]][!((1:length(mm_list2[[i]])) %in% i)])
    #mm_list2[[i]][i] = hometime
    #mm_list2[[i]][!((1:length(mm_list2[[i]])) %in% i)] = ((1 - hometime) / awaytime) * mm_list2[[i]][!((1:length(mm_list2[[i]])) %in% i)]
  }
  # precalculate the total number of person hours spent in each pixel
  # taking into account residents and visitors
  mm_list_cum = mm_list2
  for(i in 1:length(mm_list2)){mm_list_cum[[i]] = mm_list_cum[[i]] * unipix$pop[i]}
  mm_list_cum = Reduce("+", mm_list_cum)


  return(list(mm_list2, mm_list_cum))
}
