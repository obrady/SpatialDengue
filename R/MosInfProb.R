#' Calculate probability of mosquito infection
#'
#' Calculates the pixel-specific probability of infection taking into account time allocation of infectious humans across the landscape
#' @param sing_I matrix, Infectious human individuals in each pixel, rows = day of infection
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
#' betaEnv <- betaEnv.generate(bmean = 8, blogvariance = log(2), bcorrelation = -1, stim)
#' sing_I <- matrix(runif(15 * nrow(unipix), min = 0, max = 20), nrow = 15)
#' mlist <- move.matrix.gene(unipix, hometime = 0.2)
#' MIP <- mos.inf.prob(sing_I, betaEnv, unipix, mlist[[1]], mlist[[2]])
#' plot.state(MIP, sgpop, unipix)



mos.inf.prob <- function(sing_I, betaEnv, unipix, mm_list2, mm_list_cum){
  # run through each day of infectiousness (infectiousness of each human varies by day of infection)
  exp_prob<- rep(0, length(sing_I[1, ]))
  for(k in 1:15){
    # collate patches with infectious people on day k of infection
    cur_sing_I = sing_I[k, ]
    # infected pixel's details
    #inf_I_i = data.frame(unipix[cur_sing_I > 0, ], humIs = cur_sing_I[cur_sing_I > 0])
    inf_I_i = data.frame(patchID = 1:ncol(sing_I), humIs = cur_sing_I)
    inf_I_i = inf_I_i[which(cur_sing_I > 0), ]
    if(nrow(inf_I_i) > 0){
      # aggregate relevant precalculated movement matrices to calculate total
      # infectious person hours spent in other patches
      inf_exp_col = mm_list2[inf_I_i$patchID]
      inf_exp_col = Map("*", inf_exp_col, as.list(inf_I_i$humIs))
      inf_exp_col = Reduce("+", inf_exp_col)
      
      # calculate their infectiousness as a function of days since infection
      # NB scaling to make sure max(p(infection)) = 1 -i.e. perfectly infectiousness at onset of symptoms in human
      #infectiousness <- dnorm(k, rlnorm(length(inf_exp_col),  meanlog = log(5.9), sdlog = 0.045), 2) / 0.1982747
      #infectiousness <- dnorm(k, rep(5.9, nrow(unipix)), 2) / 0.1982747
      infectiousness <- dnbinom((k - 3), size = rep(29.374466, nrow(unipix)), mu = 4.171385) * 5.438625
      exp_prob = exp_prob + ((inf_exp_col / mm_list_cum) * betaEnv * infectiousness)
    }
  }
  # make sure cumulative probability of infection not > 1
  exp_prob[exp_prob > 1] = 1
  return(exp_prob)
}