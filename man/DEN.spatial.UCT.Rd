% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DENSpatialUCT.R
\name{DEN.spatial.UCT}
\alias{DEN.spatial.UCT}
\title{Simulate an urban cluster trial of dengue prophylactic drugs}
\usage{
DEN.spatial.UCT(weekdates, fitdat, pastdat, unipix, pixdistmat, UCTinfo,
  nruns, paramsList = NULL)
}
\arguments{
\item{weekdates}{Two element vector of the start and end weeks of the simulation over which the mdoel will be evaluated over}

\item{fitdat}{Data frame of the locations, numbers and timings (in weeks) of cases to fit the model to, see ?sgdat}

\item{pastdat}{Data frame of the locations, numbers and timings (in weeks) of all cases in the dataset (is used to generate the starting immunity profile), see ?sgdat}

\item{unipix}{Universal pixel lookup table, see ?make.unipix}

\item{pixdistmat}{A patch distance matrix, see example}

\item{UCTinfo}{A data frame including details on Urban Cluster Trial tart time, delay between enrollment and followup, number of treatment and control clustes and drug effective coverage}

\item{nruns}{integer, number of stochastic runs of the model over which results should be averaged}

\item{paramsList}{Optional parameter list. If not supplied returns to defaults, see tutorial for full parameter list, see ?model.run for full list and explanation of parameters}
}
\description{
Provide basic information on the UCT and the number of stochastic runs over which model runs should be averaged
}
\details{
runs "nruns" number of stochastic simulations and enrolls treatment or control clusters once the trail has begun then stops enrolling once the target
number of clusters has been reached. Returns the normal DEN.spatial sumamry results but with "treatlog" and "contlog" which detail which patches were enrolled into treatment or control arms at which point
}
\examples{
data(sgdat)
data(sgpop)
sgpop <- pop.process(sgpop, agg = 10)
unipix <- make.unipix(sgpop)
pixdistmat <- distm(cbind(unipix$x, unipix$y))
sgdat <- data.frame(sgdat, patchID = apply(cbind(sgdat[, 3:2]), 1, pix.id.find, unipix))
weekdates <- c(40, 92)

UCTinfo <- list(Tstart = 50,
               Tfollowup = 28,
               ntreat = 30,
               ncontrol = 30,
               DrugEfficacy = 0.9)

denmod_UCT = DEN.spatial.UCT(weekdates, sgdat, sgdat, unipix, pixdistmat, UCTinfo, nruns = 10)
effectiveness <- UCT.Eff.calc(denmod_UCT, unipix)
}
\keyword{UCT}
\keyword{model}
\keyword{simulation}
