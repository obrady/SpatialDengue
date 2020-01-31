#' Reported dengue cases in Singapore 2013-2016
#'
#' A dataset containing the number and location of reported dengue cases in Singapore each week for the period May 2013- June 2016. Cases
#'  within the same cluster (200m radius) are reported as one row with the number of cases detailed
#'
#' @format A data.frame composed of 32399 rows (each case cluster) and 6 columns
#' \describe{
#'   \item{Number_of_cases}{number of dengue cases identified in the week in that cluster}
#'   \item{Latitude}{Latitude of cluster}
#'   \item{Longitude}{Longitude of cluster}
#'   \item{Date}{Date of cluster report in format Y/M/D}
#'   \item{Month}{Month of cluster report in the year}
#'   \item{Week}{Week of cluster report since the beginign of the data period (23 May 2013)}
#' }
#' @source \url{http://outbreak.sgcharts.com}
"sgdat"