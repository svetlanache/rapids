#' Simulated data
#'
#' Simulated data for cross-validated adaptive signature design/cross-validated risk scores design.
#' @docType data
#'
#' @usage data(simdata)
#'
#' @format A list of 3 data frames: patients, covar, response. The list corresponds to the output of the \code{simulate.data} function.
#'
#'         patients: A data frame with one row per patient and the following columns: 
#'         FID (family ID), IID (individual ID), treat (1 for the treatment arm and 0 for the control arm), rr (probability of a binary response).
#'
#'         covar: Covariate data for L covariates, one row per patient. 
#'
#'         response: A data frame of simulated binary responses, one column per simulation run.
#'         
#'         The order of the patients in the data frames should match.
#' @keywords datasets
#'
#' @examples
#'
#' #"cvrs" method
#' data(simdata)
#' sig = 0.05
#' group.prop.sig = 0.2
#' seed = 123
#' plotrs = T
#' method = "cvrs"
#' eta = NULL
#' R = NULL
#' G = NULL 
#' \donttest{simres.cvrs = analyse.simdata(simdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)}
#'
#' #"cvasd" method
#' data(simdata)
#' sig = 0.05
#' group.prop.sig = 0.2
#' seed = 123
#' method = "cvasd"
#' seed = 123
#' plotrs = F
#' eta = c(0.01,0.02,0.03)
#' R = c(2,2.5,1.2)
#' G = c(3,2,1) 
#' \donttest{simres.cvasd = analyse.simdata(simdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)}
"simdata"

