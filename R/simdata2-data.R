#' Simulated data
#'
#' Simulated data for cross-validated risk scores design with two outcomes (CVRS2).
#' @docType data
#'
#' @usage data(simdata2)
#'
#' @format A list of 4 data frames: patients, covar, response, response2. The list corresponds to the output of the \code{simulate.data2} function.
#'
#'         patients: A data frame with one row per patient and the following columns: 
#'         FID (family ID), IID (individual ID), , sens.resp.true (true sensitivity to response),
#'         sens.resp2.true (true sensitivity to response2), 
#'         cluster.true (1 for sens.pred.true ==1 and sens.pred2.true == 1, 2 for sens.pred.true ==0 and sens.pred2.true == 1, 
#'         3 for sens.pred.true ==1 and sens.pred2.true == 0 and 4 for sens.pred.true ==0 and sens.pred2.true == 0), 
#'         sens.true (1 for cluster.true == 1, 0 otherwise),
#'         treat (1 for treatment and 0 for control), rr (probability of a binary response 1), 
#'         rr2 (probability of a binary response 2)
#'
#'         covar: Covariate data for L covariates, one row per patient. 
#'
#'         response: A data frame of simulated binary response 2, one column per simulation run.
#'
#'         response2: A data frame of simulated binary response 2, one column per simulation run.
#'         
#'         The order of the patients in the data frames should match.
#' @keywords datasets
#'
#' @examples
#'
#' #"cvrs" method
#' data(simdata2)
#' nclust = 4
#' sig = 0.05
#' group.prop.sig = 0.2
#' seed = 123
#' plotrs = T
#' \donttest{simres.cvrs2 = analyse.simdata2 (simdata2, nclust, sig, group.prop.sig, seed, plotrs)}
"simdata2"

