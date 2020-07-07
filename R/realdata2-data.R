#' Hypothetical real data
#'
#' Hypothetical real data for cross-validated risk scores design with two outcomes (CVRS2).
#' @docType data
#'
#' @usage data(realdata2)
#'
#' @format A list of 2 data frames with patients and covariates data, and a vector of responses.
#'
#'         patients: A data frame with one row per patient and the following columns: 
#'         FID (family ID), IID (individual ID), treat (1 for the treatment arm and 0 for the control arm).
#'
#'         covar: Covariate data for L covariates, one row per patient. 
#'
#'         response: A vector of hypothetical true binary response 1. 
#'        
#'         response2: A vector of hypothetical true binary response 1.
#'
#'         The order of the patients in the data frames and the vector of responses should match.
#' @keywords datasets
#'
#' @examples
#' data(realdata2)
#' sig = 0.05
#' nclust = 4
#' group.prop.sig = 0.2
#' seed = 123
#' plotrs = T
#' \donttest{realres.cvrs2 = analyse.realdata2(realdata2, nclust, sig, group.prop.sig, seed, plotrs)}
"realdata2"
