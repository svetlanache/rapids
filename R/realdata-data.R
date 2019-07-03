#' Hypothetical real data
#'
#' Hypothetical real data for cross-validated adaptive signature design/cross-validated risk scores design.
#' @docType data
#'
#' @usage data(realdata)
#'
#' @format A list of 2 data frames with patients and covariates data, and a vector of responses.
#'
#'         patients: A data frame with one row per patient and the following columns: 
#'         FID (family ID), IID (individual ID), treat (1 for the treatment arm and 0 for the control arm).
#'
#'         covar: Covariate data for L covariates, one row per patient. 
#'
#'         response: A vector of hypothetical true binary responses. 
#'
#'         The order of the patients in the data frames and the vector of responses should match.
#' @keywords datasets
#'
#' @examples
#' #"cvrs" method
#' data(realdata)
#' sig = 0.05
#' group.prop.sig = 0.2
#' method = "cvrs"
#' seed = 123
#' plotrs = T
#' eta = NULL
#' R = NULL
#' G = NULL 
#' \donttest{realres.cvrs = analyse.realdata(realdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)}
#'
#' #"cvasd" method
#' data(realdata)
#' sig = 0.05
#' group.prop.sig = 0.2
#' method = "cvasd"
#' seed = 123
#' plotrs = F
#' eta = c(0.1, 0.2, 0.3)
#' R = c(2, 1.5, 1.2)
#' G = c(3, 2,1)
#' \donttest{realres.cvasd = analyse.realdata(realdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)}


"realdata"
