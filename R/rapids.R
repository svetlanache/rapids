#' rapids: A package for cross-validated adaptive signature design (CVASD),
#'         cross-validation risk scores design (CVRS), and two outcome 
#' cross-validation risk scores design (CVRS2)
#'
#' The "rapids" package provides three main funcitons: \code{simulate.data} \code{analyse.simdata} and \code{analyse.realdata}.
#' Additional functions are \code{permutation.test} for permutation tests for the real data, \code{cvrs.plot} for plotting the risk scores, and also \code{print} and \code{plot} generic methods.
#' 
#' @section simulate.data:
#' The function simulates covariates data and binary responses to be used in the analysis of the cross-validated adaptive signature design or cross-validation risk scores design.
#'
#' @section analyse.simdata:
#' The function computes  the power of the design for the simulated data according to the input method ("cvasd" or "cvrs").
#'
#' @section analyse.realdata:
#' The function computes the p-value for the interaction effect between the treatment and the sensitivity status. 
#' The sensitivity status is predicted according to the input method ("cvasd" and "cvrs").
#'
#' @section permutation.test:
#' The function performs permutation test for the real data.
#'
#' @section cvrs.plot:
#' The function plots the risk scores for the "cvrs" method.

#' @section simulate.data2:
#' The function simulates covariates data and two binary responses to be used in the analysis of the cross-validation risk scores design for two outcomes (CVRS2)/
#'
#' @section analyse.simdata2:
#' The function computes  the power of the design for the simulated data for the CVRS2.
#'
#' @section analyse.realdata2:
#' The function computes the p-value for the interaction effect between the treatment and the sensitivity status for the CVRS2 design.

#' @docType package
#' @name rapids
NULL
