#' @title 
#' Compute the power for the cross-validate adaptive signature design/cross-validated risk scores design for simulated data.
#'
#' @description  The power of the adaptive design 
#' is computed according to the input method ("cvasd" for cross-validated adaptive signature design or "cvrs" for cross-validated risk scores design),
#' It is made up of the power for the overall
#' test and the power for the sensitive group test.
#' The significance level of the test for the sensitive group is 
#' determined by the input proportion of the significance  
#' level for the sensitive group (group.prop.sig) out of the overall significance (sig).
#' For the "cvasd" method, the input tuning set \code{eta,R,G} might be comprised of vectors. In this case,
#' the tuning is performed using nested cross-validation. For each outer cross-validation fold, only one inner cross-validation fold is used, to save the computational time.
#'
#' @param datalist A list of 3 data frames that corresponds to the output of the \code{simulate.data} function (see \code{simdata} object): patients (a data frame with patients inormation), covar (a data frame with the covariates), response (a data frame of simulated responses).
#' @param  sig An overall significance level for the adaptive design.
#' @param  group.prop.sig Proportion of significance level for the sensitive group test.
#' @param  method  "cvasd" (for cross-validated adaptive signature design) or "cvrs" (for cross-validated risk scores design).
#' @param  runs  Number of simulation runs.
#' @param  eta A significance level for the covariate-wise logistic regression for the "cvasd" method (double, or vector of doubles).
#' @param  R A threshold of the odds ratio for the "cvasd" method (double, or vector of doubles).
#' @param  G A threshold for the number of covariates for the "cvasd" method (integer, or vector of integers).
#' @param  seed A  seed for the random number generator.
#' @param  plotrs An indicator whether to plot the risk scores for the "cvrs" method (default = FALSE).
#'
#' @return An object of class \code{"rapids"}.
#' @return patients: A data frame with patients inormation.
#' @return pwr.overall: Power for the overall test.    
#' @return pwr.group: Power for the sensitive group test.
#' @return pwr.adaptive: Power for the adaptive design.
#' @return estimate.rr: Empirical response rate on the experimental arm for the sensitive group, one value per simulation run.
#' @return psens: Sensitivity of identifying the sensitive group, one value per simulation run.
#' @return pspec: Specificity of identifying the sensitive group, one value per simulation run.
#' @return sens.pred: Predicted sensitivity status (rows = patients, columns = simulations).
#' @return cvrs: A matrix of the risk scores (rows = patients, columns = simulations), for the "cvrs" method only.
#' @return eta,R,G: Significance level, threshold of the odds ratio and threshold for the number of covariates for covariate-wise logistic regression, for "cvasd" method only.
#'
#' @examples
#' #"cvrs" method
#' data(simdata)
#' sig = 0.05
#' group.prop.sig = 0.2
#' method = "cvrs"
#' seed = 123
#' plotrs = T
#' eta = NULL
#' R = NULL
#' G = NULL
#' simres.cvrs = analyse.simdata (simdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)
#'
#' #"cvasd" method
#' data(simdata)
#' sig = 0.05
#' group.prop.sig = 0.2
#' method = "cvasd"
#' seed = 123
#' plotrs = T
#' eta = c(0.01, 0.02, 0.03)
#' R = c(2.5, 2, 1.5)
#' G = c(3,2,1)
#' simres.cvasd = analyse.simdata (simdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)

#' @seealso
#' \code{\link{analyse.realdata}}, \code{\link{simulate.data}} and \code{\link{cvrs.plot}} functions; \code{\link{print}} and \code{\link{plot}} methods.
#' @author Svetlana Cherlin, James Wason
#' @export 
#'
analyse.simdata <- function(datalist, sig = 0.05, group.prop.sig = 0.2, 
                            method = c("cvrs", "cvasd"), eta = NULL, R = NULL, G = NULL, 
                            seed = NULL, plotrs = F) 
{  

   patients = datalist$patients
   covar = datalist$covar
   response = datalist$response

   if (is.null(dim(patients))) {
      stop ("Patients data is missing.")
   } else if (ncol(patients) < 4 ) {
      stop ("Patients data must have 4 columns (FID, IID, sensitivity status, treatment allocation).")
   }
   if (is.null(dim(response))) {
      stop ("Response data is missing.")
   } 

   if (is.null(dim(covar))) {
      stop ("Covariate data is missing.")
   } 
   if (!((nrow(response) == nrow(patients)) & (nrow(patients) == nrow(covar)))) {
      stop ("Patients data, covariate data and response data must have the same number of rows.")
   }
   if (!(is.numeric(sig) & (length(sig) == 1) & (sig > 0) & (sig < 1))) {
      stop ("Significance level sould be between 0 and 1.")
   }
   if (!(is.numeric(group.prop.sig) & (length(group.prop.sig) == 1) & (group.prop.sig > 0) & (group.prop.sig < 1))) {
      stop ("Proportion of significance level for the sensitive group shoud be between 0 and 1.")
   }   

   method = match.arg(method)    
   runs = ncol(response)
   sig.group = sig*group.prop.sig
   sig.overall = sig - sig.group   

   ## Power for the overall arm comparison
   N = nrow(response)
   N.treat = nrow(patients[patients$treat == 1,])
   N.con = nrow(patients[patients$treat == 0,])
   pval = apply (response, 2, function (y)
   {
            prop.test(x = c(sum(y[patients$treat == 0]), sum(y[patients$treat == 1])),
                   n = c(N.con, N.treat), alternative = "two.sided")$p.value
   })
   pwr.overall = sum(pval < sig.overall)/runs  

   ## Power for the sensitive group test
   res = apply(response, 2, function(x) 
   {
      # Find sensitive group according to the input methods
      sens = switch (method, cvasd = sens.cvasd(patients, covar, x, eta, R, G, seed),
                             cvrs = sens.cvrs(patients, covar, x, seed))                                      
      x.group = x[sens$sens.pred == 1]                  #response for sensitive group
      treat.group = patients$treat[sens$sens.pred == 1] #treatment allocation for sensitive group
 
      conf = matrix(nrow = 2, ncol = 2, 
             data = c(sum(x.group & !treat.group),  #number of responders in control
                      sum(!x.group & !treat.group), #number of non-responders in control
                      sum(x.group & treat.group),   #number of responders in treatment
                      sum(!x.group & treat.group)), #number of non-responders in treatment
             byrow = TRUE) 

      p.value = fisher.test(conf, alternative = "two.sided")$p.value
      estimate.rr = conf[2,1]/(conf[2,1] + conf[2,2])

      list(pval.group = p.value, psens = sens$psens, pspec = sens$pspec, estimate.rr = estimate.rr, sens.pred = sens$sens.pred,
          cvrs = switch(method, cvrs = sens$cvrs, cvasd = NA))
   })
   pval.group = numeric()
   psens = numeric()
   pspec = numeric()
   estimate.rr = numeric()
   cvrs = matrix(nrow = N, ncol = runs)
   rownames(cvrs) = rownames(patients)
   colnames(cvrs) = 1:ncol(response)
   sens.pred = matrix(nrow = N, ncol = runs)
   rownames(sens.pred) = rownames(patients)
   colnames(sens.pred) = 1:ncol(response)
   for (i in seq(runs)) {
      pval.group = c(pval.group, res[[i]]$pval.group)
      psens = c(psens, res[[i]]$psens)
      pspec = c(pspec, res[[i]]$pspec)
      estimate.rr = c(estimate.rr, res[[i]]$estimate.rr)
      sens.pred[,i] = res[[i]]$sens.pred
      if (method == "cvrs") cvrs[,i] = res[[i]]$cvrs
   }
   pwr.group = sum(pval.group < sig.group)/runs

   ## Overall power of the adaptive design
   pwr.adaptive = pwr.overall + (1-pwr.overall)*pwr.group

   output = switch(method, cvasd = list(patients = patients, pwr.overall = pwr.overall, pwr.group = pwr.group, pwr.adaptive = pwr.adaptive, estimate.rr = estimate.rr, psens = psens, pspec = pspec, sens.pred = sens.pred, eta = eta, R = R, G = G),
                   cvrs = list(patients = patients, pwr.overall = pwr.overall, pwr.group = pwr.group, pwr.adaptive = pwr.adaptive, estimate.rr = estimate.rr, psens = psens, pspec = pspec, sens.pred = sens.pred, cvrs = cvrs)) 

   if (method == "cvrs" & plotrs) cvrs.plot(cvrs, sens.pred)

   class(output) = "rapids"

   return(output)
}
