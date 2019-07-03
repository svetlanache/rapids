#' @title 
#' Compute p-value for the overall treatment effect and for the 
#' interaction effect between the treatment and the predicted sensitivity status, for real data.
#'
#' @description
#' The sensitive group is identified according to the input method ("cvasd" for the cross-validated adaptive signature design or "cvrs" for the cross-validated risk scores design),
#' significance level of the test for the sensitive group is 
#' determined by the input proportion of the significance  
#' level for the sensitive group (group.prop.sig) out of the overall significance (sig).
#' The p-value for the overall treatment effect is computed using logistic regression with treatment allocation as a predictor.
#' The p-value for the interaction effect between the treatment and the predicted sensitivity status is computed using logistic regression with the treatment effect, the effect for the predicted sensitivity status and the interaction effect. 
#'
#' @param datalist A list made up of 3 data frames (patients and covar) and a vector of binary responses, for real data (see \code{realdata} object): patients (a data frame with patients information), covar (a data frame with the covariates), response (a vector of responses).
#' @param  sig An overall significance level for the adaptive design.
#' @param  group.prop.sig Proportion of significance level for the sensitive group test.
#' @param  method  "cvasd" (for the cross-validated adaptive signature design) or "cvrs" (for the cross-validated risk scores design).
#' @param  eta A significance level for covariate-wise logistic regression for the "cvasd" method (double, or vector of doubles).
#' @param  R A threshold of the odds ratio for the "cvasd" method (double, or vector of doubles).
#' @param  G A threshold for the number of covariates for the "cvasd" method (integer, or vector of integers).
#' @param  seed A  seed for the random number generator.
#' @param  plotrs An indicator whether to plot the risk scores for the "cvrs" method (default = FALSE).
#'
#' @return An object of class \code{"rapids"}.
#' @return patients: A data frame with patients inormation.
#' @return pval.overall: P-value for the overall test for the treatment effect.
#' @return stat.group: Statistic for the interaction effect between the treatment and the predicted sensitivity status.
#' @return pval.group: P-value for the interaction effect between the treatment and the predicted sensitivity status.
#' @return sens.pred: Predicted sensitivity status (rows = patients, columns = simulations).
#' @return estimate.rr: Estimated response rate in the sensitive group on the treatment arm.
#' @return cvrs: A matrix of the risk scores (rows = patients, columns = simulations), for "cvrs" method only.
#' @return eta,R,G: Significance level, threshold of the odds ratio and threshold for the number of covariates for covariate-wise logistic regression, for"cvasd" method only.
#' @author Svetlana Cherlin, James Wason
#'
#' @examples
#' #"cvrs" method
#' sig = 0.05
#' group.prop.sig = 0.2
#' method = "cvrs"
#' eta = NULL
#' R = NULL
#' G = NULL
#' seed = 123
#' plotrs = T
#' realres.cvrs = analyse.realdata(realdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)
#'
#' #"cvasd" method
#' sig = 0.05
#' group.prop.sig = 0.2
#' method = "cvasd"
#' eta = c(0.01, 0.02, 0.03)
#' R = c(2.5, 2, 1.5)
#' G = c(3,2,1)
#' seed = 123
#' plotrs = F
#' realres.cvasd = analyse.realdata(realdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)
#' @seealso
#' \code{\link{analyse.simdata}}, \code{\link{simulate.data}} and \code{\link{cvrs.plot}} functions; \code{\link{print}} and \code{\link{plot}} methods.
#' @export 


analyse.realdata = function(datalist, sig = 0.05, group.prop.sig = 0.2, 
                            method = c("cvrs", "cvasd"), eta = NULL, R = NULL, G = NULL, 
                            seed = NULL, plotrs = F)
                           
{  

   patients = datalist$patients
   covar = datalist$covar
   response = datalist$response
   
   if (is.null(dim(patients))) {
      stop ("Patients data is missing.")
   } else if (ncol(patients) < 3 ) {
      stop ("Patients data must have 3 columns (FID, IID, treatment allocation).")
   }
   if (is.null(length(response))) {
      stop ("Response data is missing.")
   } 
   if (is.null(dim(covar))) {
      stop ("Covariate data is missing.")
   } 
   if (!((length(response) == nrow(patients)) & (nrow(patients) == nrow(covar)))) {
      stop ("Patients data, covariate data and response data dimensions do not match.")
   }
   if (!(is.numeric(sig) & (length(sig) == 1) & (sig > 0) & (sig < 1))) {
      stop ("Significance level sould be between 0 and 1.")
   }
   if (!(is.numeric(group.prop.sig) & (length(group.prop.sig) == 1) & (group.prop.sig > 0) & (group.prop.sig < 1))) {
      stop ("Proportion of significance level for the sensitive group shoud be between 0 and 1.")
   }   

   method = match.arg(method)
   sig.group = sig*group.prop.sig
   sig.overall = sig - sig.group   

   ## P-value for the overall arm comparison
   mod = glm(response ~ patients$treat, family = "binomial")
   theta = mod$coeff
   sum.theta = summary(mod)$coeff
   pval.overall = sum.theta[names(theta) == "patients$treat", colnames(sum.theta) == "Pr(>|z|)"]

   # P-value for the sensitive group test
      # Find sensitive group according to the input method

   sens = switch (method, cvasd = sens.cvasd(patients, covar, response, eta, R, G, seed),
                          cvrs = sens.cvrs(patients, covar, response, seed))                                      
   resp.group = response[sens$sens.pred == 1] #response for sensitive group
   treat.group = patients$treat[sens$sens.pred == 1] #treatment allocation for sensitive group
         
   mod = glm(response ~ sens$sens.pred + patients$treat + sens$sens.pred:patients$treat, family = "binomial")
   theta = mod$coeff
   sum.theta = summary(mod)$coeff
   stat.group = sum.theta[names(theta) == "sens$sens.predTRUE:patients$treat", colnames(sum.theta) == "z value"]
   pval.group = sum.theta[names(theta) == "sens$sens.predTRUE:patients$treat", colnames(sum.theta) == "Pr(>|z|)"]

   resp.treat =  sum(resp.group & treat.group)   #number of responders in the treatment arm
   nonresp.treat =  sum(!resp.group & treat.group) #number of non-responders in the treatment arm
   estimate.rr = resp.treat/(resp.treat + nonresp.treat) #estimated response rate in the sensitive group in the treatment arm
   
   output = switch(method, cvasd = list(patients = patients, pval.overall = pval.overall, stat.group = stat.group, pval.group = pval.group, estimate.rr = estimate.rr, sens.pred = sens$sens.pred, eta = eta, R = R, G = G),
            cvrs = list(patients = patients, pval.overall = pval.overall, stat.group = stat.group, pval.group = pval.group, estimate.rr = estimate.rr, sens.pred = sens$sens.pred, cvrs = sens$cvrs)) 
   if (method == "cvrs" & plotrs) cvrs.plot(sens$cvrs, sens$sens.pred)

   class(output) = "rapids"

   return(output)
}
