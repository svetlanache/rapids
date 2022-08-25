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
   stat.group = NA
   pval.group = NA
   msg = NULL

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
   if(nrow(sum.theta) > 2) {
      stat.group = sum.theta[names(theta) == "sens$sens.predTRUE:patients$treat", colnames(sum.theta) == "z value"]
      pval.group = sum.theta[names(theta) == "sens$sens.predTRUE:patients$treat", colnames(sum.theta) == "Pr(>|z|)"]
   }else {
      msg = "Sensitive group is not found"
   }

   resp.treat =  sum(resp.group & treat.group)   #number of responders in the treatment arm
   nonresp.treat =  sum(!resp.group & treat.group) #number of non-responders in the treatment arm
   estimate.rr = resp.treat/(resp.treat + nonresp.treat) #estimated response rate in the sensitive group in the treatment arm

   output = switch(method, cvasd = list(patients = patients, pval.overall = pval.overall, stat.group = stat.group, pval.group = pval.group, estimate.rr = estimate.rr, sens.pred = sens$sens.pred, eta = eta, R = R, G = G, msg = msg),
            cvrs = list(patients = patients, pval.overall = pval.overall, stat.group = stat.group, pval.group = pval.group, estimate.rr = estimate.rr, sens.pred = sens$sens.pred, cvrs = sens$cvrs, msg = msg))
   if (method == "cvrs" & plotrs) cvrs.plot(sens$cvrs, sens$sens.pred)

   class(output) = "rapids"

   return(output)
}


#' @title
#' Compute p-value for the overall treatment effect and for the
#' interaction effect between the treatment and the predicted sensitivity status, for real data.
#'
#' @description
#' significance level of the test for the sensitive group is
#' determined by the input proportion of the significance
#' level for the sensitive group (group.prop.sig) out of the overall significance (sig).
#' The p-value for the overall treatment effect is computed using bivariate logistic regression with treatment allocation as a predictor.
#' The p-value for the interaction effect between the treatment and the predicted sensitivity status is computed using bivariate logistic regression with the treatment effect, the effect for the predicted sensitivity status and the interaction effect.
#'
#' @param datalist A list made up of 3 data frames (patients and covar) and a vector of binary responses, for real data (see \code{realdata} object): patients (a data frame with patients information), covar (a data frame with the covariates), response (a vector of responses).
#' @param  nclust Number of clusters
#' @param  sig An overall significance level for the adaptive design.
#' @param  group.prop.sig Proportion of significance level for the sensitive group test.
#' @param  seed A  seed for the random number generator.
#' @param  plotrs An indicator whether to plot the risk scores (default = FALSE).
#'
#' @return An object of class \code{"rapids"}.
#' @return patients: A data frame with patients inormation.
#' @return pval.resp: P-value for the overall test for the treatment effect w.r.t. response 1
#' @return pval.resp2: P-value for the overall test for the treatment effect w.r.t. response 2
#' @return stat.resp.group: Statistic for the interaction effect between the treatment  and the predicted cluster (w.r.t. response ), a vector of length nclust
#' @return pval.resp.group: P-value for the interaction effect between the treatment and the and the predicted cluster (w.r.t. response 2), a vector of length nclust
#' @return stat.resp2.group: Statistic for the interaction effect between the treatment and the predicted cluster (w.r.t. response), a vector of length nclust
#' @return pval.resp2.group: P-value for the interaction effect between the treatment and the predicted cluster (w.r.t. response 2), a vector of length nclust

#' @return estimate.rr: Estimated rate of responsein the sensitive group on the treatment arm.
#' @return estimate.rr2: Estimated rate of response 2 in the sensitive group on the treatment arm.
#' @return cvrs: A matrix of the risk scores (rows = patients, columns = simulations).
#' @return cvrs2: A matrix of the risk scores 2(rows = patients, columns = simulations).
#' @return cluster.pred: Predicted clusters

#' @author Svetlana Cherlin, James Wason

#' @examples
#' sig = 0.05
#' nclust = 4
#' group.prop.sig = 0.2
#' seed = 123
#' plotrs = T
#' realres.cvrs2 = analyse.realdata2(realdata2, nclust, sig, group.prop.sig, seed, plotrs)
#'
#' @seealso
#' \code{\link{analyse.simdata2}}, \code{\link{simulate.data2}} and \code{\link{cvrs.plot}} functions; \code{\link{print}} and \code{\link{plot}} methods.
#' @export
#' @importFrom "ggplot2" "ggsave"


analyse.realdata2 = function(datalist, nclust, sig = 0.05, group.prop.sig = 0.2,
                            seed = NULL, plotrs = F)
{
   patients = datalist$patients
   covar = datalist$covar
   response = datalist$response
   response2 = datalist$response2

   if (is.null(dim(patients))) {
      stop ("Patients data is missing.")
   } else if (ncol(patients) < 3 ) {
      stop ("Patients data must have 3 columns (FID, IID, treatment allocation).")
   }
   if (is.null(length(response))) {
      stop ("Response data is missing.")
   }
   if (is.null(length(response2))) {
      stop ("Response2 data is missing.")
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

   sig.group = sig*group.prop.sig
   sig.overall = sig - sig.group

   ## P-value for the overall arm comparison
   mod = tryCatch({vglm(cbind(response, response2) ~ patients$treat, binom2.or)}, error = function(e) e, warning = function(w) w)
   if (is(mod, "error") | is(mod, "warning")) {
      pval.resp = 1
      pval.resp2 = 1
   } else if (is.na(logLik(mod))) {
      pval.resp = 1
      pval.resp2 = 1
   } else {
      sum.mtx = coef(summary(mod))
      pval.resp = sum.mtx[rownames(sum.mtx) == "patients$treat:1", colnames(sum.mtx) == "Pr(>|z|)"]
      pval.resp2 = sum.mtx[rownames(sum.mtx) == "patients$treat:2", colnames(sum.mtx) == "Pr(>|z|)"]
   }

   estimate.rr = numeric (length = nclust)
   estimate.rr2 = numeric (length = nclust)
   stat.resp.group = numeric (length = nclust)
   stat.resp2.group = numeric (length = nclust)
   pval.resp.group = numeric (length = nclust)
   pval.resp2.group = numeric (length = nclust)

   ## P-value for the sensitive group test
   sens = sens.cvrs2(patients, covar, response, response2, seed, nclust)

  for (i in 1:nclust) { ## loop on clusters
     sens.pred = numeric(length = nrow(patients))
     sens.pred[sens$cluster.pred == i] = 1

     mod = tryCatch({vglm(cbind(response, response2) ~ sens.pred + patients$treat + sens.pred:patients$treat, binom2.or)}, error = function(e) e, warning = function(w) w)
     if (is(mod, "error") | is(mod, "warning")) {
        stat.resp.group[i] = -999
        pval.resp.group[i] = 1
        stat.resp2.group[i] = -999
        pval.resp2.group[i] = 1
     } else if (is.na(logLik(mod))) {
        stat.resp.group[i] = -999
        pval.resp.group[i] = 1
        stat.resp2.group[i] = -999
        pval.resp2.group[i] = 1
     } else {
        sum.mtx = coef(summary(mod))
        stat.resp.group[i] = sum.mtx[rownames(sum.mtx) == "sens.pred:patients$treat:1", colnames(sum.mtx) == "z value"]
        pval.resp.group[i] = sum.mtx[rownames(sum.mtx) == "sens.pred:patients$treat:1", colnames(sum.mtx) == "Pr(>|z|)"]
        stat.resp2.group[i] = sum.mtx[rownames(sum.mtx) == "sens.pred:patients$treat:2", colnames(sum.mtx) == "z value"]
        pval.resp2.group[i] = sum.mtx[rownames(sum.mtx) == "sens.pred:patients$treat:2", colnames(sum.mtx) == "Pr(>|z|)"]
     }
     resp.group = response[sens.pred == 1] #response for sensitive group
     resp2.group = response2[sens.pred == 1] #response2 for sensitive group
     treat.group = patients$treat[sens.pred == 1] #treatment allocation for sensitive group

     resp.treat =  sum(resp.group & treat.group)   #number of responders in the treatment arm
     nonresp.treat =  sum(!resp.group & treat.group) #number of non-responders in the treatment arm
     estimate.rr[i] = resp.treat/(resp.treat + nonresp.treat) #estimated response rate in the sensitive group in the treatment arm

     resp2.treat =  sum(resp2.group & treat.group)   #number of responders2 in the treatment arm
     nonresp2.treat =  sum(!resp2.group & treat.group) #number of non-responders2 in the treatment arm
     estimate.rr2[i] = resp2.treat/(resp2.treat + nonresp2.treat) #estimated response rate in the sensitive group in the treatment arm
  } ## end loop on clusters

output = list(patients = patients, pval.resp = pval.resp, pval.resp2 = pval.resp2, stat.resp.group = stat.resp.group, pval.resp.group = pval.resp.group, stat.resp2.group = stat.resp2.group, pval.resp2.group = pval.resp2.group, estimate.rr = estimate.rr, estimate.rr2 = estimate.rr2, cvrs = sens$cvrs, cvrs2 = sens$cvrs2, cluster.pred = sens$cluster.pred)
   if (plotrs) {
	   plotcvrs = cvrs2.plot(sens$cvrs, sens$cvrs2, sens$cluster.pred, NULL, FALSE)
	   ggsave ("cvrs2.pdf", plotcvrs)
   }
   class(output) = "rapids"
   return(output)
}

