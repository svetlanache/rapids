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
#' @importFrom "VGAM" "coef"
#' @importFrom "VGAM" "summary"
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

#' @title 
#' Compute the power for the cross-validate risk scores design for simulated data with two outcomes (CVRS2)
#'
#' @description The power is computed for the two outcome variables. For each outcome variable, the power of the adaptive design is comprise of the power for the overall test and  test and the power for the sensitive group test.
#' The significance level of the test for the sensitive group is 
#' determined by the input proportion of the significance  
#' level for the sensitive group (group.prop.sig) out of the overall significance (sig).
#'
#' @param datalist A list of 3 data frames that corresponds to the output of the \code{simulate.data2} function (see \code{simdata2} object): patients (a data frame with patients inormation), covar (a data frame with the covariates), response and response2 (data frames of simulated responses).
#' @param nclust Number of clusters
#' @param  sig An overall significance level for the adaptive design.
#' @param  group.prop.sig Proportion of significance level for the sensitive group test.
#' @param  runs  Number of simulation runs.
#' @param  seed A  seed for the random number generator.
#' @param  plotrs An indicator whether to plot the risk scores(default = FALSE).
#'
#' @return An object of class \code{"rapids"}.
#' @return patients: A data frame with patients inormation.
#' @return pwr.overall: Power for the overall test.
#' @return pwr.group: Power for the sensitive group test, a vector of length nclust
#' @return pwr.adaptive: Power for the adaptive design, a vector of length nclust
#' @return estimate.rr: Empirical rate of response on the experimental arm for the sensitive group, one value per simulation run.
#' #return estimate.rr2: Empirical rate of response2 on the experimental arm for the sensitive group, one value per simulation run.
#' @return psens: Sensitivity of identifying the sensitive group, a vector of the length nclust per simulation run.
#' @return pspec: Specificity of identifying the sensitive group, a vector of the length nclust per simulation run.
#' @return cluster.pred: A matrix with predicted clusters (rows = patients, columns = simulations).
#' @return cvrs: A matrix of the risk scores for the response(rows = patients, columns = simulations).
#' @return cvrs2: A matrix for the risk scores for the response2 (rows = patients, columns = simulations).
#'
#' @examples
#' data(simdata2)
#' nclust = 4
#' sig = 0.05
#' group.prop.sig = 0.2
#' seed = 123
#' plotrs = T
#' simres.cvrs2 = analyse.simdata2 (simdata2, nclust, sig, group.prop.sig, seed, plotrs)
#'
#' @seealso
#' \code{\link{analyse.realdata2}}, \code{\link{simulate.data2}} and \code{\link{cvrs2.plot}} functions; \code{\link{print}} and \code{\link{plot}} methods.
#' @author Svetlana Cherlin, James Wason
#' @importFrom "VGAM" "vglm"
#' @importFrom "VGAM" "binom2.or"
#' @importFrom "ggplot2" "ggsave"
#' @export
#'

analyse.simdata2 <- function(datalist, nclust, sig = 0.05, group.prop.sig = 0.2,
                             seed = NULL, plotrs = F)
{
   patients = datalist$patients
   covar = datalist$covar
   response = datalist$response
   response2 = datalist$response2
   runs = ncol(response)

   if (is.null(dim(patients))) {
      stop ("Patients data is missing.")
   } else if (ncol(patients) < 4 ) {
      stop ("Patients data must have 4 columns (FID, IID, sensitivity status, treatment allocation).")
   }
   if (is.null(dim(response))) {
      stop ("Response data is missing.")
   }
   if (is.null(dim(response2))) {
      stop ("Response2 data is missing.")
   }
   if (is.null(dim(covar))) {
      stop ("Covariate data is missing.")
   }
   for (i in 1:runs) {
      if (!((nrow(response) == nrow(patients)) & (nrow(patients) == nrow(covar[i,,])) & (nrow(covar[i,,]) == nrow(response2)))) {
         stop ("Patients data, covariate data, response and response2 data must have the same number of rows.")
      }
   }
   if (!(is.numeric(sig) & (length(sig) == 1) & (sig > 0) & (sig < 1))) {
      stop ("Significance level sould be between 0 and 1.")
   }
   if (!(is.numeric(group.prop.sig) & (length(group.prop.sig) == 1) & (group.prop.sig > 0) & (group.prop.sig < 1))) {
      stop ("Proportion of significance level for the sensitive group shoud be between 0 and 1.")
   }

   N = nrow(patients)
   sig.group = sig*group.prop.sig
   sig.overall = sig - sig.group

   ## Power for the overall arm comparison 
   pval.resp = numeric(length = runs)
   pval.resp2 = numeric(length = runs)
   for (i in 1:runs) {
      mod = tryCatch({vglm(cbind(response[,i], response2[,i]) ~ patients$treat, binom2.or)}, error = function(e) e, warning = function(w) w)
      if (is(mod, "error") | is(mod, "warning")) {
         pval.resp[i] = 1
         pval.resp2[i] = 1
      } else if (is.na(logLik(mod))) {
         pval.resp[i] = 1
         pval.resp2[i] = 1
      } else {
         sum.mtx = coef(summary(mod))
         pval.resp[i] = sum.mtx[rownames(sum.mtx) == "patients$treat:1", colnames(sum.mtx) == "Pr(>|z|)"]
         pval.resp2[i] = sum.mtx[rownames(sum.mtx) == "patients$treat:2", colnames(sum.mtx) == "Pr(>|z|)"]
      }
   }

   pwr.resp.overall = sum(pval.resp < sig.overall)/runs
   pwr.resp2.overall = sum(pval.resp2 < sig.overall)/runs

   sens = list()
   pval.resp = matrix (nrow = nclust, ncol = runs)
   pval.resp2 = matrix (nrow = nclust, ncol = runs)
   estimate.rr = matrix (nrow = nclust, ncol = runs)
   estimate.rr2 = matrix (nrow = nclust, ncol = runs)
   for (i in 1:runs) {
      ## Find clusters
      sens[[i]] = sens.cvrs2(patients, covar[i,,], response[,i], response2[,i], seed, nclust)
      if (sens[[i]]$cluster.pred[1] == 0) {
         return (NA) 
      } else {
         for (j in 1:nclust) {
            ## Each cluster in turn is a sensitive group
            resp.group = response[sens[[i]]$cluster.pred == j, i] #response for sensitive group
            resp2.group = response2[sens[[i]]$cluster.pred == j, i] #response2 for sensitive group
            treat.group = patients$treat[sens[[i]]$cluster.pred == j] #treatment allocation for sensitive group

            conf.resp = matrix(nrow = 2, ncol = 2,
                data = c(sum(resp.group & !treat.group),  #number of responders in control
                      sum(!resp.group & !treat.group), #number of non-responders in control
                      sum(resp.group & treat.group),   #number of responders in treatment
                      sum(!resp.group & treat.group)), #number of non-responders in treatment
                byrow = TRUE)

            pval.resp[j, i] = fisher.test(conf.resp, alternative = "two.sided")$p.value
            estimate.rr[j, i] = conf.resp[2,1]/(conf.resp[2,1] + conf.resp[2,2])

            conf.resp2 = matrix(nrow = 2, ncol = 2,
                data = c(sum(resp2.group & !treat.group),  #number of responders in control
                      sum(!resp2.group & !treat.group), #number of non-responders in control
                      sum(resp2.group & treat.group),   #number of responders in treatment
                      sum(!resp2.group & treat.group)), #number of non-responders in treatment
                byrow = TRUE)

            pval.resp2[j, i] = fisher.test(conf.resp2, alternative = "two.sided")$p.value
            estimate.rr2[j, i] = conf.resp2[2,1]/(conf.resp2[2,1] + conf.resp2[2,2])
         } ##end loop on clusters
      } ## if cluster.pred isn't 0
   } ## end loop on runs

   cvrs = matrix(nrow = N, ncol = runs)
   cvrs2 = matrix(nrow = N, ncol = runs)
   rownames(cvrs) = rownames(patients)
   colnames(cvrs) = 1:ncol(response)
   sens.pred = matrix(nrow = N, ncol = runs)
   rownames(sens.pred) = rownames(patients)
   colnames(sens.pred) = 1:ncol(response)
   cluster.pred = matrix(nrow = N, ncol = runs)
   rownames(cluster.pred) = rownames(patients)
   colnames(cluster.pred) = 1:ncol(response)
   psens = matrix(nrow = nclust, ncol = runs)
   pspec = matrix(nrow = nclust, ncol = runs)
   for (i in seq(runs)) {
      psens[,i] = sens[[i]]$psens
      pspec[,i] =  sens[[i]]$pspec
      cvrs[,i] = sens[[i]]$cvrs
      cvrs2[,i] = sens[[i]]$cvrs2
      cluster.pred[,i] = sens[[i]]$cluster.pred
   }
   pwr.resp.group = numeric(length = nclust)
   pwr.resp2.group = numeric(length = nclust)
   for (i in 1:nclust) {
      pwr.resp.group[i] = sum(pval.resp[i,] < sig.group)/runs
      pwr.resp2.group[i] = sum(pval.resp2[i,] < sig.group)/runs
   }

   ## Overall power of the adaptive design
   pwr.resp.adaptive = pwr.resp.overall + (1-pwr.resp.overall)*pwr.resp.group
   pwr.resp2.adaptive = pwr.resp2.overall + (1-pwr.resp2.overall)*pwr.resp2.group

   output = list(patients = patients, pwr.resp.overall = pwr.resp.overall, pwr.resp2.overall = pwr.resp2.overall, pwr.resp.group = pwr.resp.group, pwr.resp2.group = pwr.resp2.group, pwr.resp.adaptive = pwr.resp.adaptive, pwr.resp2.adaptive = pwr.resp2.adaptive, estimate.rr = estimate.rr, estimate.rr2 = estimate.rr2, psens = psens, pspec = pspec, cvrs = cvrs, cvrs2 = cvrs2, cluster.pred = cluster.pred)

   if (plotrs) {
      plotcvrs = cvrs2.plot(cvrs, cvrs2, cluster.pred, patients$cluster.true, TRUE)
      ggsave ("cvrs2.png", plotcvrs)
   }
   class(output) = "rapids"
   return(output)
}

