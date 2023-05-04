#' @title 
#' Simulate data for the cross-validated adaptive signature design/cross-validated risc scores design.
#'
#' @description
#' The data is simulated assuming that the response to treatment is influenced by a subset of K unknown covariates (the sensitive covariates) through the following model:
#'
#' logit(p_i)= mu+lambda*t_i+gamma_1*t_i*x_i1+...+gamma_K*t_i*x_iK,
#'
#'where p_i is the probability of response to treatment for the i-th patient; mu is the intercept; lambda is the treatment main effect that all patients experience regardless of the values of the covariates; t_i is the treatment that the i-th patient receives (t_i = 0 for the control arm and t_i=1 for the treatment arm); x_i1,...,x_iK are the values for the K unknown sensitive covariates; gamma_1,...,gamma_K are treatment-covariate interaction effects for the K covariates.  
#'The model assumes that there is a subset of patients (the sensitive group) with a higher probability of response when treated with the new treatment, compared with the control treatment. 
#'
#'
#' @param N Number of patients.
#' @param L Overall number of covariates.
#' @param K Number of sensitive covariates.
#' @param mu1,sigma1,rho1 Mean, sd and correlation for sensitive covariates in sensitive patients.
#' @param mu2,sigma2,rho2 Correlation parameter for sensitive covariates in non-sensitive patients.
#' @param mu0,sigma0,rho Correlation parameter for non-sensitive covariates in all patients.
#' @param perc.sp  Percentage of sensitive patients.
#' @param rr.nsp.treat  Response rate on the treatment arm in non-sensitive patients.
#' @param rr.con  Response rate on the control arm.
#' @param rr.sp.treat  Response rate on the treatment arm in sensitive patients.
#' @param runs  Number of replicates to simulate.
#' @param seed  A seed for the random number generator.
#'
#' @return A list of 3 data frames: patients, covar, response.
#' @return patients: a data frame with one row per patient and the following columns: 
#'         FID (family ID), IID (individual ID), sens.true (true sensitivity indicator),
#'         treat (1 for treatment and 0 for control), rr (probability of response for a binary outcome)
#' @return covar: covariate data for L covariates
#' @return response: simulated binary responses, one column per simulation (number of columns = runs)
#' 
#' @examples
#' N = 400
#' L = 100
#' K=10
#' rho1 = 0
#' rho2 = 0
#' rho0 = 0
#' mu1 = 1
#' mu2 = 0
#' mu0 = 0
#' sigma1 = 0.5
#' sigma2 = 0.1
#' sigma0=0.5
#' perc.sp = 0.1
#' rr.nsp.treat = 0.25
#' rr.con = 0.25
#' rr.sp.treat = 0.98
#' runs = 10
#' seed = 123
#'
#' datalist = simulate.data (N , L , K, rho1, rho2, rho0, mu1, mu2, mu0, sigma1, sigma2, sigma0, perc.sp, rr.nsp.treat, rr.con, rr.sp.treat, runs, seed)
#' @seealso
#' \code{\link{analyse.simdata}} and \code{\link{cvrs.plot}} functions; \code{\link{print}} and \code{\link{plot}} methods.
#' @author Svetlana Cherlin, James Wason
#' @export simulate.data
#' @importFrom "MASS" "mvrnorm"

simulate.data <- function(N = 1000, L = 100, K=10, 
                     rho1 = 0, rho2 = 0, rho0 = 0,
                     mu1 = 1, mu2 = 0, mu0 = 0, sigma1 = 0.5,
                     sigma2 = 0.1, sigma0=0.5,
                     perc.sp = 0.1,
                     rr.nsp.treat = 0.25,
                     rr.con = 0.25,
                     rr.sp.treat = 0.6, 
		     runs = 1, seed = 123)
{

  mu = log(rr.con/(1-rr.con)) #intercept corresponding to response rate for controls
  lambda = log(rr.nsp.treat/(1-rr.nsp.treat)) - mu  #main treatment effect
  interaction.scaling = log(rr.sp.treat/(1-rr.sp.treat)) - mu - lambda #interaction scaling 

  #sensitive covariates, sensitive patients
  m1 = numeric(length = K) 
  m1[] = mu1
  Sigma1 = matrix(nrow = K, ncol = K, data = sigma1*sigma1*rho1) 
  diag(Sigma1) = sigma1^2 
  
  #sensitive covariates, non-sensitive patients
  m2 = numeric (length = K) 
  m2[] = mu2 
  Sigma2 = matrix(nrow = K, ncol = K, data = sigma2*sigma2*rho2) 
  diag(Sigma2) = sigma2^2 
  
  #non-sensitive covariates, all patients
  m0 = numeric(length = L-K) 
  m0[0] = mu0 
  Sigma0 = matrix(nrow = L-K, ncol = L-K, data = sigma0*sigma0*rho0)
  diag(Sigma0) = sigma0^2 
  
  ## Simulate covariates for sensitive patients
  scovar.sp = mvrnorm(n = N*perc.sp, m1, Sigma1, tol = 1e-6) #sensitive covariates
  nscovar.sp = mvrnorm(n = N*perc.sp, m0, Sigma0, tol = 1e-6) #non-sensitive covariates
  sp = cbind(scovar.sp, nscovar.sp) 
  sp = as.data.frame(sp)
  sp$sens.true = 1 
  
  ## Simulate covariates for non-sensitive patients
  scovar.nsp = mvrnorm(n = N*(1-perc.sp), m2, Sigma2, tol = 1e-6) #sensitiv covariates
  nscovar.nsp = mvrnorm(n = N*(1-perc.sp), m0, Sigma0, tol = 1e-6) #non-sensitive covariates
  nsp = cbind(scovar.nsp, nscovar.nsp) 
  nsp = as.data.frame(nsp)
  nsp$sens.true = 0 
  
  ## Equal randomisation to control/treatment arm for sensitive patients
  ind = 1:nrow(sp)
  control.sp= sp[ind %% 2 == 0,]
  treatment.sp = sp[ind %% 2 == 1,]
  
  ## Equal randomisation to control/treatment arm for non-sensitive patients
  ind = 1:nrow(nsp)
  control.nsp = nsp[ind %% 2 == 0,]
  treatment.nsp = nsp[ind %% 2 == 1,]
  
  ## Combine control and treatment
  control = rbind(control.sp, control.nsp)
  treatment = rbind(treatment.sp, treatment.nsp)
  control$treat = 0   #control/treatment indicator
  treatment$treat = 1 #control/treatment indicator
  patients.data  = rbind (control, treatment)
  rownames(patients.data) = seq(1:N)
  
  covar = as.data.frame(patients.data[, 1:L])
  colnames(covar) = paste ("Covar", 1:L, sep = "")
  patients = data.frame(FID = seq(1:N), IID = seq(1:N), sens.true = patients.data$sens.true, treat =  patients.data$treat)

  # Simulate response probabilities for binomail outcoume or rate for poisson outcome,  based on sensitive covariates
  sens.covar = as.matrix(covar[,1:K])
  gamma = sapply(m1, function(x) { ifelse (x==0, 0, interaction.scaling/(K*x))}) ##covariate-treatment interation
  sens.covar.levels = sens.covar %*% as.matrix(gamma, nrow = K)
  linpred =  mu + patients$treat*lambda + patients$treat*sens.covar.levels
  patients$rr = as.vector(exp(linpred)/(1+exp(linpred)))
  response = matrix(nrow = N, ncol = runs, 
     data = rbinom(N*runs, 1, rep(patients$rr, runs)), byrow = FALSE)
  rownames(response) = seq(1:N)
  colnames(response) = paste("resp", 1:runs, sep = "")

  return(list(patients = patients, covar = covar, response = response))  
}


#' @title 
#' Simulate data for the cross-validated risc scores design with 2 outcomes (CVRS2).
#'
#' @description
#' The data is simulated assuming that there are two outcomes (response and response2). Both response and response2 are influenced by a subset of K and T unknown covariates (the sensitive covariates) through the following model:
#'
#' logit(p.response_i)= mu+lambda*t_i+gamma_1*t_i*x_i1+...+gamma_K*t_i*x_iK,
#' logit(p.response2_i)= mu+lambda*t_i+gamma_1*t_i*x_i1+...+gamma_K*t_i*x_iT,
#'
#'where p.response_i is the probability of response to treatment for the i-th patient; mu is the intercept; lambda is the treatment main effect that all patients experience regardless of the values of the covariates; t_i is the treatment that the i-th patient receives (t_i = 0 for the control arm and t_i=1 for the treatment arm); x_i1,...,x_iK/x_iT are the values for the K/T unknown sensitive covariates; gamma_1,...,gamma_K/gamma_T are treatment-covariate interaction effects for the K/T covariates. The model assumes that there is a subset of patients (the resp.sensitive group) with a higher probability of response when treated with the new treatment, compared with the control treatment, and a subset of patients (the resp2.sensitive group) with a highter probability of response2 when treated with the new treatment.  The model assumes that there are 4 clusters of patients: cluster1 (high probability of both response and response2), cluster2 (low probability of response and hight probability of response2), cluster3 (high probability of response and low probability of response2) and cluster4(low probability of both response and response2). Cluster1 is considered a sensitive group (high probability of both  response and response2  when treated with the new treatment, compared with the control treatment. 
#'
#' @param N Number of patients.
#' @param L Overall number of covariates.
#' @param K Number of sensitive covariates that influence the response only.
#' @param K2 Number of sensitive covariates that influence the response2 only.
#' @param Both Number of overlapping sensitive covariates (influence both the response and  response2).
#' @param mu1,sigma1,rho1 Mean, sd and correlation for sensitive covariates in sensitive patients.
#' @param mu2,sigma2,rho2 Correlation parameter for sensitive covariates in non-sensitive patients.
#' @param mu0,sigma0,rho Correlation parameter for non-sensitive covariates in all patients.
#' @param perc.sp  Percentage of patients with higher probability of response.
#' @param rr.nsp.treat  Response rate on the treatment arm in non-resp.sensitive patients.
#' @param rr.con  Response rate on the control arm.
#' @param rr.sp.treat  Response rate on the treatment arm in resp.sensitive patients.
#' @param mu1.resp2,sigma1.resp2,rho1.resp2 Mean, sd and correlation for covariates that influence the response2
#' @param mu1.both,sigma1.both,rho1.both Mean, sd and correlation for overlaping covariates in sensitive patients.
#' @param mu2.both,sigma2.both,rho2.both Mean, sd and correlation for overlaping covariates in non-sensitive patients.
#' @param rr2.con Probability of response2 for the control arm
#' @param rr2.nsp.treat Probability of response2 for the non-resp2.sensitive patients in the treatment arm
#' @param rr2.sp.treat Probability of response2 for the resp2.sensitive patients in the treatment arm
#' @param runs  Number of replicates to simulate.
#' @param seed  A seed for the random number generator.
#' @return A list of 4 data frames: patients, covar, response, response2
#' @return patients: a data frame with one row per patient and the following columns: 
#'         FID (family ID), IID (individual ID), sens.resp.true (true sensitivity to response),
#'         sens.resp2.true (true sensitivity to response2), cluster.true (1 for sens.pred.true ==1 and sens.pred2.true == 1, 2 for sens.pred.true ==0 and sens.pred2.true == 1, 3 for sens.pred.true ==1 and sens.pred2.true == 0 and 4 for sens.pred.true ==0 and sens.pred2.true == 0), sens.true (1 for cluster.true == 1, 0 otherwise).
#'         treat (1 for treatment and 0 for control), rr (probability of response), rr2 (probability of response2)
#' @return covar: covariate data for L covariates
#' @return response: simulated binary variable "response", one column per simulation (number of columns = runs)
#' @return response2: simulated binary variable "response2", one column per simulation (number of columns = runs)
#' 
#' @examples
#' N = 1000
#' L = 100
#' K1 = 10
#' K2 = 10
#' Both = 3
#' rho1 = 0
#' rho2 = 0
#' rho0 = 0
#' mu1 = 1
#' mu2 = 0
#' mu0 = 0
#' sigma1 = 0.5
#' sigma2 = 0.1
#' sigma0 = 0.5
#' rho1.resp2 = 0
#' rho2.resp2 = 0
#' mu1.resp2 = 1
#' mu2.resp2 = 0
#' sigma1.resp2 = 0.5
#' sigma2.resp2 = 0.1
#' mu1.both = 1
#' mu2.both = 0
#' sigma1.both = 0.5
#' sigma2.both = 0.1
#' rho1.both = 0
#' rho2.both = 0
#' perc.sp = 0.1
#' perc.sp2 = 0.1
#' perc.sp.both = 0.1
#' rr.con = 0.25
#' rr.sp.treat = 0.8
#' rr2.sp.treat = 0.8
#' rr.nsp.treat = 0.25
#' rr2.nsp.treat = 0.25
#' rr2.con = 0.25
#' runs = 5
#' seed = 123
#' datalist2 = simulate.data2 (N , L , K1, K2, Both, rho1, rho2, rho0, mu1, mu2, mu0, sigma1, sigma2, sigma0, rho1.resp2, rho2.resp2, mu1.resp2, mu2.resp2, sigma1.resp2, sigma2.resp2,mu1.both, mu2.both, sigma1.both, sigma2.both, rho1.both, rho2.both, perc.sp, perc.sp2, perc.sp.both, rr.nsp.treat, rr.con, rr.sp.treat, rr2.nsp.treat, rr2.con, rr2.sp.treat, runs, seed)

#' @seealso
#' \code{\link{analyse.simdata2}} and \code{\link{cvrs2.plot}} functions; \code{\link{print}} and \code{\link{plot}} methods.
#' @author Svetlana Cherlin, James Wason
#' @importFrom "MASS" "mvrnorm"
#' @export simulate.data2

simulate.data2 <- function(N = 1000, L = 100, 
                     K1 = 10, K2 = 10, Both = 3, 
                     rho1 = 0, rho2 = 0, rho0 = 0,
                     mu1 = 1, mu2 = 0, mu0 = 0,
                     sigma1 = 0.5, sigma2 = 0.1, sigma0 = 0.5,
                     rho1.resp2 = 0, rho2.resp2 = 0,
                     mu1.resp2 = 1, mu2.resp2 = 0,
                     sigma1.resp2 = 0.5, sigma2.resp2 = 0.1,
		                 mu1.both = 1, mu2.both = 0, 
		                 sigma1.both = 0.5, sigma2.both = 0.1, 
		                 rho1.both = 0, rho2.both = 0, 
                     perc.sp = 0.1,
		                 perc.sp2 = 0.1,
		                 perc.sp.both = 0.1,
                     rr.nsp.treat = 0.25,
                     rr.con = 0.25,
                     rr.sp.treat = 0.8,
                     rr2.nsp.treat = 0.25,
                     rr2.con = 0.25,
                     rr2.sp.treat = 0.8,
                     runs = 1, seed = 123)
{
  covar = array(NA, c(runs, N, L))
  ## Response parameters
  mu = log(rr.con/(1-rr.con)) # intercept corresponding to response rate for controls
  lambda = log(rr.nsp.treat/(1-rr.nsp.treat)) - mu  # main treatment effect
  interaction.scaling = log(rr.sp.treat/(1-rr.sp.treat)) - mu - lambda #interaction scaling 

  ## Response2 parameters
  mu.resp2 = log(rr2.con/(1-rr2.con)) #intercept corresponding to the response2 rate for controls
  lambda.resp2 = log(rr2.nsp.treat/(1-rr2.nsp.treat)) - mu.resp2  #main treatment effect
  interaction.scaling.resp2 = log(rr2.sp.treat/(1-rr2.sp.treat)) - mu.resp2 - lambda.resp2 #interaction scaling 

  # resp.sensitive covariates, resp.sensitive patients
  m1 = numeric(length = K1)
  m1[] = mu1
  Sigma1 = matrix(nrow = K1, ncol = K1, data = sigma1*sigma1*rho1)
  diag(Sigma1) = sigma1^2

  # resp2.sensitive covariates, resp2.sensitive patients
  m1.resp2 = numeric(length = K2)
  m1.resp2[] = mu1.resp2
  Sigma1.resp2 = matrix(nrow = K2, ncol = K2, data = sigma1.resp2*sigma1.resp2*rho1.resp2)
  diag(Sigma1.resp2) = sigma1.resp2^2

  # resp.sensitive covariates, non-resp.sensitive patients
  m2 = numeric (length = K1)
  m2[] = mu2
  Sigma2 = matrix(nrow = K1, ncol = K1, data = sigma2*sigma2*rho2)
  diag(Sigma2) = sigma2^2

  # resp2.sens covariates, non-resp2.sens patients
  m2.resp2 = numeric (length = K2)
  m2.resp2[] = mu2.resp2
  Sigma2.resp2 = matrix(nrow = K2, ncol = K2, data = sigma2.resp2*sigma2.resp2*rho2.resp2)
  diag(Sigma2.resp2) = sigma2.resp2^2

  # overlapping sensitive covariates, sensitive patients
  m1.both = numeric (length = Both)
  m1.both[] = mu1.both
  Sigma1.both = matrix(nrow = Both, ncol = Both, data = sigma1.both*sigma1.both*rho1.both)
  diag(Sigma1.both) = sigma1.both^2

  # overlapping sensitive covariates, non-sensitive patients
  m2.both = numeric (length = Both)
  m2.both[] = mu2.both
  Sigma2.both = matrix(nrow = Both, ncol = Both, data = sigma2.both*sigma2.both*rho2.both)
  diag(Sigma2.both) = sigma2.both^2

  # non-sens covariates, all patients
  K0 = L-K1-K2-Both
  m0 = numeric(length = K0)
  m0[0] = mu0
  Sigma0 = matrix(nrow = K0, ncol = K0, data = sigma0*sigma0*rho0)
  diag(Sigma0) = sigma0^2

  rr = matrix(data = NA, nrow = N, ncol = runs)
  rr2 = matrix(data = NA, nrow = N, ncol = runs)
  response = matrix(data = NA, nrow = N, ncol = runs)
  response2 = matrix(data = NA, nrow = N, ncol = runs)

  for (i in 1:runs) { 

     sens.resp.true = numeric(length = N)
     sens.resp2.true = numeric(length = N)

     ## Response sensitive patients 
     
     if (K1 == 0) {
	      resp.covar = matrix(nrow = N*perc.sp, ncol = 0)
     } else {
        resp.covar = mvrnorm(n = N*perc.sp, m1, Sigma1, tol = 1e-6)
     }
     if (K2 == 0) {
        resp2.covar = matrix(nrow = N*perc.sp, ncol = 0)
     } else {
        resp2.covar = mvrnorm(n = N*perc.sp, m2.resp2, Sigma2.resp2, tol = 1e-6)
      }
     if (Both == 0) {
	      both.covar = matrix(nrow = N*perc.sp, ncol = 0)
     } else {
        both.covar = mvrnorm(n = N*perc.sp, m1.both, Sigma1.both, tol = 1e-6)
     }
     if (K0 == 0) {
	      other.covar = matrix(n = N*perc.sp, ncol = 0)
     } else {
        other.covar = mvrnorm(n = N*perc.sp, m0, Sigma0, tol = 1e-6)
     }

     resp.sens.patients = cbind(resp.covar, resp2.covar, both.covar, other.covar)
     sens.resp.true[1:(N*perc.sp)] = 1

     ## Response2 sensitive patients

     if (K1 == 0) {
        resp.covar = matrix(nrow = N*perc.sp2, ncol = 0)
     } else {        
        resp.covar = mvrnorm(n = N*perc.sp2, m2, Sigma2, tol = 1e-6)
     }
     if (K2 ==0) {
	      resp2.covar = matrix(nrow = N*perc.sp2, ncol = 0)
     } else {        
        resp2.covar = mvrnorm(n = N*perc.sp2, m1.resp2, Sigma1.resp2, tol = 1e-6)
     }
     if (Both == 0) {
        both.covar = matrix(nrow = N*perc.sp2, ncol = 0)
     } else {
        both.covar = mvrnorm(n = N*perc.sp2, m1.both, Sigma1.both, tol = 1e-6)
     }
     if (K0 == 0) {
        other.covar = matrix(n = N*perc.sp2, ncol = 0)
     } else {        
        other.covar = mvrnorm(n = N*perc.sp2, m0, Sigma0, tol = 1e-6)
     }
     
     resp2.sens.patients = cbind(resp.covar, resp2.covar, both.covar, other.covar)
     sens.resp2.true[(N*perc.sp+1):(N*perc.sp+N*perc.sp2)] = 1

     ## Patients sensitive to both responses

     if (K1 == 0) {
	      resp.covar = matrix(nrow = N*perc.sp.both, ncol = 0)
     } else {
        resp.covar = mvrnorm(n = N*perc.sp.both, m1, Sigma1, tol = 1e-6)
     }
     if (K2 == 0) {
	      resp2.covar = matrix(nrow = N*perc.sp.both, ncol = 0)
     } else {
        resp2.covar = mvrnorm(n = N*perc.sp.both, m1.resp2, Sigma1.resp2, tol = 1e-6)
     }
     if (Both == 0) {
        both.covar = matrix(nrow = N*perc.sp.both, ncol = 0)
     } else {
        both.covar = mvrnorm(n = N*perc.sp.both, m1.both, Sigma1.both, tol = 1e-6)
     }
     if (K0 == 0) {
	      other.covar = matrix(nrow = N*perc.sp.both, ncol = 0)
     } else {
        other.covar = mvrnorm(n = N*perc.sp.both, m0, Sigma0, tol = 1e-6)
     }

     both.sens.patients = cbind(resp.covar, resp2.covar, both.covar, other.covar)
     sens.resp.true[(N*perc.sp+N*perc.sp2+1):(N*perc.sp+N*perc.sp2+N*perc.sp.both)] = 1
     sens.resp2.true[(N*perc.sp+N*perc.sp2+1):(N*perc.sp+N*perc.sp2+N*perc.sp.both)] = 1

     ## Non-sensitive patients

     NN = N*(1 - perc.sp - perc.sp2 - perc.sp.both)
     if (K1 == 0) {
	      resp.covar = matrix(nrow = NN, ncol = 0)
     } else {
        resp.covar = mvrnorm(n = NN, m2, Sigma2, tol = 1e-6)
     }
     if (K2 == 0) {
	      resp2.covar = matrix(nrow = NN, ncol = 0)
     } else {
        resp2.covar = mvrnorm(n = NN, m2.resp2, Sigma2.resp2, tol = 1e-6)
     }
     if (Both == 0) {
        both.covar = matrix(nrow = NN, ncol = 0) 
     } else {
        both.covar = mvrnorm(n = NN, m2.both, Sigma2.both, tol = 1e-6)
     }
     if (K0 == 0) {
	      other.covar = matrix(nrow = NN, ncol = 0)
     } else {
        other.covar = mvrnorm(n = NN, m0, Sigma0, tol = 1e-6)
     }

     non.sens.patients = cbind(resp.covar, resp2.covar, both.covar, other.covar)

     ## Combine patients

     covar.data = as.data.frame(rbind(resp.sens.patients, resp2.sens.patients, both.sens.patients, non.sens.patients))
     colnames(covar.data) = paste ("Covar", 1:L, sep = "")
     covar.data = as.matrix(covar.data)

     ## Gamma

     gamma = sapply(c(m1,m1.both), function(x) { ifelse (x==0, 0, interaction.scaling/((K1+Both)*x))}) 
     gamma2 = sapply(c(m1.resp2, m1.both), function(x) { ifelse (x==0, 0, interaction.scaling.resp2/((K2+Both)*x))}) 

     ## Equal randomisation to control/treatment arm 

     treat = seq(N) %% 2 #alternates between 1 and 0, N times

     # Rate of response for all patients

     if (Both != 0) {
        if (K1 == 0) {
           resp.sens.covar.levels = covar.data[,(K2+1):(K2+Both)] %*% as.matrix(gamma, nrow = (Both))
        } else {     
           resp.sens.covar.levels = covar.data[,c(1:K1, (K1+K2+1):(K1+K2+Both))] %*% as.matrix(gamma, nrow = (K1+Both))
        }
     } else { #Both == 0 therefore K1 cannot be 0 
           resp.sens.covar.levels = covar.data[,1:K1] %*% as.matrix(gamma, nrow = (K1))
     }

     #Rate of response for response2-sens. patients
     if (K1 > 0) {
        resp.sens.covar.levels[(N*perc.sp+1):(N*perc.sp+N*perc.sp2)] = covar.data[(N*perc.sp+1):(N*perc.sp+N*perc.sp2),1:K1] %*% as.matrix(gamma[1:K1], nrow = K1) 
     }

     #Rate of response2 for all patients

     resp2.sens.covar.levels = covar.data[, (K1+1):(K1+K2+Both)] %*% as.matrix(gamma2, nrow = (K2+Both))

     # Rate of response2 for response-sens. patients
     if (K2 > 0) {
        resp2.sens.covar.levels[1:(N*perc.sp)] = covar.data[1:(N*perc.sp), (K1+1):(K1+K2)] %*% as.matrix(gamma2[1:K2], nrow = K2)
     }

     ## Response rate for treatment

     resp.linpred =  mu + treat*lambda + treat*resp.sens.covar.levels
     rr[,i] = as.vector(exp(resp.linpred)/(1+exp(resp.linpred)))

     resp2.linpred =  mu.resp2 + treat*lambda.resp2 + treat*resp2.sens.covar.levels
     rr2[,i] = as.vector(exp(resp2.linpred)/(1+exp(resp2.linpred)))

     ## Simulate response and response2

     response[,i] = rbinom(N, 1, rr[,i])
     response2[,i] = rbinom(N, 1, rr2[,i])
     
     ## Treatment/control sort (Treatment first)
     covar.data = as.matrix(covar.data[order(treat, decreasing = TRUE),])
     rownames(covar.data) = seq(nrow(covar.data))

     rr[,i] = rr[order(treat, decreasing = TRUE),i]
     rr2[,i] = rr2[order(treat, decreasing = TRUE),i]

     response[,i] = response[order(treat, decreasing = TRUE),i]
     response2[,i] = response2[order(treat, decreasing = TRUE),i]

     sens.resp.true = sens.resp.true[order(treat, decreasing = TRUE)]
     sens.resp2.true = sens.resp2.true[order(treat, decreasing = TRUE)]

     treat = treat[order(treat, decreasing = TRUE)]

     covar[i,,] = covar.data

  } ## End of loop on the number of replicates

  patients = data.frame(FID = seq(1:N), IID = seq(1:N), treat =  treat, sens.resp.true = sens.resp.true, sens.resp2.true = sens.resp2.true)
  rownames(rr) = seq(1:N)
  rownames(rr2) = seq(1:N)
  rownames(response) = seq(1:N)
  colnames(response) = paste("resp", 1:runs, sep = "")
  rownames(response2) = seq(1:N)
  colnames(response2) = paste("resp2", 1:runs, sep = "")

  ## Define clusters

  patients$cluster.true = 0
  patients$cluster.true[patients$sens.resp.true == 1 & patients$sens.resp2.true == 1] = 4 #2 ###For 2 clusters
  patients$cluster.true[patients$sens.resp.true == 0 & patients$sens.resp2.true == 1] = 2 #1 ###For 2 clusters
  patients$cluster.true[patients$sens.resp.true == 1 & patients$sens.resp2.true == 0] = 3 #1 ###For 2 clusters
  patients$cluster.true[patients$sens.resp.true == 0 & patients$sens.resp2.true == 0] = 1 #1 ###For 2 clusters

  return(list(patients = patients, covar = covar, response = response, response2 = response2, rr = rr, rr2 = rr2))
}

