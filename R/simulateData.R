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
