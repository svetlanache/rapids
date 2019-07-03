# rapids
A package for cross-validated adaptive signature design (CVASD)
and cross-validation risk scores design (CVRS).<br/>
The package provides three main funcitons:
‘simulate.data’ ‘analyse.simdata’ and ‘analyse.realdata’.
Additional functions are ‘permutation.test’ for permutation tests
for the real data, ‘cvrs.plot’ for plotting the risk scores, and
also ‘print’ and ‘plot’ generic methods.<br/>
simulate.data:<br/>
     The function simulates covariates data and binary responses to be
     used in the analysis of the cross-validated adaptive signature
     design or cross-validation risk scores design.<br/>
analyse.simdata:<br/>
     The function computes the power of the design for the simulated
     data according to the input method ("cvasd" or "cvrs").<br/>
analyse.realdata:<br/>
     The function computes the p-value for the interaction effect
     between the treatment and the sensitivity status.  The sensitivity
     status is predicted according to the input method ("cvasd" and
     "cvrs").<br/>
permutation.test:<br/>
     The function performs permutation test for the real data.<br/>
cvrs.plot:<br/>
     The function plots the risk scores for the "cvrs" method.
# installation
library('devtools')<br/>
install_github('svetlanache/rapids')
# example
```{r }
library(rapids)

## Simulate data
N = 400
L = 100
K=10
rho1 = 0
rho2 = 0
rho0 = 0
mu1 = 1
mu2 = 0
mu0 = 0
sigma1 = 0.5
sigma2 = 0.1
sigma0 = 0.5
perc.sp = 0.1
rr.nsp.treat = 0.25
rr.con = 0.25
rr.sp.treat = 0.98
runs = 10
seed = 123
     
simdata = simulate.data (N , L , K, rho1, rho2, rho0, mu1, mu2, mu0, sigma1, sigma2, sigma0, perc.sp, rr.nsp.treat, rr.con, rr.sp.treat, runs, seed)

## Analyse simulate data with the "cvrs" methods
sig = 0.05
group.prop.sig = 0.2
method = "cvrs"
seed = 123
plotrs = T
eta = NULL
R = NULL
G = NULL

simres.cvrs = analyse.simdata (simdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)

## Analyse simulated data with the "cvasd" method
sig = 0.05
group.prop.sig = 0.2
method = "cvasd"
seed = 123
plotrs = T
eta = c(0.01, 0.02, 0.03)
R = c(2.5, 2, 1.5)
G = c(3,2,1)
simres.cvasd = analyse.simdata (simdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)
 
```
   
