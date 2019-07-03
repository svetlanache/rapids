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
