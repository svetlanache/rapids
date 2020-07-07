

#' @title Get subgroup of sensitive patiens according to the "cvasd" method.
#'
#' @description
#'
#' For a vector of binary responses, get sensitive patients by CVASD
#' according to the input tuning set (eta.in, R.in, G.in).      
#' In the training subset, fit logistic regression for each covariate.          
#' In the testing subset, compute the new vs control arm odds ratio for 
#' covariates that have treatment-covariate interaction significant at a level eta.in
#' in the training subset. A patient is sensitive if the odds ration > R.in for at least G.in covariates.
#' If the length of each parameter in the tuning set (eta.in, R.in and G.in) > 1 then the tuning
#' set is found by a nested cross-validation where only one inner fold is used (get.tuning.param function). 
#' 
#' @param  patients - a data frame of patients inormation 
#'         covar - a data frame of covariates
#'         y - a vector of responses 
#'         eta.in - significance level for covariate-wise logistic regression (double or a vector of doubles)
#'         R.in - a threshold of the odds ratio (double or a vector of doubles)
#'         G.in - a threshold for the number of covariates (integer or a vector of integers)
#'         seed - a seed for random number generator  
#'
#' @return  A list of 3 :
#'          psens - sensitivity of identifying the sensitive group, one value per simulation run 
#'          pspec - specificity of identifying the sensitive group, one value per simulation run 
#'          sens.pred - predicted sensitivity status (rows = patienst, columns = simulations)
#'
#' @author Svetlana Cherlin, James Wason

sens.cvasd = function(patients, covar, y, eta.in, R.in, G.in, seed)
{
        if (is.null(R.in) | is.null(eta.in) | is.null(G.in)) {
           stop ("Parameters eta, R, and G should be specified.")
        }
        if (!(length(R.in) == length(eta.in) & length(eta.in) == length(G.in))) {
           stop("Parameters eta, R, G should have the same length")
        }        
        if (!is.null(seed)) {
           set.seed(seed)
        } 

        ## Divide to folds amd keep prevalence of responders/non-responders within folds
        nfolds = 10       
        patients.nr = patients[y==0,]
        patients.r = patients[y==1,]
        covar.nr = covar[y==0,]
        covar.r = covar[y==1,]
        y.nr = y[y==0]
        y.r = y[y==1]
        foldid.nr = sample(rep(seq(nfolds), length = length(y.nr))) #non-responders
        foldid.r = sample(rep(seq(nfolds), length = length(y.r))) #responders     
        sens.df = data.frame()
        
        ## CV loop
        for (i in seq(nfolds)) 
        {
           which.nr = foldid.nr == i
           which.r = foldid.r == i
           patients.train = rbind(patients.nr[!which.nr,], patients.r[!which.r,])
           patients.test = rbind(patients.nr[which.nr,], patients.r[which.r,])
           covar.train = rbind(covar.nr[!which.nr,], covar.r[!which.r,])
           covar.test = rbind(covar.nr[which.nr,], covar.r[which.r,])
           y.train = c(y.nr[!which.nr], y.r[!which.r])
           y.test = c(y.nr[which.nr], y.r[which.r])
   
           ## Fit single covariate regression model for the training data
   	   reg.res = apply (as.matrix(covar.train), 2, function (x)
           {

            if (with(patients, exists('sens.true'))) {
               mod = glm(y.train ~ patients.train$treat + patients.train$treat:x, family = "binomial")  #for sim data
	    } else {
               mod = glm(y.train ~ patients.train$treat + x + patients.train$treat:x, family = "binomial")  #for real data
	    }
           
            if (is.na(mod$coeff[names(mod$coeff) == "patients.train$treat:x"])) {
                     pval = 1
                     lambda.hat = 0
                     beta.hat = 0
             } else {
                     mod.sum = summary(mod)$coefficients
                     pval = mod.sum[rownames(mod.sum) == "patients.train$treat:x",colnames(mod.sum) == "Pr(>|z|)"]
                     lambda.hat = mod.sum[rownames(mod.sum) == "patients.train$treat",colnames(mod.sum) == "Estimate"]
                     beta.hat = mod.sum[rownames(mod.sum) == "patients.train$treat:x",colnames(mod.sum) == "Estimate"]
             }
    	     cbind (pval, lambda.hat, beta.hat)
	   })

	   reg.res = t(as.matrix(reg.res))
	   colnames(reg.res) = c("pval", "lambda.hat", "beta.hat")
	   reg.res = as.data.frame(reg.res)       

           eta = eta.in
           R = R.in
           G = G.in

           ## Get tuning parameters using inner CV
           if (length(R.in) > 1) {
              tuning.param = get.tuning.param(patients.train, covar.train, y.train, eta.in, R.in, G.in, seed)
              eta = tuning.param$eta
              R = tuning.param$R
              G = tuning.param$G
           } 
	   ## Find sensitivity covar from training data
	   covar.sig = rownames(reg.res[reg.res$pval < eta,]) #names of sensitivity covariates
	   ind = match(covar.sig, colnames(covar.train))
           
           ## Compute odds ratio for the testing data
 	   odds.ratio = apply(as.matrix(covar.test), 1, function (x) 
           {
  	      exp(reg.res$lambda.hat[ind] +
                  reg.res$beta.hat[ind]*x[ind])
           })
	   odds.ratio = t(as.matrix(odds.ratio))         
	   
           ## Predict sensitivity status for the test data
	   sens.pred = apply(odds.ratio, 1, function(x)
           {
  	      length(x[x > R]) > G 
           })
	   sens.fold = data.frame(FID = patients.test$FID, IID = patients.test$IID, sens.pred = sens.pred)  
           sens.df = rbind(sens.df, sens.fold)   
        } ## End of the CV loop

        ## Order according to the patients data
        m = match (paste(as.character(patients$FID), as.character(patients$IID), sep = ":"), 
                   paste(as.character(sens.df$FID), as.character(sens.df$IID), sep = ":"))
        sens.df = sens.df[m,]
       

        ## Compute sensitivity and specificity of the sensitive group selection algorithm for simulated data
        if (with(patients, exists('sens.true'))) {
           conf = matrix(nrow = 2, ncol = 2,
              data = c(sum(!sens.df$sens.pred & !patients$sens.true), #predicted non.sens, true non.sens [1,1]
                   sum(!sens.df$sens.pred & patients$sens.true),      #predicted non.sens, true sens [1,2]
                   sum(sens.df$sens.pred & !patients$sens.true),      #predicted sens, true non. sens [2,1]
                   sum(sens.df$sens.pred & patients$sens.true)),      #predicted sens, true sens [2,2]
                   byrow = TRUE)
           psens = conf[2,2]/(conf[2,2] + conf[1,2])
           pspec = conf[1,1]/(conf[1,1] + conf[2,1])
           ret = list(psens = psens, pspec = pspec, sens.pred = sens.df$sens.pred, cvrs = sens.df$cvrs)
        } else {ret  = list(sens.pred = sens.df$sens.pred, cvrs = sens.df$cvrs) }

        return (ret)


}

#' @title Get tuning parameters for CVASD by nested cross-validation.
#'
#' @description
#' An internal function for tuning the parameters.
#' In one inner fold only, compute p-value for the test of treatment effect comparison in the 
#' selected sensitive group, for each combination of (eta, R, G). Output the (eta, R, G) combination that
#' corresponds to the smallest p-value.
#'           
#' @param  patients - a data frame of patients inormation 
#'         covar - a data frame of covariates
#'         y - a vector of responses 
#'         eta - significance level for covariate-wise logistic regression for "cvasd" method (double of vector of doubles)
#'         R - a threshold of the odds ratio for "cvasd" method (double of vector of doubles)
#'         G - a threshold for the number of covar for "cvasd" method (integer of vector of integers)
#'         seed - a seed for random number generator  
#'
#' @return A data frame with the optimal eta, R, and G 
#'
#' @author Svetlana Cherlin, James Wason

get.tuning.param = function(patients.train, covar.train, y.train, eta, R, G, seed)
{
          if (!is.null(seed)) {
             set.seed(seed)
          }
          nfolds=10
          M = data.frame(eta, R, G)
          p.value = numeric(length = nrow(M))

          ## Loop for eta, R, G combination
          for (j in 1:nrow(M))
          {   
             ## Divide to inner folds, keep prevalence of responsers/non-responders within the inner folds
             patients.inner.nr = patients.train[y.train==0,]
             patients.inner.r = patients.train[y.train==1,]
     
             covar.inner.nr = covar.train[y.train==0,]
             covar.inner.r = covar.train[y.train==1,]
     
             y.inner.nr = y.train[y.train==0]
             y.inner.r = y.train[y.train==1]

             foldid.inner.nr = sample(rep(seq(nfolds), length = length(y.inner.nr))) #non-responders
             foldid.inner.r = sample(rep(seq(nfolds), length = length(y.inner.r)))   #responders

             ## Use first fold only in the inner CV
             which.inner.nr = foldid.inner.nr == 1
             which.inner.r = foldid.inner.r == 1

             patients.train.inner = rbind(patients.inner.nr[!which.inner.nr,], patients.inner.r[!which.inner.r,])
             patients.test.inner = rbind(patients.inner.nr[which.inner.nr,], patients.inner.r[which.inner.r,])

             covar.train.inner = rbind(covar.inner.nr[!which.inner.nr,], covar.inner.r[!which.inner.r,])
             covar.test.inner = rbind(covar.inner.nr[which.inner.nr,], covar.inner.r[which.inner.r,])

             y.train.inner = c(y.inner.nr[!which.inner.nr], y.inner.r[!which.inner.r])
             y.test.inner = c(y.inner.nr[which.inner.nr], y.inner.r[which.inner.r])

             ## Fit single-covariate regression model to the train data
     	     reg.res.inner = apply (as.matrix(covar.train.inner), 2, function (x)
             {
                  if (with(patients.train, exists('sens.true'))) {
                     mod = glm(y.train.inner ~ patients.train.inner$treat + patients.train.inner$treat:x, family = "binomial") #for sim data
	          } else {
                     mod = glm(y.train.inner ~ patients.train.inner$treat + x + patients.train.inner$treat:x, family = "binomial") #for real data
		  }
                  if (is.na(mod$coeff[names(mod$coeff) == "patients.train.inner$treat:x"])) {
                     pval = 1
                     lambda.hat = 0
                     beta.hat = 0
                  } else {
                     mod.sum = summary(mod)$coefficients
                     pval = mod.sum[rownames(mod.sum) == "patients.train.inner$treat:x",colnames(mod.sum) == "Pr(>|z|)"]
                     lambda.hat = mod.sum[rownames(mod.sum) == "patients.train.inner$treat",colnames(mod.sum) == "Estimate"]
                     beta.hat = mod.sum[rownames(mod.sum) == "patients.train.inner$treat:x",colnames(mod.sum) == "Estimate"]
                  }
    	          cbind (pval, lambda.hat, beta.hat)
	     })
	     reg.res.inner = t(as.matrix(reg.res.inner))
	     colnames(reg.res.inner) = c("pval", "lambda.hat", "beta.hat")
	     reg.res.inner = as.data.frame(reg.res.inner)
             
             ## Get sensitive covar from the train data
	     covar.sig.inner = rownames(reg.res.inner[reg.res.inner$pval < M$eta[j],]) #names of sensitivity covar
	     col.ind.inner = match(covar.sig.inner, colnames(patients.train.inner))
             row.ind.inner = match(covar.sig.inner, rownames(reg.res.inner))                   

   	     ## Compute odds ratio for sensitive covar in the test data
	     odds.ratio.inner = apply(as.matrix(patients.test.inner), 1, function (x) 
             {
  	          exp(reg.res.inner$lambda.hat[row.ind.inner] + 
                          reg.res.inner$beta.hat[row.ind.inner]*x[col.ind.inner])
             })
	     odds.ratio.inner = t(as.matrix(odds.ratio.inner))        

	     ## Predict sensitivity status in the test data
	     sens.pred = apply(odds.ratio.inner, 1, function(x)
             {
  	        length(x[x > M$R[j]]) > M$G[j] 
             })

	     ## Compare treatment arms for sensitive patinets in the test data
             y.test.inner.sens = y.test.inner[sens.pred]
             patients.test.inner.sens = patients.test.inner[sens.pred,]

             conf = matrix(nrow = 2, ncol = 2, 
             data = c(sum(y.test.inner.sens & !patients.test.inner.sens$treat),  #number of responders in control
                         sum(!y.test.inner.sens & !patients.test.inner.sens$treat), #number of non-responders in control
                         sum(y.test.inner.sens & patients.test.inner.sens$treat),   #number of responders in treatment
                         sum(!y.test.inner.sens & patients.test.inner.sens$treat)), #number of non-responders in treatment
                       byrow = TRUE) 
             p.value[j] = fisher.test(conf, alternative = "two.sided")$p.value
          } # End of eta, R, G combinations loop 
          w = match(min(p.value), p.value)
          return(data.frame(eta = M$eta[w], R = M$R[w], G = M$G[w]))
}

#' @title Get subgroup of sensitive patiens according to the "cvrs" method.
#'
#' @description
#' For one vector of binary responses, get sensitive patients by risk scores and k-means.      
#' In the testing subsets, the risk scores are computed based on the treatment-covariate interaction
#' effects from the training subsets.
#' The risk scores are divided into 2 clusters by a k-means clustering (within each fold).
#'
#' @param  patients - a data frame with patients inormation 
#'         covar - a data frame with covariates
#'         y - a vector of responses 
#'         seed - a seed for random number generator
#'
#' @return A list of 4 :
#'          psens - sensitivity of identifying the sensitive group, one value per simulation run 
#'          pspec - specificity of identifying the sensitive group, one value per simulation run 
#'          sens.pred - predicted sensitivity status (rows = patienst, columns = simulations)
#'          cvrs - a matrix of the risk scores (rows = patients, columns = simulations).  
#' @author Svetlana Cherlin, James Wason

sens.cvrs = function(patients, covar, y, seed) 
{
        if (!is.null(seed)) {
           set.seed(seed)
        }

        ## Divide to folds and keep prevalence of responders/non-responders within folds
        nfolds = 10
        patients.nr = patients[y==0,]
        patients.r = patients[y==1,]

        covar.nr = covar[y==0,]
        covar.r = covar[y==1,]

        y.nr = y[y==0]
        y.r = y[y==1]

        foldid.nr = sample(rep(seq(nfolds), length = length(y.nr))) #non-responders
        foldid.r = sample(rep(seq(nfolds), length = length(y.r))) #responders   
        sens.df = data.frame() 

        ## CV loop
        for (i in seq(nfolds)) {
           which.nr = foldid.nr == i
           which.r = foldid.r == i
           patients.train = rbind (patients.nr[!which.nr,], patients.r[!which.r,])
           patients.test = rbind(patients.nr[which.nr,], patients.r[which.r,])
           covar.train = rbind (covar.nr[!which.nr,], covar.r[!which.r,])
           covar.test = rbind(covar.nr[which.nr,], covar.r[which.r,])                    
           y.train = c(y.nr[!which.nr], y.r[!which.r])
           y.test = c(y.nr[which.nr], y.r[which.r])          

           ## Fit a single-covariate regression model for the training data
           beta.hat = apply (as.matrix(covar.train), 2, function (x)
           {	     
		if (with(patients, exists('sens.true'))) {   
                   mod = glm(y.train ~ patients.train$treat:x, family = "binomial") # for sim data
		} else {
                   mod = glm(y.train ~ patients.train$treat + x + patients.train$treat:x, family = "binomial")  # for real data
		}
                sum.mod = summary(mod)$coefficients
                if (is.na(mod$coeff[names(mod$coeff) == "patients.train$treat:x"])) {0
                } else { mod$coeff[names(mod$coeff) == "patients.train$treat:x"]}

            })

           beta.hat[is.na(beta.hat)] = 0

	   ## Compute  risk scores (CVRS) in the testing data 
           cvrs = apply(as.matrix(covar.test), 1, function (x)  
           {
  	      t(as.matrix(x)) %*% (as.matrix(beta.hat))  
           })           
           sens.fold = data.frame(FID = patients.test$FID, IID = patients.test$IID, cvrs = cvrs) 

           ## Divide testing data to 2 clusters by applying k-means to CVRS within each fold
           sens.fold$sens.pred = FALSE	
           km = kmeans(sens.fold$cvrs, 2)
           if (km$centers[1] > km$centers[2]) {
              sens.fold$sens.pred[km$cluster == 1] = TRUE
           } else {
            sens.fold$sens.pred[km$cluster == 2] = TRUE
           }    
           sens.df = rbind(sens.df, sens.fold)               
        } # End of CV loop
        m = match (paste(as.character(patients$FID), as.character(patients$IID), sep = ":"), 
                   paste(as.character(sens.df$FID), as.character(sens.df$IID), sep = ":"))
        sens.df = sens.df[m,]

	## Compute sensitivity and specificity of the sensitive group selection algorithm for simulated data
        if (with(patients, exists('sens.true'))) {
	   conf = matrix(nrow = 2, ncol = 2, 
              data = c(sum(!sens.df$sens.pred & !patients$sens.true), #predicted non.sens, true non.sens [1,1]
                   sum(!sens.df$sens.pred & patients$sens.true),      #predicted non.sens, true sens [1,2]
                   sum(sens.df$sens.pred & !patients$sens.true),      #predicted sens, true non. sens [2,1]
                   sum(sens.df$sens.pred & patients$sens.true)),      #predicted sens, true sens [2,2]
                   byrow = TRUE)
	   psens = conf[2,2]/(conf[2,2] + conf[1,2])
	   pspec = conf[1,1]/(conf[1,1] + conf[2,1])      
           ret = list(psens = psens, pspec = pspec, sens.pred = sens.df$sens.pred, cvrs = sens.df$cvrs)
        } else {ret  = list(sens.pred = sens.df$sens.pred, cvrs = sens.df$cvrs) }

        return (ret)
}

#' @title Get subgroup of sensitive patiens according to the "cvrs" method for 2 outcomes.
#'
#' @description
#' For two vectors of binary responsess, get sensitive patients by risk scores and k-means.      
#' In the testing subsets, the risk scores are computed based on the treatment-covariate interaction
#' effects from the training subsets.
#' The two sets of the risk scores are divided into 4 clusters by a k-means clustering.
#' cluster 1 corresponsds to low values of both response and response2,
#' cluster 2 corresponss to low values of response and high values of response2,
#' cluster 3 corresponds to high values of response and low values of response2,
#' cluster 4 corresponds to high values of both response and response2.
#'
#' @param  patients - a data frame with patients inormation 
#'         covar - a data frame with covariates
#'         response - a first vector of responses 
#'         response2 - a second vector of binary response
#'         seed - a seed for random number generator
#'
#' @return A list of 4 :
#'          psens - sensitivity of identifying the sensitive group (w.r.t. every cluster), a vector of the length nclust per simulation run 
#'          pspec - specificity of identifying the sensitive group(w.r.t. every cluster), a vector of the length nclust per simulation run 
#'          sens.pred - predicted sensitivity status (rows = patienst, columns = simulations)
#'          cvrs - a matrix of the risk scores (rows = patients, columns = simulations).  
#'          cvrs2 - a matrix of the risk scores for response2 (row = patients, colmns = simulations)
#' @author Svetlana Cherlin, James Wason
#' @importFrom "VGAM" "vglm"
#' @importFrom "VGAM" "binom2.or"


sens.cvrs2 = function(patients, covar, response, response2, seed, nclust)
{
        if (!is.null(seed)) {
           set.seed(seed)
        }
        ## Divide to folds 
        nfolds = 10
        foldid = sample(rep(seq(nfolds), length = nrow(patients)))
        sens.df = data.frame()

        ## CV loop
       for (i in seq(nfolds)) {
           which = foldid == i
           patients.train = patients[!which,]
           patients.test = patients[which,]
           covar.train = covar[!which,]
           covar.test = covar[which,]
           response.train = response[!which]
           response.test = response[which]
           response2.train = response2[!which]
           response2.test = response2[which]

           ## Fit a single-covariate regression model for the training data for response and response2
           reg = apply (as.matrix(covar.train), 2, function (x)
           {
                dat = data.frame(response.train = response.train, response2.train = response2.train, treat.train = patients.train$treat, covarx.train = x)
                if (with(patients, exists('cluster.true'))) { #sim. data  
                   mod = tryCatch({vglm(cbind(response.train, response2.train) ~ treat.train:covarx.train, data = dat, binom2.or)}, error = function(e) e, warning = function(w) w)
                } else { #real data
                   mod = tryCatch({vglm(cbind(response.train, response2.train) ~ treat.train + covarx.train + treat.train:covarx.train, data = dat, binom2.or)}, error = function(e) e, warning = function(w) w)
                }

                if (is(mod, "error") | is (mod, "warning")) {
                      beta.resp = 0
		      beta.resp2 = 0
                } else if (is.na(logLik(mod))){
		      beta.resp = 0
                      beta.resp2 = 0
		} else { 
                   beta.resp = coefficients(mod)[names(coefficients(mod)) == "treat.train:covarx.train:1"]
                   beta.resp2 = coefficients(mod)[names(coefficients(mod)) == "treat.train:covarx.train:2"]
		}
                cbind(beta.resp, beta.resp2)

           })
           reg = t(as.matrix(reg))
           colnames(reg) = c("beta.resp", "beta.resp2")
           reg = as.data.frame(reg)
           ## compute risk scores for response and response2 
           risk.scores = apply(as.matrix(covar.test), 1, function (x)
           {
                 cvrs = t(as.matrix(x)) %*% reg$beta.resp
                 cvrs2 = t(as.matrix(x)) %*% reg$beta.resp2
                 cbind(cvrs, cvrs2)
           })
           risk.scores = t(as.matrix(risk.scores))
           colnames(risk.scores) = c("cvrs", "cvrs2")
           risk.scores = as.data.frame(risk.scores)
           sens.fold = data.frame(FID = patients.test$FID, IID = patients.test$IID, cvrs = risk.scores$cvrs, cvrs2 = risk.scores$cvrs2)

           ## Divide testing data to clusters by applying k-means to CVRS within each fold
	   km = tryCatch ({kmeans(sens.fold[, 3:4], nclust, iter.max = 25, nstart = nrow(sens.fold))}, error = function(e) e, warning = function(w) w)
           if (is(km, "error") | is (km, "warning")) {
              message(km)
              sens.fold$cluster.pred = 0
              return(sens.fold)
           } else {
              centers.sorted = km$centers[order(km$centers[,1], km$centers[,2]),] #sort the clusters according to the centers
              cluster.sorted = numeric(length = length(km$cluster))
	      for (i in 1:nclust) {
                 cluster.sorted[km$cluster == rownames(centers.sorted)[i]] = i
              }
              sens.fold$cluster.pred = cluster.sorted
	   }
	   sens.df = rbind(sens.df, sens.fold)     
        } # End of CV loop

        m = match (paste(as.character(patients$FID), as.character(patients$IID), sep = ":"),
                   paste(as.character(sens.df$FID), as.character(sens.df$IID), sep = ":"))
        sens.df = sens.df[m,]

	## For simulated data, compute sensitivity and specificity with respect to every cluster
        if (with(patients, exists('cluster.true'))) {
	   psens = numeric(length = nclust)
	   pspec = numeric(length = nclust)
	   for (i in 1:nclust) {
	      sens.true = numeric(length = nrow(patients))
	      sens.pred = numeric(length = nrow(patients))
	      sens.true[patients$cluster.true == i] = 1
	      sens.pred[sens.df$cluster.pred == i] = 1 
              conf = matrix(nrow = 2, ncol = 2,
                 data = c(sum(!sens.pred & !sens.true), #predicted sens. group = 0, true sens. group = 0 [1,1]
                      sum(!sens.pred & sens.true),      #predicted sens. group = 0, true sens. group = 1 [1,2]
                      sum(sens.pred & !sens.true),      #predicted sens. group = 1, true sens. group = 0 [2,1]
                      sum(sens.pred & sens.true)),      #predicted sens. group = 1, true sens. group = 1 [2,2]
                   byrow = TRUE)
              psens[i] = conf[2,2]/(conf[2,2] + conf[1,2])
              pspec[i] = conf[1,1]/(conf[1,1] + conf[2,1])
	   }
	   ret = list(psens = psens, pspec = pspec, cvrs = sens.df$cvrs, cvrs2 = sens.df$cvrs2, cluster.pred = sens.df$cluster.pred)
        } else {
           ret  = sens.df
        }
        return (ret)
} 
