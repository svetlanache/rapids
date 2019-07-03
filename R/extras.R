#' @title Plot the risk scores.
#' @description The function plots the risk scores.
#' @param cvrs A data frame of the  risk scores.
#' @param sens.pred A data frame of the predicted sensitivity statuses.
#'
#' @details A plot of the risk scores is produced and saved into a "CVRSplot.png" file.
#' @author Svetlana Cherlin, James Wason
#' @examples
#' #Analyse simulated data and plot the risk scores
#' #"cvrs" method
#' data(simdata)
#' sig = 0.05
#' group.prop.sig = 0.2
#' seed = 123
#' plotrs = T
#' method = "cvrs"
#' eta = NULL
#' R = NULL
#' G = NULL
#' simres.cvrs = analyse.simdata(simdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)
#'  #Plot the risk scores
#' cvrs.plot(simres.cvrs$cvrs, simres.cvrs$sens.pred)
#' @seealso
#' \code{\link{analyse.simdata}} and \code{\link{analyse.realdata}} functions; \code{\link{plot}} and \code{\link{print}} methods.
#' @export
#' @importFrom "ggplot2" "ggplot" "aes" "geom_point" "labs" "theme" "scale_colour_continuous" "guides" "guide_legend" "ggsave"

cvrs.plot = function(cvrs, sens.pred) 
{


   cvrs = as.data.frame(cvrs)
   sens.pred = as.data.frame(sens.pred)
   N = nrow(cvrs)
   cvrs.df = data.frame()
   for (i in 1:ncol(cvrs)) {
      cvrs.df = rbind(cvrs.df, data.frame(ids = seq(1:N), cvrs = cvrs[,i], sens.pred = as.numeric(sens.pred[,i]), run = i))
   }
   ggplot(cvrs.df, aes(x = ids, y=cvrs, color = sens.pred)) +
   geom_point(fill = "transparent", shape = 21, size = 2, stroke = 0.5) +
   labs(x = "Patients", y = "CVRS") +
   theme(legend.position = "right") +
   scale_colour_continuous(name  ="Patients",
                            breaks=c(0, 1),
                            labels=c("Non-sensitive", "Sensitive")) +
   guides(colour = guide_legend(override.aes = list(shape = 15)))
   
   ggsave("CVRSplot.png")
}


#' @title Print an object of class "rapids".
#' @param x An object of class "rapids" 
#' @details For simulated data, the following operating characteristics are printed: (i) the power for the overall test; (ii) the power for the sensitive group test; (iii) the power for the adaptive design; (iv) mean sensitivity and specificity of identifying the sensitive group; (v) mean response rate in the sensitive group on the treatment arm. For real data, the following operating characteristics are printed: (i) p-value for the overall test; (ii) p-value for the  sensitivity group test; (iii) estimated response rate in the sensitive group on the treatment arm.
#' @author Svetlana Cherlin, James Wason
#' @examples
#' #Analyse simulated data and print the results
#' #"cvrs" method
#' data(simdata)
#' sig = 0.05
#' group.prop.sig = 0.2
#' seed = 123
#' plotrs = T
#' method = "cvrs"
#' eta = NULL
#' R = NULL
#' G = NULL
#' simres.cvrs = analyse.simdata(simdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)
#' #Print the results
#' print(simres.cvrs)
#' @seealso
#' \code{\link{cvrs.plot}} function and \code{\link{print}} method.
#' @export

print.rapids=function (x,...)
{

   if (with(x, exists('psens'))) { #Simulated data
      if (length(x$estimate.rr) == 1) {	   
         sens.print = "Sensitivity of identifying the sensitive group: "
         spec.print = "Specificity of identifying the sensitive group: "
         rr.print = "Response rate in the sensitive group on the treatment arm: "
      } else {
         sens.print = "Mean sensitivity of identifying the sensitive group: "
         spec.print = "Mean specificity of identifying the sensitive group: "
	 rr.print = "Mean response rate in the sensitive group on the treatment arm: "
      }
      toprint = paste ("\nPower for the overall test: ", round(x$pwr.overall, 3),
          "\nPower for the sensitive group: ", round(x$pwr.group, 3),
          "\nPower for adaptive design: ", round(x$pwr.adaptive,3),
          "\n", sens.print, round(mean(x$psens), 3),
          "\n", spec.print, round(mean(x$pspec), 3), 
	  "\n", rr.print, round(mean(x$estimate.rr, na.rm = T), 3), "\n", sep = "")
   } else { #Real data
      toprint = paste("\nP-value for the overall test: ", x$pval.overall,
          "\nP-value for the sensitive group test: ", x$pval.group,
	  "\nResponse rate in the sensitive group on the treatment arm: ", round (x$estimate.rr, 3), "\n", sep = "")
   }
   cat(toprint)
}

#' @title Plot risk scores, sensitivity and specificity from a "rapids" object.
#' @param x An object of class "rapids" 
#' @param   ... Other parameters to plot 
#' @details For simulated data, the following is plotted: ROC plots for predicting the sensitivity status, boxplots for the AUC, sensitivity and specificity, risk scores for the "cvrs" method (the risk scores from all of the simulation runs are shown on the same plot).
#' @details For real data, risk scores are plotted for the "cvrs" method. For real data  analysed with the "cvasd" method, a message "nothing to plot" is printed.
#' @author Svetlana Cherlin, James Wason
#' @examples
#' #Analyse simulated data and plot the risk scores
#' #"cvrs" method
#' data(simdata)
#' sig = 0.05
#' group.prop.sig = 0.2
#' seed = 123
#' plotrs = T
#' method = "cvrs"
#' eta = NULL
#' R = NULL
#' G = NULL
#' simres.cvrs = analyse.simdata(simdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)
#' #Plot the results
#' plot(simres.cvrs)
#' @seealso
#' \code{\link{cvrs.plot}} function and \code{\link{print}} method.
#' @export
#' @importFrom "graphics" "plot" "boxplot" "legend"
#' @importFrom "pROC" "roc"

plot.rapids = function(x, ...)  
{
   x$sens.pred = as.matrix(x$sens.pred)
   x$sens.pred = apply(x$sens.pred, 2, as.numeric)
   x$sens.pred = as.matrix(x$sens.pred)
   N = nrow(x$sens.pred)
   runs = ncol(x$sens.pred)
 

   if (with(x, exists('cvrs')) & with(x, exists('psens'))) { #CVRS, sim.data 
       x$cvrs = as.matrix(x$cvrs)
       if (runs > 1) {
          par(mfrow = c(2, 2))
       } else {
	  par(mfrow = c(2, 1))
       }
   } else if (with(x, exists('psens'))) { #CVASD, sim.data
       if (length(x$psens) == 1) {
          par(mfrow = c(1, 1))
       } else { 
          par(mfrow = c(1, 2))
       }
   } else if (with(x, exists('cvrs'))) { #CVRS, real data
       x$cvrs = as.data.frame(x$cvrs)
       par(mfrow = c(1, 1), mar = par('mar') + c(0,0,0,6)) 
   } else { #CVASD, real data
      print ("Nothing to plot")
   }

   ## Plot ROC curve for identifying the sensitve group for simulate dated
   if (with(x, exists('psens'))) {
      res = apply(x$sens.pred, 2, function(y)
      { 
         if (length(table(y)) == 2) {  
            roc = roc (y, x$patients$sens.true, levels=base::levels(as.factor(x$patients$sens.true)), direction="<")
         }
      })
      auc = res[[1]]$auc
      par(pty = "s")
      plot(res[[1]], main = "Sensitive group", ...)
      if (runs == 1) {
	 legend("bottomright", paste("AUC: ", signif(auc, 2), sep = ""), box.lwd = 0)
      } else {
         for (i in 2:runs) {
            if (with(res[[i]], exists('auc'))) {
               plot(res[[i]], add = TRUE)
               auc = c(auc, res[[i]]$auc)
            }
         }
         legend("bottomright", paste("Mean AUC: ", signif(mean(auc), 2), sep = ""), box.lwd = 0)
         boxplot(x$psens, x$pspec, auc, names = c("Sensitivity", "Specificity", "AUC"), ylim = c(0,1), boxwex = 0.8, las = 2, ylim = c(0,1), main = "Sensitive group", ...)
      }
   }

   ## Plot CVRS
   if (with(x, exists('cvrs'))) {   
      cvrs.df = data.frame()
      for (i in 1:ncol(x$cvrs)) {
         cvrs.df = rbind(cvrs.df, data.frame(ids = seq(1:N), cvrs = x$cvrs[,i], sens.pred = as.numeric(x$sens.pred[,i]), run = i))
      }
      plot(cvrs.df$ids, cvrs.df$cvrs, col = x$sens.pred + 1, xlab = "Patients", ylab = "CVRS", main = "CVRS", ...)
      
      legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, c("Non-sensitive", "Sensitive"), pch=c(1), title="Patients", col = c(1,2))   
   }
}

#' @title  Perform permutation test for the real data.
#' @description
#' For each permutation run, a test statistic for the interaction effect between the treatment and the predicted sensitivity status is computed.
#' A permutation p-value is computed as the proportion of simulations where the test statistic is larger than the test statistic for the original data.
#' @param  realdata Input data (see \code{realdata} object).
#' @param  res An object of class "rapids" containing the results for the design for the original data.
#' @param  method An input method ("cvasd" or "cvrs").
#' @param  reps Number of permutations.
#' @details Perform permutation test for the real data.
#' @return A list of 5:
#' @return stat0.group:  Statistic for the interaction effect between the treatment and the sensitivity status, for the original data.
#' @return treat: A matrix of treatment allocations, one column per permutation run.
#' @return sens.pred: A matrix of sesnsitivy statuses, one column per permutation run.
#' @return stat.group:  A matrix of permutation statistics, one column per permutation run.
#' @return ppval.group: P-value from the permutation test.
#' @author Svetlana Cherlin, James Wason
#' @examples
#' #Analyse data with "cvrs" method and perform permutation test
#' data(realdata)
#' sig = 0.05
#' group.prop.sig = 0.2
#' method = "cvrs"
#' seed = 123
#' plotrs = T
#' eta = NULL
#' R = NULL
#' G = NULL 
#' realres.cvrs = analyse.realdata(realdata, sig, group.prop.sig, method, eta, R, G, seed, plotrs)
#' #Permutation test
#' reps = 10
#' res.perm = permutation.test(realdata, realres.cvrs, method, reps)
#' @seealso
#' \code{\link{analyse.realdata}} function.
#' @export 

permutation.test <- function(realdata, res, method, reps)
{

  treat0 = realdata$patients$treat #original treatment allocation
  stat0.group = res$stat.group #statistic for the treatment effect in the subgroup

  st <- replicate(reps, perm(realdata, res, method, treat0))

  treat = as.matrix(as.data.frame(st[rownames(st) == "treat",]))
  colnames(treat) = seq(1:reps)
  rownames(treat) = seq(1:nrow(realdata$patients))
  sens.pred = as.matrix(as.data.frame(st[rownames(st) == "sens.pred",]))
  colnames(sens.pred) = seq(1:reps)
  rownames(sens.pred) = seq(1:nrow(realdata$patients))
  stat.group = st[rownames(st) == "stat.group",]

  ppval.group = (1+sum(stat.group > stat0.group))/(1+reps) #permutation p-value for the treatment-sensitivity interaction effect

  res.perm = list(stat0.group = as.numeric(stat0.group), treat = treat, sens.pred = sens.pred,
             stat.group = as.numeric(stat.group),
             ppval.group = as.numeric(ppval.group))

  return (res.perm)
}  

  perm = function(realdata, res, method, treat0) {

     ## permute treatment labels
     realdata$patients$treat = sample(treat0, replace=FALSE) 

     ## get sensitive group
     sens = switch (method, cvasd = sens.cvasd(realdata$patients, realdata$covar, realdata$response, res$eta, res$R, res$G, NULL),
                             cvrs = sens.cvrs(realdata$patients, realdata$covar, realdata$response, NULL))

     ## fit glm
     mod = glm(realdata$response ~ sens$sens.pred + realdata$patients$treat + sens$sens.pred:realdata$patients$treat, family = "binomial")

     ## statistic for the treatment effect in the sensitive group
     theta = mod$coeff

     if (is.na(theta[names(theta) == "sens$sens.predTRUE:realdata$patients$treat"])) { #no sensitive group identified
        stat.group = -999
     } else {
        sum.theta = summary(mod)$coeff 
        stat.group = sum.theta[names(theta) == "sens$sens.predTRUE:realdata$patients$treat", colnames(sum.theta) == "z value"]
     }

     return(list(stat.group = stat.group, treat = realdata$patients$treat, sens.pred = sens$sens.pred))
  }
