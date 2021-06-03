########################################################
## Perform Score Bootstrap for Cox Model Parameter
##
## scoreBSE:   computes score bootstrapped se's
## getScore:   compute score vector from a coxph object 
##
## Written by:  Ayesha Ali
##
## Date started: July 13, 2012
## Last modified: September 2, 2012
##
#########################################################


getScore <- function(cfit, id){
  ############################################
  # Compute the score vector (i.e.
  # observation specific contributions to
  # the first derivative of the partial
  # likelihood at bhat).  The sum of the
  # contributions at bhat should be zero.  
  #
  # This code was hijacked from the example
  # given in ?coxph.detail() for computing 
  # the Schoenfeld residuals.  
  #
  # Function returns score contributions
  # and id of each contribution.
  #
  # Code exploits that p=1 in simulated data.
  # test <- getScore(coxfit, SimData$ID)
  ############################################

  cfitd  <- coxph.detail(cfit)
  events <- cfit$y[,3]==1           # assumes y = (start, stop, event)
  W      <- cfit$weights[events]    # updated sept 2, 2012
  id.fail<- id[events]
  etime  <- cfit$y[events,2]   #the event times --- may have duplicates
  indx   <- match(etime, cfitd$time)
  U.tilde<- W*(cfit$x[events] - cfitd$means[indx])
  
  out <- NULL
  out$id.fail <- id.fail
  out$U.tilde <- U.tilde
  return(out)
} # end of getScore

NewtRaph <- function(beta, U, J, id){
  ###########################################
  ## Needed for Score Bootstrap procedure.
  ## Take one Newtown-Raphson step from beta
  ## using score contributions updated by id
  ############################################
  U.boot <- id*U             # score for beta estimate
  U.star <- sum(U.boot)
  b.new  <- beta + U.star/J  # take a Newton-Raphson step

  return(b.new)  
}  

getVar <- function(U.V, J, wts){
  ##########################################
  ## Needed for Score Bootstrap procedure
  ## use bbotstrapped weights to update
  ## the variance-covariance amtrix for beta
  ## exploits p=1 element in beta
  ###########################################

  U.bootV<- wts*U.V          # score for wald statistic
  Var.U <- sum(U.bootV^2)    # variance of the (uni-variate) score
  covB  <- (Var.U/(J^2))     # var-covar matrix
  return(covB)
}  
  
 
processSBparams <- function(Bw, B, T, a){
  #######################################################
  ## Process outputs from Score Bootstrap
  ## i.e. save empirical SES, compute bootstrapped CI's
  ## and compute total run time.
  #######################################################
  
  ## compute bootstrapped CI's 
  #lrank <- ceiling(.025*B) # old
  #urank <- floor(.975*B)   # old
  lrank <- round(.025*(B+1))
  urank <- (B+1)-lrank
  Bbeta <- sort(Bw[,1])
  Bwald <- sort(Bw[,2])

  ## save standard error, 95% CI, and wald statistic 
  out      <- NULL
#  out$bse  <- sqrt(var(Bbeta-Bw[B,1])) 
  out$bse<- sqrt(var(Bbeta)) 
#  out$mcbse<- sqrt(var(Bbeta)) 
  out$wald <- sum(Bwald > T)/B       # rejection rate
  out$mwald<- mean(Bwald)            # added sept2, 2012
  out$waldsum <- summary(Bwald)
  out$waldci  <- (c(Bwald[lrank], Bwald[urank]))
  out$bci  <- (c(Bbeta[lrank], Bbeta[urank]))
  out$cpu  <- Sys.time() - a

  return(out)
}
 
scoreBSEfunc <- function(coxfit, id, B=-1, Radem="T"){
#################################################### 
# Score bootstrap for censored data
# Input is: coxfit (coxph object), observation id's,
# B= number bootstrap samples, and 
# Radem=T means Rademacher weights, else
# Mammen (1993) weights for skewed data used.
# Uses wild bootstrap to resample score
# contributions, and then recopmutes one
# Newton-Raphson update to get bootstrapped
# betas, with the Hessian from the original 
# full data.
# See Kline and Santos (2011) for details
#  
# This uses lappy() to do the bootstrap loop
# Code exploits that p=1 for Young's simulated data. 
#################################################### 

  if (B < 1) {
    stop("Number of bootstraps needs to be at least 1 in scoreBSEfunc()");
  }

  start.time  <- Sys.time()          # needed for tracking cpu time
    
  ## compute observation contributions to the score at coeff0
  ## this is used for each bootstrapped sample
  update  <- getScore(coxfit, id)
  id.fail <- update$id.fail          # id of failed subjects
  Nfails  <- length(id.fail)         # number subjects who failed

  U.tilde <- update$U.tilde          # score contributions by failure times
  U.tildeV<- residuals(coxfit, weighted = T, type = "score", collapse=id) # by case
  J.tilde <- coxph.detail(coxfit)$imat
  J.star  <- sum(J.tilde)  

  ## DEFINE WEIGHTING SCHEMES: RADEMACHER OR WILD?
  if(Radem=="T"){ 
     wwts <- c(-1,1)     # wild weights
	 pwts <- c(0.5, 0.5) # prob. of wild weights
  } else {
    ## THESE NEED TO BE TESTED!! BUT NOT NEEDED FOR NOW.
    wwts <- c((1-sqrt(5))/2, (1+sqrt(5))/2)
	pwts <- c((5+sqrt(5))/2, 1-(5+sqrt(5))/2) 
  }	 

  b.hat  <- coxfit$coeff[1]
  
  ## CREATE (B-1) BOOTSTRAP SAMPLES OF: Nsub id's sampled with replacement
  ## randomly generate Rademacher weights
  ## and apply to observation contribution to score
  id.uni <- unique(id)
  Nsub   <- length(id.uni) 

  ## index.f gives the index to the weights in bwts for the failed subjects 
  ## it first takes creates a vector of Nfails 0's, concatenated with
  ## Nfails 1's, ..., Nfails (B-2)'s, then multiplies by Nsub to compute
  ## constant to add to id.fail and finally compute the index for bwts.
#  browser()
  
  # bwts   <- matrix(sample(wwts, Nsub*(B-1), replace=T, prob=pwts), nrow=B-1, byrow=T)
  bwts   <- matrix(sample(wwts, Nsub*B, replace=T, prob=pwts), nrow=B, byrow=T) # maali 2012-10-21

  id.bwt <- bwts[,id.fail] 

 
  do.one <- function(b)
  {
    ##########################################################
    ## Build and analyze one bootstrap sample. 
    ## INPUT: b = handle to row of id.bwt to be processed. 
    ## OUTPUT: vector containing 1 beta estimate and 
    ##                           1 wald statistic value
    ##
    ## 'U.tilde', 'U.tildeV', 'J.star', 'b.hat' and 'id.bwt' are the 
	## same across workers when parallelized. We define function 
	## internally so that these paramters need not be passed in 
	## as arguments to do.one().
    ##########################################################

	  ## update beta estimate	  
	  ## update the variance estimator
	  ## reutrn new beta estimate and new wald statistic
	  b.new     <- NewtRaph(b.hat, U.tilde, J.star, id.bwt[b,])
      b.var     <- getVar(U.tildeV, J.star, bwts[b,])
      bootparam <- c(b.new, ((b.new^2)/b.var))

      return(bootparam)
  }
  
  # 1- define result holder  

#  my.ncpus <- detectCores()                   # can run 2 threads per cpu detected
#  M  <- min(my.ncpus,(B-1))
#  cl <- makeCluster(M)
#  clusterExport(cl, c("NewtRaph", "getVar"))  # fn's needed in do.one()
  Bwald <- NULL

  # 2- serial list apply

#  Bwald <- do.call(c, parLapply(cl, 1:(B-1), do.one) )
#  Bwald <- lapply(1:(B-1), do.one)
  Bwald <- lapply(1:B, do.one)
  # Bwald = 4xB vector containing (b1, b2, w1, w2)xB 
  # cat("length of Bwald: ", length(Bwald), "\n")

# 3- destroy cluster 
#  stopCluster(cl)

  # 3- convert Bwald into an output matrix with each column
  # holding results from respective bootstrap samples

#  Bwald  <- matrix(c(unlist(Bwald), coxfit$coeff[1], coxfit$wald.test), nrow=B, byrow=T)  
  Bwald  <- matrix(unlist(Bwald), nrow=B, byrow=T)  
  #print(tail(Bwald))
 
  output <- processSBparams(Bwald, B, coxfit$wald.test, start.time)
  return(output)
}
