### ###############################################
###
###  Fitting Time-dependent weighted Cox Models
###  Started by: Zhe Wei
###  Modified by: Ayesha Ali
###
###  This file contains the functions to fit a 
###  TD-weighted Cox model with the proper standard
###  error computations for marginal structural 
###  Cox model estimation.
### 
###  Last Updated:  April 23, 2012
### 
### ###############################################
### ###############################################


source("coxtdw.utils.R")

########################################################
### cox function starts here
########################################################
### maali
### - added robust parameter so we can use interface from the
###   survival.coxph package

## XXX max.iter, eps should be in a param called 'control' as
##    in coxph()
coxtdw <-
  function(formula = formula(data), data = parent.frame(), weights,
           max.iter = 10, eps=1E-4, method="breslow",
           robust=TRUE, init = NULL)
{

  
### 1- GET DATA FROM PARAMS ###
  
  ## Copied from coxph directly ##
  call   <- match.call()
  m      <- match.call(expand.dots = FALSE)
  temp   <- c("", "formula", "data", "weights", "subset")
  m      <- m[match(temp, names(m), nomatch = 0)]
  special<- c("strata", "cluster")
  
  Terms  <- if (missing(data))           # get terms
    terms(formula, special)
  else terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]]    <- as.name("model.frame")
  m         <- eval(m, parent.frame())
  if (NROW(m) == 0) 
    stop("No (non-missing) observations")
  
  Y         <- model.extract(m, "response")  # get response
  if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")
  
  ## browser()
  
  attr(Terms, "intercept") <- 1
  strats  <- attr(Terms, "specials")$strata
  cluster <- attr(Terms, "specials")$cluster #get id
  dropx   <- NULL
  if (length(cluster)) {
    if (missing(robust)) {
      robust <- TRUE
    }
    tempc <- untangle.specials(Terms, "cluster", 1:10)
    ord   <- attr(Terms, "order")[tempc$terms]
    if (any(ord > 1)) {
      stop("Cluster can not be used in an interaction")
    }
    cluster <- strata(m[, tempc$vars], shortlabel = TRUE)
    dropx   <- tempc$terms
  }
  
  ## maali
  id <- cluster
  if(is.null(id)) {
    id <- 1:nrow(Y)
  }
  ## XXX check if there was a "cluster" clause in
  ##     input formula
  
  subjects <- unique(id) # match id.uni in coxtdw.utils::ui.w
  nsubject <- length(subjects)
  
   
  if (length(strats)) {
    temp  <- untangle.specials(Terms, "strata", 1)
    dropx <- c(dropx, temp$terms)
    if (length(temp$vars) == 1) {
      strata.keep  <- m[[temp$vars]]
    } else {
      strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
    }
    strats <- as.numeric(strata.keep)
  }
  
  if (length(dropx)) {
    newTerms  <- Terms[-dropx]
  } else {
    newTerms <- Terms
  }
  
  X      <- model.matrix(newTerms, m)
  assign <- lapply(attrassign(X, newTerms)[-1], function(x) x - 1)
  X      <- X[, -1, drop = FALSE]
  type   <- attr(Y, "type")
  if (type != "right" && type != "counting") {
    stop(paste("Cox model doesn't support \"",
               type, "\" survival data", 
               sep = ""))
  }
  
  weights <- model.extract(m, "weights")   # get weights
  
  
  ##cat("blah")                         
  ##browser()
  
  
  
  ## ################################################
  ## compute preliminary calcs for S.tilde values using Breslow
  ## approximation for the Anderson-Gill form
  ## ################################################
  ## browser()
  start  <- Y[, 1]
  stopp  <- Y[, 2]
  event  <- Y[, 3]


  ## ###
  
  ## XXX need to check that events don't re-occur, you only die once
  ## XXX need to check that intervals don't overlap for a given subject
  ## XXX print if we have case weights or time-dependent ones
  ## XXX need to check if any events occur at all
  ## XXX check for coxph()-specific parameters such as: ‘cluster’,
  ##     ‘strata’, ‘Surv’, ‘survfit’, ‘pspline’, ‘frailty’,
  ##     ‘ridge’, 'tt' and report error to user
  ## XXX compare likelihood values to thse in coxph()
  ##
  ## XXX check length of input vectors
  ## XXX check weights are >0
  ## XXX 
  

### 2 - NOW FILL IN DEFAULT VALUES ###
  
  ## maali - use uniform weights if none are provided
  if ( is.null(weights) ) {
    weights <- rep(1,length(stopp)); 
  }

  
### 3 - CALCULATE HANDY VALUES, PREPROCESS INPUTS ###
  
  ## ###################################
  ## need to identify the failure time
  ## of each subject, and the covariate
  ## values of all subjects at risk
  ## at those times
  ## ###################################

  ## XXX preconditions
  ## - start & stopp times are given +ve values
  ## - start < stopp for all observations
  ## - event (i.e. death) only occurs once per id (i.e. subject)
  ## - after event occurs, subject leaves study
  
  n   <- nrow(X) # number of observations
  p   <- ncol(X) # number of parameters (number of covariates or "factors")

  ## ################################
  ## subtract mean from each covar
  ## as this makes regression more
  ## stable.
  ## #################################
  X.fit <- X
  temp <- apply(X, 2, mean)
  for(i in 1:n) {X.fit[i,] <- X[i,] - temp}
  
  ## notes to me
  d.t  <- array(0,dim=c(1,n))  # failure times for those who failed
  d.wt <- array(0,dim=c(1,n))  # weight at failure times
  d.Xt <- array(0,dim=c(n,p))  # covariates at failure times

  ## browser()                           
  for (i in 1:n) {
    if (event[i]==1) {
      d.t[i]   <- stopp[i]
      d.Xt[i,] <- X.fit[i,]
      d.wt[i]  <- weights[i]
    }
  }

  ## browser()
  d.nidx <- (d.t!=0) # index in n obs of events (deaths)
  d.Xt   <- d.Xt[d.nidx,]
  d.wt   <- d.wt[d.nidx]
  d.t    <- d.t[d.nidx]
  d.id   <- id[d.nidx]  # index id of subjects who failed

  N <- length(d.wt);                 # number of failures (deaths)

  ## browser() # d.Xt should be (<=n) x p

    fid  <- array(dim=c(nsubject,1)) # final stop time by id
    sid  <- array(dim=c(nsubject,1)) # start time by id
    efid <- array(dim=c(nsubject,1)) # event at final stop by id
    Xfid <- array(dim=c(nsubject,p)) # covariates at final time by id
    wfid <- array(dim=c(nsubject,1)) # weight at final time by id

    i <- 1
    for (s in subjects) {
      fid[i]  <- max(stopp[id==s])
      sid[i]  <- min(start[id==s])
      
      ## vector of logicals, all false except for i's final observation
      ## fidx = final time index into parallel arrays
      fidx <- ( stopp == fid[i] ) & ( id == s )
      
      
      efid[i]  <- event[ fidx ]
      Xfid[i,] <- X.fit[ fidx, ]
      wfid[i,] <- weights[ fidx ]

      i <- i+1;
    }

    ##
    ## Ykti is Y_k(t_i) in the paper, which gives the survival status
    ## of suject 'k' at the final stop time for (last observation of)
    ## subject 'i'.
    ## 
    ## Ykti is initially size dim=c(nsubject,nsubject), but is stored
    ## more briefly as a nsubject-by-N matrix afterward, since in the
    ## expression for likelihood, the delta_i term cancels out the
    ## 'i' terms where the subjects survive the last observation
    ## 
    ## Xkti is the same but for X_k(t_i) each entry of which is
    ## a column vector of length p (size $p \times 1$)
    ##
    
    ##browser()
    Ykti <- array(dim=c(nsubject,nsubject))
    Xkti <- array(dim=c(nsubject,nsubject,p))
    wkti <- array(dim=c(nsubject,nsubject))

    i <- 1
    for (sub.i in subjects) {
      fi <- fid[i]                     # final stop time of i
      k <- 1
      for (sub.k in subjects) {
        fk <- fid[k]                   # final stop time of k
        sk <- sid[k]                   # first start time of k
        
        ## if
        ##   subject k's stop time is greater than subject i
        ## then
        ##   at the stop time of subject i, subject k
        ##   is still alive, so survival status $Y_k(t_i) = 1$
        ## else
        ##   subject k has left the study, so we take the last
        ##   observation of subject k, using efid -- event at
        ##   final time by id -- and set survival status to
        ##   Y = 1-efid[k]
        ##   NO IT'S JUST 0!!!
        ## endif         
        
        ## if
        ##   subject k's stop time is greater than subject i
        ## then
        ##   subject k is still alive & in the study, so 
        ##   at the stop time of subject i, the value of
        ##   subject k's covariates X is given by taking
        ##   k's maximum stop time that is less than fi.
        ##   the expression here uses the fact that using R's
        ##   which.max preserves index numbers as the name,
        ##   then uses that to index into X which is a parallel
        ##   array to 'stopp'. Same for weight.
        ##   
        ## else
        ##   subject k has left the study, so we take the last
        ##   observation of subject k
        ## endif
        
        ## A:      fi
        ##         |
        ## -----]   
        
        ## B:      fi
        ##         |
        ## --------]   
        
        ## C:      fi
        ##         |
        ## -----]    (------------
        
        ## D:      fi
        ##         |
        ## --------]    (---------
        
        ## E:      fi
        ##         |
        ## (----------]   
                
        ## F:      fi
        ##         |
        ##         (--------------
        
        ## G:      fi
        ##         |
        ##              (---------
        
        
        ## A, C takes the last interval's value at end of interval
        ## B, D, & E takes the current interval's value at end of interval
        ## F & G are the same as each other -- arbitrary
        
        
        if (fk >= fi) {
          kb <- start[id==sub.k]  # k's set of stop times (begin)
          ke <- stopp[id==sub.k]  # k's set of stop times (end)
          
          ## test if there's an intersecting interval 
          bf <- (kb < fi) & (ke >= fi)  # "best fit" of index (cases BDE)
          
          ## try last interval ending before fi
          if (!any(bf)) {
            bf <- ke < fi # new "best fit" of index (cases AC)
            
            ## if all intervals start AFTER fi, arbitrarily
            ## take the first observation of the covariates
            if (!any(bf)) {
              bf[1] <- T # worst "best fit" of index (cases FG)
            }
          }
          
          Ykti[k,i] <- ifelse(sk<fi, 1, 0)
          bfidx <- as.numeric(names(which.max(ke[bf])))
          Xkti[k,i,] <- X.fit[bfidx,]
          wkti[k,i] <- weights[bfidx]
        } else {
          Ykti[k,i] <- 0 # 1-efid[k]
          Xkti[k,i,] <- Xfid[k,]
          wkti[k,i] <- wfid[k]
        }
        
        k <- k+1
      } # end for k
      i <- i+1
    } # end for i 
        
    ## survival indicator Y and covariate value X for subject
    ## k at the i-th failure time ('i' is counting only failed
    ## subjects for d.* variables)
    
    ## notes to me
    d.Ykti <- array(dim = c(nsubject, N)) 
    d.Xkti <- array(dim = c(nsubject, N, p)) # cov vals
    d.wkti <- array(dim = c(nsubject, N)) 
    
    d.Ykti <- Ykti[,efid==1]
    d.Xkti <- Xkti[,efid==1,]
    d.wkti <- wkti[,efid==1]
    
  if(p==1){
    d.Xkti.N <- d.Xkti[efid==1,]  # X values for those who failed
  } else{
    d.Xkti.N <- d.Xkti[efid==1,,]  # X values for those who failed
  }

### 4 - INITIALIZE ITERATIVE PARAMS ###
  
  if(is.null(init)) {
    init <- matrix(0, nrow=p,ncol=1)
  }
  
  if(length(init)!=p) {
    stop(cat("init value has wrong length() ",
             length(init), ", should be of length ",
             p, ".\n"));
  }
  beta.old <- matrix(init,nrow=p,ncol=1)
  iter.num <- 1  # iteration index

  ## The following terms are used in the main iteration to 
  ## compute the log-likelihood value for convergence, but
  ## are also needed later for variance calculations
  W    <- array(dim = c(N,1)) 
  C    <- array(dim = c(nsubject,N))  
  Xbar <- array(dim = c(N,p))          # S1.tilde/S0.tilde in MSCM paper
  C.wt <- array(dim = c(nsubject, N))  # C/colsSums(C), used in I.tilde and for variance

  I.tilde.beta <- array(dim=c(p,p))
  
  W <- d.wt


### 5 - BEGIN ITERATIONS ###

  while(iter.num <= max.iter) {
    
    b <- beta.old
    
    
    ## notes to me
    B <- array(dim = c(N,1)) 
    B <- d.Xt %*% b      # N x 1 vector
    
    expval <- array(dim = c(nsubject, N)) # note to me
    i <- 1
    for (s in subjects) {
      for (j in 1:N) {
        if (p==1) {
          expval[i,j] <- d.Xkti[i,j] * b
        } else {
          expval[i,j] <- array(d.Xkti[i,j,],dim=c(1,p)) %*% b
        }
      }
      i <- i+1
    }
    
    C    <- d.Ykti * d.wkti * exp( expval )
    if(N==1){
	  sumC <-  sum(C)         # S0.tilde in MSCM paper
	} else sumC <-  colSums(C)

    A <- B - log(colSums(C)); # sum over k
    LL <- sum(W * A);         # sum over i
    
    cat("LL: ", format(LL, digits=6), "\n")
    
    ## CONVERGENCE CRITERIA
    ##
    
    ## LL log likelihood
    ## d.wt[k]
    ## l_{beta} =                                     # LL
    ##   sum_{i=1}^n {
    ##     w_i(t_i)*\delta_i*                         # W 
    ##     [
    ##       \beta^t * X_i(t_i) -                     # B
    ##       log sum_{k=1}^n (                      
    ##         Y_k(t_i) * w_k(t_i) *                  # ---- C    
    ##           exp( \beta^t * X_k(t_i) )            # _|
    ##       )
    ##     ]
    ##   }
    ##
    ##  LL = sum W*(A)           
    ##     = sum W*(B - log sum C)                    # expand A
    ##
   
    ## browser()


    ## loglik, newlik - scalars
    ## iter, - int
    ## beta, newbeta 1xp
    
    
    ## u is score px1, I is pxp
    ## I x = u
    ## x = I^-1 u
    
	## AYESHA'S CODE ADDED APRIL23, 2012 
	## BEGIN
    ## Next we compute Xbar and C.wt, both of which will be used
	## for variance estimation at beta-hat.  We also compute the
	## score and information matrix for the Newton-Raphson update.
    ## Need to accommodate tied deaths...
	
	for(j in 1:N){
		C.wt[,j] <- C[,j]/sumC[j]
        if (p==1){
		  Xitj    <- array(d.Xkti[,j],dim=c(nsubject,1))
 	      Xbar[j] <- sum(Xitj*C.wt[,j])
        } else {
		  for(k in 1:p){
		    Xitj      <- array(d.Xkti[,j,k],dim=c(nsubject,1))
            Xbar[j,k] <- sum(Xitj*C.wt[,j])
		  }
        }
	}
    ## END OF AYESHA'S NEW CODE 
	
	## COMPUTATION OF U.tilde and I.tilde below
	## were updated using the notation/quantities
	## Ayesha computed above:  C.wt, Xbar, Xitj
	
	## update score vector U.tilde
    U.tilde <- array(0, dim=c(p, 1))

	if(p==1){
	   U.tilde <- sum(W*(diag(d.Xkti.N) - Xbar))
	   } else {
	   for(k in 1:p){
		    U.tilde[k] <- sum(W*(diag(d.Xkti.N[,,k]) - Xbar[,k]))
	   }
	}
 
#  	cat("C: \n")
#	print(C)
#	cat("sumC: \n")
#	print(colSums(C))
#	cat("C.wt: \n")
#	print(C.wt)

 ## update information matrix I.tilde
    I.tilde <- array(0, dim=c(p, p))
    for(i in 1:N){
	  if(p==1){
	    X.ti     <- array(d.Xkti[,i], dim=c(nsubject,p))
	    Xbar.ti   <- array(Xbar[i], dim=c(p,1))
	    for(j in 1:nsubject){
	      I.tilde <- I.tilde + W[i]*((X.ti[j])^2)*C.wt[j,i]
        }
        I.tilde <- I.tilde - W[i]*((Xbar.ti)^2)
	  } else {
	    X.ti     <- t(array(d.Xkti[,i,], dim=c(nsubject,p)))
	    Xbar.ti   <- array(Xbar[i,], dim=c(p,1))
		for(j in 1:nsubject){
		  I.tilde <- I.tilde + W[i]*(aaT(X.ti[,j])*C.wt[j,i])
		}
        I.tilde <- I.tilde - W[i]*(aaT(Xbar.ti))
	  }
	}
    ## END of using Ayesha's code to update U.tilde, I.tilde
    beta.new   <- beta.old + solve(I.tilde) %*% U.tilde
    
    ##  cat(beta.new, "\n")

    ## 1. check for final convergence condition
    
    ## 2. check for local divergence 
    ##    2a. was this a step in the wrong direction? then half step.
    ## 3. otherwise, we are converging, so don't half-step
    ## 4. check for convergence
    
    
    ## old verison
    ## cc - convergence criteria

    cc <- max(abs(beta.new - beta.old)/(abs(beta.old + 1e-16) ))
    if (cc < eps) {
	  I.tilde.beta <- I.tilde 
      break
    } else{
      if(iter.num >= max.iter) {
        cat("estimates did not converge in ",
            max.iter, "iterations...stopping\n")
        stop
      }
    }
    
    beta.old <- beta.new              # update beta estimate
    iter.num <- (iter.num + 1)        # update iteration number
    
  } # end while (iter<maxiter)
  
  beta.hat <- beta.new
  
  cat("Estimates converged...computing standard errors...\n\n")
  
  
### 5 - BEGIN VARIANCE ESTIMATION ###
##  THIS SECTION WAS RE-CODED BY AYESHA APRIL 23, 2012
##  USING NOTATION DEVELOPED IN THE VALIDATE COXTDW FILE
##  (i.e. no S0, S1 or S2, but W, C, Xbar, etc.
## so utility file ui.w is no longer used at all.
  
  ## ###########################################
  ## need to compute covariates for ALL
  ## subjects at failure times (eg. X_i(t_j))
  ## ###########################################

 
  C.wt.N        <- C.wt[efid==1,]
  
  U.tilde.beta <- array(0, dim=c(nsubject,p))
  for(i in 1:nsubject){

    term1 <- rep(0,p)
	if(efid[i]==1){
	  # d.id lists id number of those who failed
	  irow  <- d.id==i
	  if(p==1){
        term1 <- W[irow]*(d.Xt[irow] - Xbar[irow])
	  } else {
        term1 <- array(W[irow]*(d.Xt[irow,] - Xbar[irow,]), dim=c(p,1))
	  }
    }  
	
    term2 <- rep(0,p)
	for(j in 1:N){
	  if(p==1){
	   Xitj   <- array(d.Xkti[i,j], dim=c(p,1))
	   Xbartj <- array(Xbar[j], dim=c(p,1))
	  } else {
	    Xitj   <- array(d.Xkti[i,j,], dim=c(p,1))
	    Xbartj <- array(Xbar[j,], dim=c(p,1))
      }
      term2  <- term2 - W[j]*C.wt[i,j]*(Xitj - Xbartj)
	}
#	cat("term1: \n")
#	print(term1)
#	cat("term2: \n")
#	print(term2)
    U.tilde.beta[i,] <- term1 + term2
  }
###  END OF AYESH'S RE-CODING OF VARIANCE ESTIMATOR ###

	cat("U.tilde.beta: \n")
	print(U.tilde.beta)
  
	cat("I.tilde.beta: \n")
	print(I.tilde.beta)

  Var.U   <- aTa(U.tilde.beta) 
  J.inv   <- solve(I.tilde.beta) 
  Cov.mat <- J.inv %*% Var.U %*% J.inv  
  se      <- sqrt(diag(Cov.mat))
  rr      <- exp(beta.hat)

  

### 6 - OUTPUT RESULTS

  cat("\n")
  cat("-----------------------------------------------------------------","\n")
  
  cat ("RESULTS SUMMARY: \t estimates converged in ", iter.num, " iterations.")
  cat("\n")
  cat("-----------------------------------------------------------------","\n")

  cat("\n")
  cat(format("Variable",width=10),
      "\tCoeff   \t exp(coef) \t Robust SE \t 95% CI for exp(coef)" )
  cat("\n")
  cat("-----------------------------------------------------------------","\n")
  cat("\n")
  for (i in 1: p) {
    cat(format(names(m)[i+1],width=10), " \t",
        format(round(beta.hat[i], 4), digits=4,width=10), "\t",
        format(rr[i],digits=3,width=10),
        format(se[i],digits=4,width=12),
        "\t (",format(exp(beta.hat[i]-1.96*se[i]),digits=4), ",",
        format(exp(beta.hat[i]+1.96*se[i]),digits=4),")","\n")
  }
  cat("-----------------------------------------------------------------","\n")

  out <- NULL
  out$beta  <- beta.hat
  out$se    <- se
  out$score <- U.tilde.beta
  out$info  <- I.tilde
  out$covb  <- Cov.mat
  return(out)
}

## end of file ##
