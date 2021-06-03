aaT <- function(a) {return(a %*% t(a))} 
aTa <- function(a) {return(t(a) %*% a)} 


S.tilde <- function(b, ftime, X, start1, stop1, weight1, get.s2=FALSE)
{
  ##-######################################################
  ## this function computes the s.tilde terms needed to   #
  ## compute the score vector and information matrix      #
  ## input:  1xp beta vector        b                     #
  ##         nxp covariate matrix   X                     #
  ##         nx1 weight vector      wt                    #
  ## output: s0, s1, s2, as in MSCM paper                 #
  ##-######################################################
  
  n       <- length(start1)
  p       <- ncol(X)

  mybeta  <- array(b, dim=c(1,p))
  myweight<- rep(0,n)
  y       <- rep(0,n)

  for (j in 1:n){
    y[j] <- 0
    if ((start1[j]<ftime) && (ftime <=stop1[j])) {
      y[j]<-1
      myweight[j]<-weight1[j]
    }
  }

  X.t<-array(dim=c(n,p))
  for (k in 1:n){
    X.t[k,]<- X[k,]
  }

  w.Y.e.betaX <- array(myweight*y* exp(mybeta %*% t(X.t)), dim=c(1,n))
  out         <- NULL
  out$s0      <- sum(w.Y.e.betaX)                       # a number (scalar)
  out$s1      <- array(w.Y.e.betaX %*% X.t, dim=c(p,1)) # px1 vector

  if(get.s2) {
    out$s2 <- array(0, dim=c(p,p))                   # pxp matrix
    for (k in 1:n) {
      out$s2 <- out$s2 + w.Y.e.betaX[k]*(X.t[k,]%*%t(X.t[k,]))
    }
    out$s2 <- array(out$s2, dim=c(p,p))
  }  

  out$s0      <- as.real(out$s0)
  out$s1      <- as.matrix(out$s1,nrow=n,ncol=1) # p*1 vector

  return(out)
}


getLLkhd <- function(b, d.Xt, d.Xkti, d.Ykti, d.wkti, W, nsubject, N, p)
{
	## This function computes the loglikelihood value
    ## given the current beta estimate, and the covariate
    ## information at the failure times.
	
	## This function returns the log-lkhd value as well as
	## matrix C, which is needed to update the score and information 
	## matrix.
	
	
    ## notes to me
    B <- array(dim = c(N,1)) 
    C <- array(dim = c(nsubject,N))  

    expval <- array(dim = c(nsubject, N)) # note to me
    i <- 1
    for (s in nsubject) {
      for (j in 1:N) {
        if (p==1) {
          expval[i,j] <- d.Xkti[i,j] * b
        } else {
          expval[i,j] <- array(d.Xkti[i,j,],dim=c(1,p)) %*% b
        }
      }
      i <- i+1
    }
    
    ## calculate log-likelihood	
    B <- d.Xt %*% b      # N x 1 vector
    C <- d.Ykti * d.wkti * exp( expval )
    A <- B - log(colSums(C)); # sum over k

	out    <- NULL
    out$LL <- sum(W * A);   # sum over i
	out$C  <- C

	out$LL <- as.real(out$LL)
	out$C  <- as.matrix(out$C, nrow=nsuject, ncol=N)
    return(out)
}

ui.w <- function(mybeta, Xmat, event, t0, t1, id, wts, wkti)
{
  ##-######################################################
  ## this function computes the weighted score residuals  #
  ## matrix.  It is needed for variance estimation in     #
  ## an MSCM analysis.                                    # 
  ##-######################################################

  p      <- ncol(Xmat)
  id.uni <- unique(id)     # subject id's 
  nsub   <- length(id.uni)       # number subjects
  rname  <- 1:nrow(Xmat)         # index rows of X

  delta <- rep(0,nsub)           # failure/censored indicator
  wt.ft <- rep(0,nsub)           # wt at failure time
  f.t   <- rep(0,nsub)           # follow-up time
  X.t   <- array(0,dim=c(nsub,p))# covars at failure times

  for (j in 1:nsub) {           # for each subject
    pt    <- id.uni[j]          # get id
    r.pt  <- rname[id==pt]      # rows of X with their obs
    t1.pt <- t1[id==pt]         # end time of intervals
    wt.pt <- wts[id==pt]        # weights
    Y.pt  <- event[id==pt]      # response (failed or not)
    nobs  <- length(t1.pt)      # number visits on subject
    X.pt  <- array(Xmat[id == pt,], dim=c(nobs,p)) 
    
    f.t[j]  <- t1.pt[nobs]      # time at end of fol-up
    wt.ft[j]<- wt.pt[nobs]      # wt at end of fol-up
    X.t[j,] <- X.pt[nobs,]      # covars at end of fol-up

    if (any(Y.pt==1)) delta[j]<-1 # failed or censored?
  }


  U.beta <- matrix(0,nrow=p,ncol=nsub)
  id.uni <- unique(id)          # subject id's unsorted

  cat("Calculating weighted score residuals matrix. It may take a while.\n")
#  cat("k: ")
  for (k in 1:nsub) {
    cat(k)

    S.k  <- S.tilde(mybeta,f.t[k],Xmat,t0,t1,wts)
    S0.k <- as.numeric(S.k$s0)
    S1.k <- array(S.k$s1, c(p,1))

    U.beta[,k] <- delta[k]* (array(X.t[k,],c(p,1)) - S1.k/S0.k)
    term1      <- matrix(0, nrow=p,ncol=1)
    term2      <- matrix(0, nrow=p,ncol=1)

    for (j in 1:nsub) {
      
      Yk.tj <- ifelse(f.t[k]>=f.t[j], 1, 0)

      # if subj. j is a failure and in risk set at t_k
      if ((delta[j]*Yk.tj) == 1 )
      {  
        pt    <- id.uni[k]
        t0.k  <- t0[id == pt]
        t1.k  <- t1[id == pt]
        index <- rank(t1.k)[t1.k == min((t1.k[(t1.k >= f.t[j])]))]
        const <- wt.ft[j]
        Xk.tj <- matrix(0, p, 1)
        
        if(t0.k[index] < f.t[j]){
          temp <- matrix(Xmat[id == pt,])
          if(length(t0.k) > 1){
            Xk.tj  <- array(matrix(temp[index,]), c(p,1)) 
          } else  Xk.tj <- temp
        }
        comment <- "
in original paper, we have factor
w_j(t_i)

but this should read (in corrected paper)
w_j(t_j) w_i(t_j)

so we make the correction here, where in this program
we used:
  wt.ft[j] <- weight of subject j at failure time of subject j

this was a mistake, but a correct one.

in this prog,
  paper i -> prog k
  paper j -> prog j
so
  paper w_j(t_j) * w_i(t_j) ->
   prog     wjtj * wktj

we already have wjtj in wt.ft[j]

now make wktj[k,j]

";
        S.j  <- S.tilde(mybeta,f.t[j],Xmat,t0,t1,wts)
        S0.j  <- as.numeric(S.j$s0)
        S1.j  <- array(S.j$s1, c(p,1))

        const <- as.real(wkti[k,j]*wt.ft[j]*exp(t(mybeta) %*% Xk.tj))
        term1 <- term1 - const*Xk.tj/S0.j
        term2 <- term2 + const*S1.j/(S0.j^2)
        
      } # end of if delta stmt
    } # end of j loop
    U.beta[,k] <- U.beta[,k]+term1+term2
  } # end for k
  cat("\n")
  U.beta[,k]

  return(U.beta)
}
