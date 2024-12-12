#######################################################################################
###
### File name : 02_HMM_HMM_functions.R
### Program developer: aardepi1, Piet Aarden based on work of all below:
###                    H. Ganjgahi (Oxford BDI, habib.ganjgahi@bdi.ox.ac.uk)
### Date: 10-03-2023
### Project/Trial code: NO.MS - HMM
### Description: Function script which accompanies 02_hmm/* scripts
### Application: R 4.1.0
### Libraries used: tidyverse_1.3.2, haven_2.5.1, reshape2_1.4.4
### Source table locations: /vob/CNOMSANON/anon/anon_2/anon_source/
###                         /funstorage/NSGDD_MS_MRI/HMM/00_source_data/
### Input: 
### Output: 
###
#######################################################################################

library(mclust)
library(MCMCpack)
library(Matrix)
library(expm) 
library(MASS)
library(mvtnorm)
library(msm)
library(cluster)

Init_clara = function(y,K,Metric){
  clara.res <- clara(y, K, samples = 50, pamLike = TRUE,metric = Metric)
  #ll          = order(mod1$centers)
  ll          = clara.res
  init        = list()
  init$mu     = t(ll$medoids)
  #init$sigma  = matrix(0,ncol(y),K)
  init$sigma	   = array(0,c(ncol(y),ncol(y),K))
  for (i in 1:K) {
    init$sigma[,,i]=cov(y[ll$clustering==i,])
  }
  init$A  = t(replicate(K, runif_simplex(K)))
  init$pi = runif_simplex(K)
  return(init)
}

Init_clara2 = function(y,K,USUBJID,Time,Metric){
  clara.res <- clara(y, K, samples = 50, pamLike = TRUE,metric = Metric)
  #ll          = order(mod1$centers)
  ll=clara.res
  init        = list()
  init$mu     = t(ll$medoids)
  #init$sigma  = matrix(0,ncol(y),K)
  init$sigma	   = array(0,c(ncol(y),ncol(y),K))
  for (i in 1:K) {
    init$sigma[,,i]=cov(y[ll$clustering==i,])
  }
  a=statetable.msm(ll$clustering,USUBJID)
  init$A  =a/rowSums(a)
  init$pi = table(ll$clustering[Time==0])/length(ll$clustering[Time==0])
  #init$A  = t(replicate(K, runif_simplex(K)))
  #init$pi = runif_simplex(K)
  return(init)
}


InitKM_mv <- function(y, K){
  mod1 <- kmeans(y, K, nstart = 50, iter.max = 400, algorithm="MacQueen")
  ll = mod1
  init <- list()
  init$mu <- t(mod1$centers)
  init$sigma <- array(0, c(ncol(y), ncol(y), K))
  for (i in 1:K) {
    init$sigma[ , , i] <- cov(y[mod1$cluster==i, ])
  }
  init$A <- t(replicate(K, runif_simplex(K)))
  init$pi <- runif_simplex(K)
  return(init)
}

InitKM_mv2 = function(y,K,USUBJID,Time){
  mod1        = kmeans(y,K,nstart = 20,iter.max = 30)
  while (!mod1$ifault==0) {mod1        = kmeans(y,K,nstart = 20,iter.max = 30,algorithm="MacQueen") }
  #ll          = order(mod1$centers)
  ll=mod1
  init        = list()
  init$mu     = t(mod1$centers)
  #init$sigma  = matrix(0,ncol(y),K)
  init$sigma       = array(0,c(ncol(y),ncol(y),K))
  for (i in 1:K) {
    init$sigma[,,i]=cov(y[mod1$cluster==i,])
  }
  a=statetable.msm(mod1$cluster,USUBJID)
  init$A  =a/rowSums(a)
  init$pi = table(mod1$cluster[Time==0])/length(mod1$cluster[Time==0])
  #init$A  = t(replicate(K, runif_simplex(K)))
  #init$pi = runif_simplex(K)
  return(init)
}

hmmfilter_dtApp<-function(pi, A, emission){
  
  #initialisation
  K <- ncol(A)
  N <- nrow(emission)
  alpha <- matrix(0, K, N)
  scl <- vector('numeric', N)
  beta <- matrix(0, K, N)
  emission <- t(emission)
  
  #Forward algorithm log p(z_t = j | x_{1:t}) 
  #logalpha[1,] = log(pi) + dnorm(y[1], mu, sigma,log=TRUE)
  alpha[, 1] <- pi * emission[, 1]
  scl[1] <- sum(alpha[, 1])
  if(scl[1] == 0){scl[1] = 1}
  alpha[, 1] <- alpha[,1] / scl[1]
  
  for (t in 2:N) {
    #make sure A[,,t] is the correct emission probability (double check the derivation)
    alpha[ , t] <- (t(A[ , , t]) %*% alpha[ , t - 1]) * emission[ , t]
    scl[t] <- sum(alpha[ , t])
    if(scl[t] == 0){scl[t] = 1}
    alpha[ , t] <- alpha[ , t] / scl[t]
  }
  
  loglik <- sum(log(scl + 1e-8))
  beta[ , N] <- 1
  for (t in seq(N-1, 1, -1)) {
    #make sure A[,,t+1] to make sure it's the correct emission probability (double check the derivation)
    beta[ , t] <- A[ , , t+1] %*% (beta[ , t+1] * emission[ , t + 1])
    dd <- sum(beta[ , t])
    if(dd == 0){dd = 1}
    beta[ , t] = beta[ , t] / dd
  }
  
  return(list(alpha = t(alpha), beta = t(beta), scl = scl, loglik = loglik))
}

dthmmTrans<-function(deltaT,A){
  nu	  = length(deltaT)
  K	  = nrow(A)
  Ptdelta = array(0,c(K,K,nu))
  for ( t in 1:nu) {
    Ptdelta[,,t] = as.matrix(A%^%deltaT[t])
  }
  return(Ptdelta)
}

hmmtwoslice_dtApp <- function(A, alpha, beta, emission){
  
  K <- ncol(emission)
  N <- nrow(emission)
  xi <- array(0, c(K, K, N - 1))
  xi_summed <- matrix(0, K, K)
  
  for (t in (N - 1):1){
    b <- beta[t + 1, ] * emission[t + 1, ]
    
    #make sure A[,,t+1] is the correct transition matrix (double check the derivation)
    tmpXi <- A[ , , t + 1] * (alpha[t, ] %*% t(b))
    kk <- sum(tmpXi)
    if(kk == 0){kk = 1}
    xi_summed <- xi_summed + tmpXi / kk # inlined call to normalize
    xi[ , , t] <- tmpXi
  }
  
  return(list(xi_sum = xi_summed, xi = xi))
}

runif_simplex <- function(T) {
  x <- -log(runif(T))
  x / sum(x)
}


ctdthmm_viterbi <- function(pi, Ptdelta, emission){
  
  K <- ncol(emission)
  N <- nrow(emission)
  emission <- t(emission)
  delta <- matrix(0, K, N)
  psi <- matrix(0, K, N)
  path <- matrix(0, 1, N)
  t <- 1 
  delta[ , t] <- (pi * emission[ , t]) / sum(pi * emission[ , t])
  
  for (t in 2:N) {
    for (j in 1:K) {
      ind <- which.max(delta[ , t - 1] * Ptdelta[ , j, t])
      psi[j, t] <- ind
      delta[j, t] <- max(delta[, t - 1] * Ptdelta[ , j, t]) * emission[j, t]
    }
    delta[, t] <- delta[, t] / sum(delta[, t])
  }
  
  #Traceback
  path[N] <- which.max(delta[ , N])
  for (t in (N - 1):1) {
    path[t] <- psi[path[t + 1], t + 1]
  }
  
  return(path)
}

ctdthmm_MultSubj_viterbi <- function(hmmDT, y, seq, Time){
  
  A <- hmmDT$A
  pi <- hmmDT$pi
  mu <- hmmDT$mu
  sigma <- hmmDT$sigma
  
  nT <- nrow(y)
  nS <- length(seq)
  K <- nrow(A)
  seqidx <- cumsum(c(1, seq))
  path <- vector('numeric', nT)
  
  #calculate emission
  emission <- matrix(0, nT, K)
  for (k in 1:K) {emission[ , k] <- dmvnorm(y, mean = mu[ , k], sigma = sigma[ , , k])}
  
  #Find the most-probable (Viterbi) path 
  for (n in 1:nS) {
    ndx <- seqidx[n]:(seqidx[ n + 1] - 1)
    if (length(ndx) > 1) {
      deltaT <- Time[ndx]
      Ptdelta <- dthmmTrans(deltaT, A)
      Bi <- emission[ndx, ]
      path[ndx] <- ctdthmm_viterbi(pi, Ptdelta, Bi)
    }
  }
  
  return(path)
}


hmm_DTMultSubj_MultV_Bayes<-function(y, init, seq, K, Time, thresh){
  
  #set the initial values
  mu <- init$mu
  sigma <- init$sigma
  pi <- init$pi
  A <- init$A
  
  nT <- nrow(y)
  nS <- length(seq)
  P <- ncol(y)
  seqidx <- cumsum(c(1, seq));
  weights <- matrix(0, nT, K)
  
  #prior
  b0 <- apply(y, 2, median)
  R <- c()
  for (iii in 1:P) {
    R <- c(R, (diff(range(y[ , iii]))^2))
  }
  B0 <- diag(1 / R)
  c0 <- 2.5 + (P - 1) / 2
  g0 <- 0.5 + (P - 1) / 2
  G0 <- (100 * g0 / c0) * diag(1 / R)
  C0 <- diag(P)
  
  #find NA indices in each column
  IND <- matrix(FALSE, nT, P)
  for (i in 1:P) {
    IND[ , i] <- !is.na(y[ , i])
  }
  
  conv <- TRUE
  i <- 0
  loglike <- c()
  logPopst <- c()
  while (conv) {
    i = i + 1
    print(i)
    #E-step
    #calculate the emission for the EM algorithm
    emission <- matrix(0, nT, K)
    for (k in 1:K) {emission[ , k] = dmvnorm(y, mean = mu[ , k], sigma = sigma[ , , k], log = TRUE)}
    xi <- array(0, c(K, K, nT))
    ENk1 <- 0
    ENjk <- 0
    LogL <- 0
    for (n in 1:nS) {
      ndx <- seqidx[n]:(seqidx[n + 1] - 1)
      if(length(ndx) > 1){
        Bi <- emission[ndx, ]
        deltaT <- Time[ndx]
        Ptdelta <- dthmmTrans(deltaT, A)
        
        E_out <- hmmfilter_dtApp(pi, Ptdelta, exp(Bi))
        LogL <- LogL + E_out$loglik
        alpha <- E_out$alpha
        beta <- E_out$beta 
        gamma <- alpha * beta
        
        z <- rowSums(gamma)
        z[z == 0] = 1
        #replace the for loop with rep.row
        #gamma  = gamma/rep.row(z,n)
        for(t in 1:length(ndx)){
          gamma[t, ] <- gamma[t, ] / z[t]
        }
        weights[ndx, ] <- gamma
        E_out2 <- hmmtwoslice_dtApp(Ptdelta, alpha, beta, exp(Bi))
        #initial state at t=0
        ENk1 <- ENk1 + gamma[1, ]
        #number of times moves from state j to k
        ENjk <- ENjk + E_out2$xi_sum
      }
    }
    loglike <- c(loglike, LogL)
    #multivariate M-step
    #update state mean
    m <- matrix(0, P, K)
    dd <- colSums(weights)
    mm <- t(weights) %*% y
    logPrior <- 0
    for(k in 1:K){
      mm2 <- solve(dd[k] * solve(sigma[ , , k]) + B0)
      mm3 <- (mm[k, ] %*% solve(sigma[ , , k]) + b0 %*% B0) %*% mm2
      m[, k] <- mm3
    }
    logPrior <- logPrior + sum(dmvnorm(t(m), mean = b0, sigma = solve(B0), log = TRUE))

    si <- array(0, c(P, P, K))
    for(k in 1:K){
      for (nn in 1:nT) {
        si[ , , k] = si[ , , k] + weights[nn, k] * ((y[nn, ] - m[ , k]) %*% t(y[nn, ] - m[ , k]))
      }
      Ck <- si[ , , k] / 2 + C0
      si[,,k] <- (si[ , , k] + 2 * solve(C0))/(2 * c0 + dd[k] + (P + 1))
    }
    print(max(abs(mu - m)))
    ff <- 0
    for (k in 1:K) {
      ff <- ff + solve(sigma[ , , k])
    }
    C0 <- (ff + G0) / (g0 + K * c0 - (P + 1) / 2)
    logl <- 0
    for (k in 1:K) {
      logl <- logl - (c0 + (P + 1) / 2) * log(det(si[ , , k])) - sum(diag(solve(C0) %*% solve(si[ , , k])))
    }
    logl <- logl - (g0 + K * c0 - (P + 1) / 2) * log(det(C0)) - sum(diag(solve(C0) %*% G0))
    logPrior <- logPrior + logl
    mu <- m
    print(max(abs(sigma - si)))
    sigma <- si
    ENjk <- ENjk + 1 / K - 1
    ENjk[ENjk < 0] <- .00001
    print(max(abs(A -(ENjk /rowSums(ENjk)))))
    A <- ENjk / rowSums(ENjk)
    for (k in 1:K) {
      logPrior <- logPrior + log(ddirichlet(A[k, ], rep(1 / K, K)))
    }
    
    #pi with drichlet weights
    alpha <- 1 / K
    print(max(abs(pi - ((ENk1 - 1 + alpha) / sum(ENk1 - 1 + alpha)))))
    pi <- (ENk1 - 1 + alpha)
    pi[pi < 0] <- .000001
    pi <- pi / sum(pi)
    logPrior <- logPrior + log(ddirichlet(pi, rep(1 / K, K)))
    logPopst <- c(logPopst, logPrior + LogL)
    
    ###evaluate convergence based on likelihood
    if(i > 1){
      thr <- loglike[i] - loglike[i - 1]
      thr2 <- logPopst[i] - logPopst[i - 1]
      thr3 <- (loglike[i] - loglike[i - 1])/loglike[i - 1]
      if(abs(thr) < thresh){conv = FALSE}
      print(list(thr3=thr3,thr2 = thr2,thr = thr))
    }
    print(list(mu = mu,sA = A,pi = pi))
  }
  return(list(sigma = sigma, mu = mu, pi = pi, A = A, loglik = loglike, logPost = logPopst, C0 = solve(C0)))
}


