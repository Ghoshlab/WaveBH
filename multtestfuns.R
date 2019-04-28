#  R functions for wavelet-based multiple testing


# multtestdep.r


#  Created by Debashis Ghosh on 5/4/11.
# Copyright 2011 Penn State. All rights reserved.
# Requires waveslim 

require("waveslim")
bonffun = function(pvals, alpha = 0.05) {
  n = length(pvals)
  crit = alpha/n
  nrej = sum(pvals <= crit)
  
  return(list(nrej=nrej))
  
}

byfun = function(pvals,alpha = 0.05) {
  n = length(pvals)
  pord = sort(pvals)
  crit = (1:n)*alpha/(n*cumsum(1/(1:n)))
  rej = (1:n)[pord <= crit]
  nrej = ifelse(sum(pord <= crit) == 0,0,max(rej))
  
  return(list(nrej=nrej))
  
  
}

fourierbonffun = function(pvals,alpha=0.05) {
  phat = sort(pvals)
  m = length(pvals)
  phat = c(0,phat,1)
  tmp = diff(phat,lag=1)
  n = length(tmp)
  
  # apply fourier transform 
  #	m = (log(n,2))
  #	pad = ifelse(ceiling(m) != floor(m),2^(ceiling(m)) - n,0)
  #    evec = c(vec,rep(0,pad))
  
  tmp.fft = fft(tmp)/sqrt(n)
  tmp.i = sqrt(as.numeric(tmp.fft*Conj(tmp.fft)))
  # These components are approximately independent, so get
  # p-values using chi-squared for the 
  
  tmp.pval = 1-pchisq(tmp.i,2)
  nrej = sum(tmp.pval <= alpha/n)
  
  return(list(nrej=nrej,tmp.pval=tmp.pval))
  
  
}


fourierbhfun = function(pvals,alpha=0.05) {
  phat = sort(pvals)
  m = length(pvals)
  phat = c(0,phat,1)
  tmp = diff(phat,lag=1)
  n = length(tmp)
  
  # apply fourier transform 
  #	m = (log(n,2))
  #	pad = ifelse(ceiling(m) != floor(m),2^(ceiling(m)) - n,0)
  #    evec = c(vec,rep(0,pad))
  
  tmp.fft = fft(tmp)/sqrt(n)
  tmp.i = sqrt(as.numeric(tmp.fft*Conj(tmp.fft)))
  # These components are approximately independent, so get
  # p-values using chi-squared for the 
  
  tmp.pval = 1-pchisq(tmp.i,2)
  crit = (1:n)*alpha/n
  nrej = sum(sort(tmp.pval) <= crit)
  
  return(list(nrej=nrej,tmp.pval=tmp.pval))
  
  
}

wavebhfun =  function(pvals,alpha=0.05) {
  
  phat = sort(pvals)
  m = length(pvals)
  phat = c(0,phat,1)
  tmp = diff(phat,lag=1)
  n = length(tmp)
  vec = cumsum(tmp)/(1:n)
  crit = 1/n
  
  # apply wavelet transform 
  #	m = (log(n,2))
  #	pad = ifelse(ceiling(m) != floor(m),2^(ceiling(m)) - n,0)
  #    evec = c(vec,rep(0,pad))
  
  tvec = mra(tmp,method="modwt",J=1)
  mvec = cumsum(tvec$S1)/(1:n)
  
  prej = (1:n)[mvec <= alpha*crit]
  nrej = ifelse(sum(mvec <= alpha*crit) == 0,0,max(prej,na.rm=T)-1)
  #hrej = ifelse(nrej == 0, NA, order(pvals)[1:nrej])
  
  return(list(nrej=nrej,crit=crit))
}

waveadabhfun = function(pvals,q=0.05,alpha=0.05,levels=1) {
  
  phat = sort(pvals)
  # m is number of hypotheses
  m = length(pvals)
  
  pord = c(0,phat,1)
  tmp = diff(pord,lag=1)
  # n is length of spacings, n=m+1
  n = length(tmp)
  vec = cumsum(tmp)/(1:n)
  crit = 1/n
  
  # apply wavelet transform 
  #	m = (log(n,2))
  #	pad = ifelse(ceiling(m) != floor(m),2^(ceiling(m)) - n,0)
  #    evec = c(vec,rep(0,pad))
  
  tvec = mra(tmp,method="modwt",J=1)
  mvec = cumsum(tvec$S1)/(1:n)
  
  # First stage: apply BH at q/(1+q)
  
  st1rej = (1:n)[mvec <= crit*q/(1+q)]
  nst1rej = ifelse(sum(mvec <= crit*q/(1+q)) == 0,0,max(st1rej,na.rm=T)-1)
  # Estimate n0
  n0 = m-nst1rej
  if (n0 > 0) {
    # Second stage: apply BH at alpha with n0
    st2rej = (1:n)[mvec <= alpha/n0]
    nst2rej = ifelse(sum(mvec <= alpha/n0) == 0,0,max(st2rej,na.rm=T)-1)
  }
  else
  {
    nst2rej = 0
  }
  #hrej = ifelse(nrej == 0, NA, order(pvals)[1:nrej])
  
  return(list(nst1rej=nst1rej,n0=n0,nst2rej=nst2rej))
}
