
### Functions for multivariate Bernoulli QMLE


### Multivariate Bernoulli quasi-likelihood

BerQML = function(theta, y){
  
  y = t(y)
  d = nrow(y)
  n = ncol(y)
  
  w = exp(theta[1:d])/(1 + exp(theta[1:d]))
  th = theta[-c(1:d)]
  
  A = diag(exp(th[1:d])/(1+exp(th[1:d])))
  th = th[-c(1:d)]
  
  B = matrix(exp(th[1:d^2])/(1+exp(th[1:d^2])), d, d)
  
  
  la = matrix(NA, nrow=d, ncol=n)
  la[,1] = apply(y[,1:10],1,mean)
  
  for (t in 2:n){
    la[,t] = w +  A %*% la[,t-1] + B %*% y[,t-1]
  }
  
  if(all(la<1)==FALSE){
    l = 10e6
  }else{
    
    y0 = y[,-1]
    la0 = la[,-1]
    lpt = diag( t( y0 ) %*% log( la0 ) ) + diag( t( 1 - y0 ) %*% log( 1 - la0 ) )
    l =  - sum( lpt )
  }
  
  l
}


