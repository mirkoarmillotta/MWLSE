

### Functions Two-stage WLSE of binary data


### First stage Least Squares (LS)

BerLS = function(theta, y){
  
  y = t(y)
  d=nrow(y)
  n=ncol(y)
  
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
    l=10e6
  }else{
    
    e0 = y[,-1] - la[,-1]
    lpt = diag( t( e0 ) %*% e0  ) 
    l =  sum( lpt )
  }
  
  l
}



###  Estimate optimal weights using Pearson's residuals mean moment estimation

W_est <- function(theta, y){
  
  y = t(y)
  d=nrow(y)
  n=ncol(y)
  
  w = theta[1:d]
  th = theta[-c(1:d)]
  A = diag(th[1:d])
  th = th[-c(1:d)]
  B = matrix(th[1:d^2], d, d)
  
  la = matrix(NA, nrow=d, ncol=n)
  la[,1] = apply(y[,1:10],1,mean)
  var = matrix(NA, nrow=d, ncol=n)
  var[,1] = la[,1] * ( 1 - la[,1] )
  K = matrix(0, nrow=d, ncol=d)
  
  for (t in 2:n){
    la[,t] = w +  A %*% la[,t-1] + B %*% y[,t-1]
    var[,t] = la[,t] * ( 1 - la[,t] )
    K = K + (y[,t] - la[,t])%*%t(y[,t] - la[,t])
  }
  
  k = diag(K)
  S = diag(1/sqrt(k))
  R = S %*% K %*% S
  p = mean(R[lower.tri(R)])
  
  a <- (1+(d-2)*p)/((1-p)*(1+(d-1)*p))
  b <- -p/((1-p)*(1+(d-1)*p))
  J <- matrix(1,d,d)
  P <- (a-b)*diag(d)+b*J
  
  W = list()
  for(t in 2:n){
    W[[t]] = diag(1/sqrt(var[,t])) %*% P %*%  diag(1/sqrt(var[,t]))
  }
  
  return(list(lambda = la, W = W, var = var, p = p, P = P))
 
}


### Second stage Weighted Least Squares (WLS)

BerWLS = function(theta, y, W){
  
  y = t(y)
  d=nrow(y)
  n=ncol(y)
  
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
    l=10e6
  }else{
    
    l = 0
    e = y - la
    for(t in 2:n){
      
      l = l + t(e[,t]) %*%  W[[t]] %*% e[,t]
    }
  }
  
  l
}


## Compute conditional means (lambdas)

la_est <- function(theta, y){
  
  y = t(y)
  d=nrow(y)
  n=ncol(y)
  
  w = theta[1:d]
  th = theta[-c(1:d)]
  A = diag(th[1:d])
  th = th[-c(1:d)]
  B = matrix(th[1:d^2], d, d)
  
  la = matrix(NA, nrow=d, ncol=n)
  la[,1] = apply(y[,1:10],1,mean)
  
  for (t in 2:n){
    la[,t] = w +  A %*% la[,t-1] + B %*% y[,t-1]
  }
  
  la
  
}
