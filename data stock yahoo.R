
# load required packages (install them before)

library(tseries)
library(zoo)
library(quantmod)
library(vars)

# load other files for the estimation

source("Bernoulli QMLE.R")
source("MWLSE.R")


# download stock prices

ticker<-c('AAPL', 'MSFT', 'INTC')
nShares<-length(ticker)

start<-'1987-01-01'
end<-'2023-10-11'

prices <- function(ticker, start, end) {
  
  x=get.hist.quote(instrument = ticker,
                   start = start,
                   end = end, quote = "AdjClose",
                   retclass = "zoo")
  
  dimnames(x)[[2]] <- as.character(ticker)
  x
} 

stocks <- do.call(cbind, lapply(ticker, prices, start=start, end=end))
head(stocks)


# convert daily prices to quarterly prices 
z = to.quarterly(stocks, indexAt = "last", OHLC = FALSE)
head(z)
dim(z)


# compute returns

z1 = as.matrix(z[-1,])
z0 = as.matrix(z[-nrow(z),])

returns = log(z1) - log(z0)
head(returns)


# compute return signs

y = (returns>0)+0
View(y)

# time series plots

par(mfrow=c(3,1), mar=c(4.1,4.1,1.1,2.1))
time <- seq(from=1987,to=2023,length.out = nrow(y))
plot(time, returns[,1], type="l", ylab="AAPL", xlab="Time")
plot(time, returns[,2], type="l", ylab="MSFT", xlab="Time")
plot(time, returns[,3], type="l", ylab="INTC", xlab="Time")



### Estimation

# compute starting values from a VAR model in the demeaned series

mu = colMeans(y)
x = t(t(y) - mu)

var = VAR(x, type="none")
Bvar = as.matrix(Acoef(var)[[1]])
Bvar[Bvar<0] = 0.001

a0 = rowMeans(Bvar)

d = ncol(y) 
d0 = (diag(d) - Bvar - diag(a0))%*%mu

b0 = as.vector(Bvar)

x0 = c(d0, a0, b0)
x = log(x0/(1-x0))
mx = length(x)

# Bernoulli QMLE

opt_B = optim(par=x, fn=BerQML, method = "BFGS", y=y,
              control = list(ndeps=rep(1e-5,mx), reltol=1e-8, maxit=500))
opt_B

est_B = exp(opt_B$par)/(1+exp(opt_B$par))
w_B = est_B[1:d]
A_B = diag(est_B[(d+1):(d+d)])
B_B = matrix(est_B[(d+d+1):(d+d+d^2)], d, d)
round(w_B,3)
round(A_B,3)
round(B_B,3)


# compute standard errors and z-tests

hes <- optimHess(par = opt_B$par, fn = BerQML, y=y)
SE_B = sqrt(diag(solve(hes)))
z_B = abs(opt_B$par/SE_B)
z_B

round(est_B,2)
round(SE_B, 2)
round(z_B, 2)


# First step LS

opt_LS = optim(par=x, fn=BerLS, method = "BFGS", y=y,
              control = list(ndeps=rep(1e-5,mx), reltol=1e-8, maxit=500))
opt_LS

est_LS = exp(opt_LS$par)/(1+exp(opt_LS$par))
w_LS = est_LS[1:d]
A_LS = diag(est_LS[(d+1):(d+d)])
B_LS = matrix(est_LS[(d+d+1):(d+d+d^2)], d, d)
round(w_LS,3)
round(A_LS,3)
round(B_LS,3)


# compute standard errors and z-tests

hes <- optimHess(par = opt_LS$par, fn = BerLS, y=y)
SE_LS = sqrt(diag(solve(hes)))
z_LS = abs(opt_LS$par/SE_LS)
z_LS

round(est_LS,2)
round(SE_LS, 2)
round(z_LS, 2)


# Second step WLS

W = W_est(est_LS, y)

opt_WLS = optim(par=x, fn=BerWLS, method = "BFGS", y=y, W=W$W,
               control = list(ndeps=rep(1e-5,mx), reltol=1e-8, maxit=500))
opt_WLS

est_WLS = exp(opt_WLS$par)/(1+exp(opt_WLS$par))
w_WLS = est_WLS[1:d]
A_WLS = diag(est_WLS[(d+1):(d+d)])
B_WLS = matrix(est_WLS[(d+d+1):(d+d+d^2)], d, d)
round(w_WLS,3)
round(A_WLS,3)
round(B_WLS,3)


# compute standard errors and z-tests

hes <- optimHess(par = opt_WLS$par, fn = BerWLS, y=y, W=W$W)
SE_WLS = sqrt(diag(solve(hes)))
z_WLS = abs(opt_WLS$par/SE_WLS)
z_WLS

round(est_WLS,2)
round(SE_WLS, 2)
round(z_WLS, 2)

# compute probability of a positive return

lab = la_est(est_B,y)
lal = la_est(est_LS,y)
lah = la_est(est_WLS,y)


# compute mean absolut errors

mean(abs(t(y)-lab))
mean(abs(t(y)-lal))
mean(abs(t(y)-lah))

rowMeans(abs(t(y)-lab))
rowMeans(abs(t(y)-lal))
rowMeans(abs(t(y)-lah))


# plot probability of a positive return

par(mfrow=c(3,1), mar=c(4.1,4.1,1.1,2.1))
time <- seq(from=1987,to=2023,length.out = nrow(y))
plot(time[-1], lah[1,-1], type="l", col="red", ylab="AAPL", xlab="Time")
plot(time[-1], lah[2,-1], type="l", col="red", ylab="MSFT", xlab="Time")
plot(time[-1], lah[3,-1], type="l", col="red", ylab="INTC", xlab="Time")


# plot acf of Pearson residuals

eh = (t(y) - lah)/sqrt(lah*(1-lah))

par(mfrow=c(3,1), mar=c(4.1,4.1,1.1,2.1))
for(j in 1:d){
  acf(eh[j,-1], lag.max = 25)
}

