# Outros

#################################################
#
# plotting binomial and bernoulli likelihoods
# and log-likelihoods
#
# from Aug 23, 2012 Ohio poll:
#
#           est'd n's    reported %'s
# Obama:          418              49
# Romney:         390              46
# Subtotal:       808
# Other:           39               4 ( =  2 "other" + 2 "don't know")
#                 -------------------
# Total:          847 "likely voters"

### from Sept 12-14 2016 Suffolk Poll
### 500 respondents, 401 expressed Clinton or Trump preference.
### 208 for Trump, 193 for Clinton

n <- 401
np <- 208
nn <- 193

binom.lik <- function(x) {exp(lchoose(n,np)+np*log(x)+nn*log(1-x))}
# just calculating directly underflows.  So calculate log and then exp()

log.binom.lik <- function(x) {lchoose(n,np)+np*log(x)+nn*log(1-x)}

bern.lik <- function(x) {x^np*(1-x)^nn}
# this function underflows to zero, so need a cleverer plotting
# scheme...
#
# so, plot the binomial likelihood (same shape, no underflow), and
# adjust the labels on the tick marks to match the bernoulli...

log.bern.lik <- function(x) {np*log(x)+nn*log(1-x)}



par(mfcol=c(2,2))

curve(binom.lik(x),xlab="p [parameter]",
      ylab=expression(L[bin](p)))
#     ylab="Binomial Likelihood: P[data|p]")
title("Binomial Likelihood")

curve(binom.lik(x),xlab="p [parameter]",
      ylab=expression(L[bin](p)),axes=F)
#     ylab="Bernoulli Likelihood: P[data|p]",axes=F)
axis(1)
y.max <- optimize(binom.lik,lower=0,upper=1,maximum=T)$maximum/choose(n,np)
ys <- ((0:4)/4)*y.max
expos <- c(0,floor(log10(ys[-1])))
ys <- round(ys/(10^expos),2)*(10^expos) # round mantissas to 2 places
axis(2,at=seq(0,.025,length=5),labels=
       paste(ys))
box()
title("Bernouli Likelihood")

curve(log.binom.lik(x),xlab="p [parameter]",
      ylab=expression(log(L[bin](p))))
#     ylab="Binomial log-Likelihood: log P[data|p]")
title("Binomial Log-likelihood")

curve(log.bern.lik(x),xlab="p [parameter]",
      ylab=expression(log(L[ber](p))))
#     ylab="Bernoulli log-Likelihood: log P[data|p]")
title("Bernouli Log-likelihood")

##############################################
#
# point estimate and CI for 2004 Bush/Kerry poll in Ohio
#

(phat <- np/n)
# [1] 0.5187032

(se <- sqrt(phat*(1-phat)/n))
# [1] 0.02495133

round(phat + c(-2,2)*se,2)
# [1] 0.47 0.57

###############################################
#
# Bayes' Rule and Terrorists
#

10000/(561.9*10^6)
# [1] 1.779676e-05

0.9999*(1.779676e-05) + .0001*(1- 1.779676e-05 )
# [1] 0.0001177932

0.9999*1.78e-15 /0.0001177932
# [1] 1.510972e-11

(.9999)*(1.5E-11) + (1-.9999)*(1-1.5E-11)
# [1] 0.0001000000

(0.9999)*(1.5E-11)/(0.0001)
# [1] 1.49985e-07


###############################################
#
# Graphs of beta densities
#

par(mfrow=c(4,4))

avals <- c(0.5,1,2,4)
bvals <- c(0.5,1,2,4)

for (alpha in avals) {
  for (beta in bvals) {
    myfun <- function(x) dbeta(x,alpha,beta)
    curve(myfun(x),xlab="p",ylab=paste("dbeta(p,",alpha,",",beta,")",sep=""))
  }
}

################################################
#
# posterior with prior = beta(a=1,b=1)
#

prior <- function(x) {dbeta(x,1,1)}

lik <- function(x) {dbeta(x,np+1,nn+1)}
# binomial likelihood normalized to be a density

posterior <- function(x) {dbeta(x,np+1,nn+1)}

par(mfrow=c(1,1))

curve(posterior(x),from=0.4,to=0.6,n=1000,lty=1,xlab="p",ylab="")
curve(prior(x),from=0.4,to=0.6,n=1000,lty=3,add=T)
curve(lik(x),from=0.4,to=0.6,n=1000,lty=2,add=T)

legend(0.40,15,lty=c(2,1),legend=c("Prior","Posterior = Likelihood"))


#################################################
#
# prev survey had 942 prefer Obama and 1008 prefer Romney
# so set alpha = 942, beta = 1008
# posterior with prior = beta(a=942,b=1008)
#


prior <- function(x) {dbeta(x,942,1008)}

lik <- function(x) {dbeta(x,np+1,nn+1)}
# binomial likelihood normalized to be a density

posterior <- function(x) {dbeta(x,942+np,1008+nn)}

par(mfrow=c(1,1))

curve(posterior(x),from=0.4,to=0.6,n=1000,lty=1,xlab="p",ylab="")
curve(prior(x),from=0.4,to=0.6,n=1000,lty=3,add=T)
curve(lik(x),from=0.4,to=0.6,n=1000,lty=2,add=T)

legend(0.40,20,lty=c(3,2,1),legend=c("Prior","Likelihood","Posterior"))

#################################################
#
# simulation-based equal-tailed interval...
#

nsim <- 10000
p <- rbeta(nsim,942+np,1008+nn)
quantile(p,c(0.025,.5,.975))

########################################################################################
###########################################################
#
# Example 1
#
###########################################################	
library(tmvtnorm)
# Draw from multi-t distribution without truncation
X1 <- rtmvt(n=10000, mean=rep(0, 2), df=2)
X2 <- rtmvt(n=10000, mean=rep(0, 2), df=2, lower=c(-1,-1), upper=c(1,1))

###########################################################
#
# Example 2
#
###########################################################	

df = 2
mu = c(1,1,1)
sigma = matrix(c(  1, 0.5, 0.5,
                   0.5,   1, 0.5,
                   0.5, 0.5,   1), 3, 3)
lower = c(-2,-2,-2)
upper = c(2, 2, 2)

# Rejection sampling
X1 <- rtmvt(n=10000, mu, sigma, df, lower, upper)

# Gibbs sampling without thinning
X2 <- rtmvt(n=10000, mu, sigma, df, lower, upper, 
            algorithm="gibbs")

# Gibbs sampling with thinning
X3 <- rtmvt(n=10000, mu, sigma, df, lower, upper, 
            algorithm="gibbs", thinning=2)	

plot(density(X1[,1], from=lower[1], to=upper[1]), col="red", lwd=2,
     main="Gibbs vs. Rejection")
lines(density(X2[,1], from=lower[1], to=upper[1]), col="blue", lwd=2)
legend("topleft",legend=c("Rejection Sampling","Gibbs Sampling"), 
       col=c("red","blue"), lwd=2)

acf(X1)  # no autocorrelation in Rejection sampling
acf(X2)  # strong autocorrelation of Gibbs samples
acf(X3)  # reduced autocorrelation of Gibbs samples after thinning	
