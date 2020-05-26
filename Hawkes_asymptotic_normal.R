######  Check the asymptotic normality for MLE ###### 

params <- c(0.5, 1, 1.3, 0) #(mu, alpha, beta, x=0)

############################################################################
############################################################################

# Check the asymptotic normality of MLE.
# params : list(mu, alpha, beta)
# given the inverse of the Fisher Info matrix as g (see another code)

############################################################################
############################################################################
############################################################################
############################################################################

set.seed(0)

n <- 3000 # num of MC
t_max <- 300

############################################################################
############################################################################

## Log-Likelihood function of exp Hawkes
loglik <- function(params, t_max, arrivals){
  mu <- params[1]
  alpha <- params[2]
  beta <- params[3]
  omega <- arrivals[2:(length(arrivals)-1)]
  
  f <- function(z) {
    sum(exp(-beta*(omega[z] - omega[1:(z-1)])))
  }
  
  Bi <- c(0, sapply(2:length(omega), f))
  term_1 <- sum(log(mu + alpha*Bi))
  term_2 <- mu*t_max 
  term_3 <- -sum(alpha/beta*(exp(-beta*(t_max - omega)) - 1))
  
  return(term_1 - (term_2 + term_3))
}

############################################################################
############################################################################

## MC simulation
#set.seed(123)

# list of MLE
mu<- rep(0,n)
alpha<- rep(0,n)
beta<- rep(0,n)

# Simulate MLE
# simulate_uni_hawkes (defined in other code)
for(i in 1:n){
  t <- simulate_uni_hawkes(params, t_max)
  mle <- optim(params[1:3], loglik, arrivals = t, t_max=t_max, method = "L-BFGS-B", control=list(fnscale=-1),
               lower=c(0.001,  0, 0.001), upper=c(100, 100, 100))
  mu[i]<-mle$par[1]
  alpha[i]<-mle$par[2]
  beta[i]<-mle$par[3]
  print(i/n*100)
}

# error of MLE
res_mu0<- (mu-params[1])*sqrt(t_max)
res_alpha0<- (alpha-params[2])*sqrt(t_max) 
res_beta0<- (beta-params[3])*sqrt(t_max) 

############################################################################
############################################################################

## plot
par(mfrow=c(1,3))
  
hist(res_mu0, freq=F, breaks=c(-Inf, seq(-5*sqrt(g[1,1]), 5*sqrt(g[1,1]), 0.5), Inf), ylim=c(0,dnorm(0, 0, sqrt(g[1,1]))*1.1), xlim=c(-5*sqrt(g[1,1]), 5*sqrt(g[1,1])), 
     col = "#0fffff40", border = "#0fffff", xlab="", main=paste("T=",t_max,", mu"))
curve(dnorm(x, 0, sqrt(g[1,1])),type="l",lty=2, add=T)
curve(Vectorize(asymp_Mu)(x), type="l", add=T, col="red")

hist(res_alpha0, freq=F, breaks=c(-Inf, seq(-5*sqrt(g[2,2]), 5*sqrt(g[2,2]), 1), Inf), ylim=c(0,dnorm(0, 0, sqrt(g[2,2]))*1.1), xlim=c(-5*sqrt(g[2,2]), 5*sqrt(g[2,2])), 
     col = "#0fffff40", border = "#0fffff", xlab="", main=paste("T=",t_max,", alpha"))
curve(dnorm(x, 0, sqrt(g[2,2])),type="l",lty=2, add=T)
curve(Vectorize(asymp_Alpha)(x), type="l", add=T, col="red")

hist(res_beta0, freq=F, breaks=c(-Inf, seq(-5*sqrt(g[3,3]), 5*sqrt(g[3,3]), 1.2), Inf), ylim=c(0,dnorm(0, 0, sqrt(g[3,3]))*1.1), xlim=c(-5*sqrt(g[3,3]), 5*sqrt(g[3,3])), 
     col = "#0fffff40", border = "#0fffff", xlab="", main=paste("T=",t_max,", beta"))
curve(dnorm(x, 0, sqrt(g[3,3])),type="l",lty=2, add=T)
curve(Vectorize(asymp_Beta)(x), type="l", add=T, col="red")

