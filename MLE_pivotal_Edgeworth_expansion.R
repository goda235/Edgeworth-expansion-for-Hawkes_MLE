library(pracma)
library(tcltk)
########################################################################################################
########################################################################################################
########################################################################################################

# Functions

## Simulator of Hawkes process (Ogata's method) 
simulate_uni_hawkes <- function(params, t_max){
  mu     <- params[1]
  alpha  <- params[2]
  beta   <- params[3]
  arrivals <- 0 #t=0
  s <- -log(runif(1))/mu #interarrival time
  t <- s #jump time
  lambda <- mu + alpha #new value of lambda
  arrivals<- c(arrivals,t)
  while(t < t_max) {
    D <- 1+(beta*log(runif(1)))/(lambda - mu)
    if(D>=0){
      s <- min(-log(D)/beta,-log(runif(1))/mu)
    }else{
      s<- -log(runif(1))/mu
    }
    lambda <- mu + (lambda-mu)*exp(-beta*s) + alpha
    t <- t + s
    if(t < t_max){
      arrivals <- c(arrivals,t)}
  }
  return(c(arrivals,t_max))
}

## log-likelihood function
loglik <- function(params, arrivals, t_max){
  mu <- params[1]
  alpha <- params[2]
  beta <- params[3]
  x <- 0
  omega <- arrivals[2:(length(arrivals)-1)]
  Ai <- exp(-beta*(omega))
  f <- function(z) {
    sum(exp(-beta*(omega[z] - omega[1:(z-1)])))
  }
  Bi <- c(0, sapply(2:length(omega), f))
  term_1 <- sum(log(mu + x*Ai + alpha*Bi))
  term_2 <- mu*t_max - x/beta*(exp(-beta*t_max) - 1)
  term_3 <- -sum(alpha/beta*(exp(-beta*(t_max - omega)) - 1))
  
  return(term_1 - (term_2 + term_3))
}

## Generator of MC path of Z_T
generate_Z <- function(params, t_max, MC){
  t <- simulate_uni_hawkes(params, t_max)
  loglik2 <- function(params) loglik(params, t, t_max) # log-likelihood function given (arrivals, t_max)
  grad <- grad(loglik2, params)
  hessian <- hessian(loglik2, params)
  Z <- data.frame(t(c(grad, hessian[1,1:3], hessian[2,1:3], hessian[3,1:3])))
  pb <- txtProgressBar(min = 1, max = MC, style = 3) # Progress bar
  print("MonteCalro simulate for Z_T")
  while(length(Z$X1)<MC){
    setTxtProgressBar(pb, length(Z$X1))
    t <- simulate_uni_hawkes(params, t_max)
    logik2 <- function(params) loglik(params, t, t_max)
    grad <- grad(loglik2, params)
    hessian <- hessian(loglik2, params)
    new_Z <- data.frame(t(c(grad, hessian[1,1:3], hessian[2,1:3], hessian[3,1:3])))
    Z <- rbind(Z, new_Z)
    Z <- na.omit(Z)
  }
  for(i in 4:12){
    Z[,i] <- Z[,i] - mean(Z[,i]) # Centering
  }
  return(Z)
}

## Coefficients
Coefficients <- function(params, t_max, MC){
  Z <- generate_Z(params, t_max, MC)
  ### Fisher infomation
  Sigma <- var(Z)/t_max
  FI <- Sigma[1:3,1:3]
  ### Square root of Fisher infomation
  U <- svd(FI)$u
  V <- svd(FI)$v
  D <- diag(sqrt(svd(FI)$d))
  sqrt_FI <- U %*% D %*% t(V)
  inv_sqrt_FI <- solve(sqrt_FI)
  # Matrix C
  G11 <- diag(3)*inv_sqrt_FI[1,1]
  G12 <- diag(3)*inv_sqrt_FI[1,2]
  G13 <- diag(3)*inv_sqrt_FI[1,3]
  G21 <- diag(3)*inv_sqrt_FI[2,1]
  G22 <- diag(3)*inv_sqrt_FI[2,2]
  G23 <- diag(3)*inv_sqrt_FI[2,3]
  G31 <- diag(3)*inv_sqrt_FI[3,1]
  G32 <- diag(3)*inv_sqrt_FI[3,2]
  G33 <- diag(3)*inv_sqrt_FI[3,3]
  ZERO <- matrix(0, nrow=3, ncol=3)
  C <- rbind(cbind(inv_sqrt_FI,ZERO,ZERO,ZERO),
             cbind(ZERO,G11,G12,G13), 
             cbind(ZERO,G21,G22,G23), 
             cbind(ZERO,G31,G32,G33))
  modi_Sigma <- C %*% Sigma %*% t(C)
  modi_Z <- t(C %*% t(Z)) / sqrt(t_max) # modified Z_T
  # Coefficients
  Mu <- array(rep(0, 27), dim = c(3, 3, 3))
  kappa <- array(rep(0, 27), dim = c(3, 3, 3))
  for(i in 1:3){
    for(j in 1:3){
      for(k in 1:3){
        for(l in 1:3){
          Mu[i,j,k] <- Mu[i,j,k] + 0.5*inv_sqrt_FI[l,j]*modi_Sigma[3*i+l,k]
        }
        kappa[i,j,k] <- mean(modi_Z[,i]*modi_Z[,j]*modi_Z[,k]) * sqrt(t_max)
      }
    }
  }
  res <- list(Mu=Mu, kappa=kappa)
  return(res)
}

## Marginal distribution
### Coef are given
### 1st-dim
margdist1 <- function(x){
  FirstHermitCoef <- 0
  for(i in 1:3){
    FirstHermitCoef <- FirstHermitCoef + Coef$Mu[1,i,i]
  }
  density <- (1+((Coef$kappa[1,1,1]/6 + Coef$Mu[1,1,1])*(x^3-3*x) + FirstHermitCoef*x)/sqrt(t_max))*dnorm(x)
  density <- max(density, 0)
  return(density)
}

### 2nd-dim
margdist2 <- function(x){
  FirstHermitCoef <- 0
  for(i in 1:3){
    FirstHermitCoef <- FirstHermitCoef + Coef$Mu[2,i,i]
  }
  density <- (1+((Coef$kappa[2,2,2]/6 + Coef$Mu[2,2,2])*(x^3-3*x) + FirstHermitCoef*x)/sqrt(t_max))*dnorm(x)
  density <- max(density, 0)
  return(density)
}

### 3rd-dim
margdist3 <- function(x){
  FirstHermitCoef <- 0
  for(i in 1:3){
    FirstHermitCoef <- FirstHermitCoef + Coef$Mu[3,i,i]
  }
  density <- (1+((Coef$kappa[3,3,3]/6 + Coef$Mu[3,3,3])*(x^3-3*x) + FirstHermitCoef*x)/sqrt(t_max))*dnorm(x)
  density <- max(density, 0)
  return(density)
}

# MonteCalro simulate for MLE
MC_MLE <- function(params, MC){
  res <- list("First"=rep(0,MC), "Second"=rep(0,MC), "Third"=rep(0,MC))
  pb <- txtProgressBar(min = 1, max = MC, style = 3) # Progress bar
  print("MonteCalro simulate for MLE")
  for(i in 1:MC){
    setTxtProgressBar(pb, i)
    t <- simulate_uni_hawkes(params, t_max)
    mle <- optim(params, loglik, arrivals = t, t_max=t_max, method = "L-BFGS-B", control=list(fnscale=-1),
                    lower=c(0.001,  0.001, 0.001), upper=c(100, 100, 100))$par
    ### Observed Fisher infomation
    loglik2 <- function(params) loglik(params, t, t_max) # log-likelihood function given (arrivals, t_max)
    FI <- (-1)*hessian(loglik2, mle)/t_max
    ### Square root of the observed Fisher infomation
    U <- svd(FI)$u
    V <- svd(FI)$v
    D <- diag(sqrt(svd(FI)$d))
    sqrt_FI <- U %*% D %*% t(V)
      
    modi_mle <- sqrt_FI%*%(mle-params)*sqrt(t_max)
    res$First[i] <- modi_mle[1]
    res$Second[i] <- modi_mle[2]
    res$Third[i] <- modi_mle[3]
  }
  return(res)
}

## Quantiles for Edgeworth expansion
### Coef are given
### N_Q : Number of divisions for Quantiles
### N_V : Number of divisions for Values
### maxi : width of grid 
### mini : width of grid 
Edge_Quant <- function(N_Q, N_V=1000, upper=10, lower=-10){
  points <- seq(lower, upper, (upper-lower)/N_V)
  
  CDF <- list("First"=c(), "Second"=c(), "Third"=c())
  for(x in points){
    CDF$First <- append(CDF$First, integrate(Vectorize(margdist1), -Inf, x)$value)
    CDF$Second <- append(CDF$Second, integrate(Vectorize(margdist2), -Inf, x)$value)
    CDF$Third <- append(CDF$Third, integrate(Vectorize(margdist3), -Inf, x)$value)
  }
  
  Quant <- list("First"=c(), "Second"=c(), "Third"=c())
  for(i in 1:N_Q){
    Quant$First <- append(Quant$First, points[sum(CDF$First<i/(N_Q+1))])
    Quant$Second <- append(Quant$Second, points[sum(CDF$Second<i/(N_Q+1))])
    Quant$Third <- append(Quant$Third, points[sum(CDF$Third<i/(N_Q+1))])
  }
  return(Quant)
}

# Bootstrap distribution density function
Boot_pdf <- function(x, data, upper=10, lower=-10, break_num=30){
  len <- (upper-lower)/break_num
  breaks <- c(-Inf, seq(lower, upper, length=break_num), Inf)
  for(n in 1:length(breaks)){
    if(x>=breaks[n] & x<=breaks[n+1]){
      up <- breaks[n+1]
      down <- breaks[n]
      break
    }else if(x <= lower){
      up <- lower
      down <- -Inf
      break
    }else if(x >= upper){
      up <- Inf
      down <- upper
      break
    }
  }
  return(sum((data <= up & data > down))/length(data)/len)
}

########################################################################################################
########################################################################################################
########################################################################################################

## Setting
set.seed(0)
params <- c(0.3, 1., 1.4) #(mu, alpha, beta)
t_max <- 100
MC_coef <- 1000
MC_boot <- 1000
MC_hist <- 1000

## MLE of a Hawkes process
t <- simulate_uni_hawkes(params, t_max)
mle  <- optim(params, loglik, arrivals = t, t_max = t_max, method = "L-BFGS-B", control=list(fnscale=-1),
             lower=c(0.001,  0.001, 0.001), upper=c(100, 100, 100))$par
print(mle)
## Compute coefficients of Edgeworth expansion
Coef <- Coefficients(mle, t_max, MC_coef)
edge_quant <- Edge_Quant(MC_coef)

## Compute Bootstrap distribution
boot <- MC_MLE(mle, MC_boot)

## Compute MLEs for histgrams
mc_mle <- MC_MLE(params, MC_hist)

## Histgrams: Standard Normal VS Edgeworth Exp
par(mfrow=c(1,3))
hist(mc_mle$First, freq=F, breaks=c(-Inf, seq(-5, 5, length=30), Inf), ylim=c(0, 0.5), xlim=c(-5, 5), 
     col = "gray", border = "darkgray", xlab="", main=paste("T=",t_max,", 1st marginal"))
curve(dnorm(x),type="l",lty=2, add=T)
curve(Vectorize(margdist1)(x), type="l", add=T, col="red")

hist(mc_mle$Second, freq=F, breaks=c(-Inf, seq(-5, 5, length=30), Inf), ylim=c(0, 0.5), xlim=c(-5, 5),  
     col = "gray", border = "darkgray", xlab="", main=paste("T=",t_max,", 2nd marginal"))
curve(dnorm(x),type="l",lty=2, add=T)
curve(Vectorize(margdist2)(x), type="l", add=T, col="red")

hist(mc_mle$Third, freq=F, breaks=c(-Inf, seq(-5, 5, length=30), Inf), ylim=c(0, 0.5), xlim=c(-5, 5), 
     col = "gray", border = "darkgray", xlab="", main=paste("T=",t_max,", 3rd marginal"))
curve(dnorm(x),type="l",lty=2, add=T)
curve(Vectorize(margdist3)(x), type="l", add=T, col="red")

## Histgrams: Bootstrap VS Edgeworth Exp
boot_pdf <- function(x) Boot_pdf(x, data=boot$First, break_num=60, upper=10, lower=-10)
hist(mc_mle$First, freq=F, breaks=c(-Inf, seq(-5, 5, length=30), Inf), ylim=c(0, 0.5), xlim=c(-5, 5), 
     col = "gray", border = "darkgray", xlab="", main=paste("T=",t_max,", 1st marginal"))
curve(Vectorize(boot_pdf)(x), type="s", col="black", add=T)
curve(Vectorize(margdist1)(x), type="l", add=T, col="red")

boot_pdf<- function(x) Boot_pdf(x, data=boot$Second, break_num=60, upper=10, lower=-10)
hist(mc_mle$Second, freq=F, breaks=c(-Inf, seq(-5, 5, length=30), Inf), ylim=c(0, 0.5), xlim=c(-5, 5),  
     col = "gray", border = "darkgray", xlab="", main=paste("T=",t_max,", 2nd marginal"))
curve(Vectorize(boot_pdf)(x), type="s", col="black", add=T)
curve(Vectorize(margdist2)(x), type="l", add=T, col="red")

boot_pdf <- function(x) Boot_pdf(x, data=boot$Third, break_num=60, upper=10, lower=-10)
hist(mc_mle$Third, freq=F, breaks=c(-Inf, seq(-5, 5, length=30), Inf), ylim=c(0, 0.5), xlim=c(-5, 5), 
     col = "gray", border = "darkgray", xlab="", main=paste("T=",t_max,", 3rd marginal"))
curve(Vectorize(boot_pdf)(x), type="s", col="black", add=T)
curve(Vectorize(margdist3)(x), type="l", add=T, col="red")

# Q-Q Plot
par(oma = c(0, 0, 2.5, 0))
par(mfrow=c(3,3))

norm_quant <- qnorm(seq(0,1,1/(MC_coef+1)))[2:(MC_coef+1)]
qqplot(x=mc_mle$First, y=norm_quant, col="red", cex = 0.8, xlim=c(-4,3), ylim=c(-5,4), xlab = "MLE", ylab = "Normal Distribution", main="Normal Q-Q Plot: 1st marginal")
abline(0,1, col="blue",lwd = 2)
grid()
qqplot(x=mc_mle$Second, y=norm_quant, col="red", cex = 0.8, xlim=c(-4,3), ylim=c(-5,4), xlab = "MLE", ylab = "Normal Distribution", main="Normal Q-Q Plot: 2nd marginal")
abline(0,1, col="blue",lwd = 2)
grid()
qqplot(x=mc_mle$Third, y=norm_quant, col="red", cex = 0.8, xlim=c(-4,3), ylim=c(-5,4), xlab = "MLE", ylab = "Normal Distribution", main="Normal Q-Q Plot: 3rd marginal")
abline(0,1, col="blue",lwd = 2)
grid()

qqplot(x=mc_mle$First, y=edge_quant$First, col="red", cex = 0.8, xlim=c(-4,3), ylim=c(-5,4), xlab = "MLE", ylab = "Edgeworth Expantion", main="Edgeworth Q-Q Plot: 1st marginal")
abline(0,1, col="blue",lwd = 2)
grid()
qqplot(x=mc_mle$Second, y=edge_quant$Second, col="red", cex = 0.8, xlim=c(-4,3), ylim=c(-5,4), xlab = "MLE", ylab = "Edgeworth Expantion", main="Edgeworth Q-Q Plot: 2nd marginal")
abline(0,1, col="blue",lwd = 2)
grid()
qqplot(x=mc_mle$Third, y=edge_quant$Third, col="red", cex = 0.8, xlim=c(-4,3), ylim=c(-5,4), xlab = "MLE", ylab = "Edgeworth Expantion", main="Edgeworth Q-Q Plot: 3rd marginal")
abline(0,1, col="blue",lwd = 2)
grid()

qqplot(x=mc_mle$First, y=boot$First, col="red", cex = 0.8, xlim=c(-4,3), ylim=c(-5,4), xlab = "MLE", ylab = "Bootstrap", main="Bootstrap Q-Q Plot: 1st marginal")
abline(0,1, col="blue",lwd = 2)
grid()
qqplot(x=mc_mle$Second, y=boot$Second, col="red", cex = 0.8, xlim=c(-4,3), ylim=c(-5,4), xlab = "MLE", ylab = "Bootstrap", main="Bootstrap Q-Q Plot: 2nd marginal")
abline(0,1, col="blue",lwd = 2)
grid()
qqplot(x=mc_mle$Third, y=boot$Third, col="red", cex = 0.8, xlim=c(-4,3), ylim=c(-5,4), xlab = "MLE", ylab = "Bootstrap", main="Bootstrap Q-Q Plot: 3rd marginal")
abline(0,1, col="blue",lwd = 2)
grid()

mtext(side=3, line=0, outer=T, text = paste("Q-Q Plot, T=",t_max), cex=1.2)

