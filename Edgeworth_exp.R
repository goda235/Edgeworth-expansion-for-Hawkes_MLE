###### Density of Edgeworth expansion ###### 

############################################################################
############################################################################

# We compute the Variation of Z_T by "pracma".
# calcultate a gradient and a hessian of Log-Likelihood l_T(\theta_0).
# params : list(mu, alpha, beta)
# MC : number of MC
# given the inverse of the Fisher Info matrix as g (see another code)
# given the 3rd moments of l_T(\theta_0) as Nu (see another code)

############################################################################
############################################################################
############################################################################
############################################################################

set.seed(0)

library(pracma)
par(mfrow=c(1,3))

params <- c(0.5, 1, 1.3, 0)
t_max <- 300
MC <- 5000

loglik2 <-function(params){
  return(loglik(params,t_max,t))
}

t <- simulate_uni_hawkes(params, t_max)
grad <- grad(loglik2, params[1:3])
hessian <- hessian(loglik2, params[1:3])
Z <- data.frame(t(c(grad, hessian[1:3], hessian[2,1:3], hessian[3,1:3])))

for(i in 1:MC-1){
  t <- simulate_uni_hawkes(params, t_max)
  grad <- grad(loglik2, params[1:3])
  hessian <- hessian(loglik2, params[1:3])
  now_Z <- data.frame(t(c(grad, hessian[1:3], hessian[2,1:3], hessian[3,1:3])))
  Z <- rbind(Z, now_Z)
  print(i/MC*100)
}

#FI <- -hessian/t_max
#(FI)
Sigma <- var(Z)/t_max
# hessianの平均として算出してみる
mean_Sigma <- -matrix(c(mean(Z[,4]),mean(Z[,5]),mean(Z[,6])
                       ,mean(Z[,7]),mean(Z[,8]),mean(Z[,9])
                       ,mean(Z[,10]),mean(Z[,11]),mean(Z[,12])),3,3)/t_max

# Fisherの逆行列
g <- solve(Sigma[1:3,1:3])
g <- solve(mean_Sigma)

for(i in 4:12){
  Z[,i] <- Z[,i] - mean(Z[,i])
}

G11 <- diag(3)*g[1,1]
G12 <- diag(3)*g[1,2]
G13 <- diag(3)*g[1,3]
G21 <- diag(3)*g[2,1]
G22 <- diag(3)*g[2,2]
G23 <- diag(3)*g[2,3]
G31 <- diag(3)*g[3,1]
G32 <- diag(3)*g[3,2]
G33 <- diag(3)*g[3,3]
ZERO <- matrix(0, nrow=3, ncol=3)
C <- rbind(cbind(g,ZERO,ZERO,ZERO),
       cbind(ZERO,G11,G12,G13), 
       cbind(ZERO,G21,G22,G23), 
        cbind(ZERO,G31,G32,G33))
barSigma <- C %*% Sigma %*% t(C)

M <- rbind(cbind(diag(3), matrix(0, nrow=3, ncol=9)),
           cbind((-1)*barSigma[4:12,1:3] %*% Sigma[1:3,1:3], diag(9)))

round(M %*% C %*% Sigma %*% t(C) %*% t(M),1)
#bar_Z <- t(C %*% t(Z))/sqrt(t_max) 
tilde_Z <- t(M %*% C %*% t(Z)) / sqrt(t_max)

V <- array(rep(0, 27), dim = c(3, 3, 3))
tilde_kappa <- array(rep(0, 27), dim = c(3, 3, 3))
Mu <- array(rep(0, 27), dim = c(3, 3, 3))
Nu_modi <- array(rep(0, 27), dim = c(3, 3, 3))

for(i in 1:3){
  for(j in 1:3){
    for(k in 1:3){
      for(l in 1:3){
        V[i,j,k] <- V[i,j,k] + barSigma[3*i+j,l]*Sigma[l,k]
        #V[i,j,k] <- V[i,j,k] + mean(bar_Z[,l]*bar_Z[,3*i+j])*Sigma[l,k]
        Nu_modi[i,j,k] <- Nu_modi[i,j,k] + g[i,l]*Nu[l,j,k]
      }
      tilde_kappa[i,j,k] <- mean(tilde_Z[,i]*tilde_Z[,j]*tilde_Z[,k]) * sqrt(t_max)
    }
  }
}


for(i in 1:3){
  for(j in 1:3){
    for(k in 1:3){
      Mu[i,j,k] <- (V[i,j,k] + V[i,k,j] + Nu_modi[i,j,k])/2
    }
  }
}



hermite_mu_3 <- array(rep(0, 27*2), dim = c(3, 3, 3, 2))

for(i in 1:3){
  for(j in 1:3){
    for(k in 1:3){
      
      a3 <- -Sigma[i,1]*Sigma[j,1]*Sigma[k,1] -
        
       (Sigma[i,1]*Sigma[j,1]*(Sigma[k,2]*g[2,1] + Sigma[k,3]*g[3,1])
         +Sigma[j,1]*Sigma[k,1]*(Sigma[i,2]*g[2,1] + Sigma[i,3]*g[3,1])
         +Sigma[k,1]*Sigma[i,1]*(Sigma[j,2]*g[2,1] + Sigma[j,3]*g[3,1]))/g[1,1] -
        
      (Sigma[i,1]*Sigma[j,2]*Sigma[k,2] 
        + Sigma[j,1]*Sigma[k,2]*Sigma[i,2]
        + Sigma[k,1]*Sigma[i,2]*Sigma[j,2])*g[2,1]^2/g[1,1]^2 -
        
      (Sigma[i,1]*Sigma[j,3]*Sigma[k,3] 
        + Sigma[j,1]*Sigma[k,3]*Sigma[i,3]
        + Sigma[k,1]*Sigma[i,3]*Sigma[j,3])*g[3,1]^2/g[1,1]^2 -
        
      (Sigma[i,1]*Sigma[j,2]*Sigma[k,3] 
        + Sigma[j,1]*Sigma[k,2]*Sigma[i,3]
        + Sigma[k,1]*Sigma[i,2]*Sigma[j,3]
        + Sigma[i,1]*Sigma[k,2]*Sigma[j,3] 
        + Sigma[j,1]*Sigma[i,2]*Sigma[k,3]
        + Sigma[k,1]*Sigma[j,2]*Sigma[i,3])*g[3,1]*g[2,1]/g[1,1]^2 -
        
        (Sigma[i,2]*Sigma[j,2]*Sigma[k,2])*g[2,1]^3/g[1,1]^3 -
        
        (Sigma[i,2]*Sigma[j,2]*Sigma[k,3] +
           Sigma[i,2]*Sigma[j,3]*Sigma[k,2] +
           Sigma[i,3]*Sigma[j,2]*Sigma[k,2])*g[3,1]*g[2,1]^2/g[1,1]^3 -
        
        (Sigma[i,2]*Sigma[j,3]*Sigma[k,3] +
           Sigma[i,3]*Sigma[j,2]*Sigma[k,3] +
           Sigma[i,3]*Sigma[j,3]*Sigma[k,2])*g[3,1]^2*g[2,1]/g[1,1]^3 -
        
        (Sigma[i,3]*Sigma[j,3]*Sigma[k,3])*g[3,1]^3/g[1,1]^3
      
      
      a1 <- -(Sigma[i,1]*Sigma[j,2]*Sigma[k,2] 
              + Sigma[j,1]*Sigma[k,2]*Sigma[i,2]
              + Sigma[k,1]*Sigma[i,2]*Sigma[j,2])*
       (g[2,2] - g[2,1]^2/g[1,1]) -
        
       (Sigma[i,1]*Sigma[j,3]*Sigma[k,3] 
        + Sigma[j,1]*Sigma[k,3]*Sigma[i,3]
        + Sigma[k,1]*Sigma[i,3]*Sigma[j,3])*
       (g[3,3] - g[3,1]^2/g[1,1]) -
        
       (Sigma[i,1]*Sigma[j,2]*Sigma[k,3] 
        + Sigma[j,1]*Sigma[k,2]*Sigma[i,3]
        + Sigma[k,1]*Sigma[i,2]*Sigma[j,3]
        + Sigma[i,1]*Sigma[k,2]*Sigma[j,3] 
        + Sigma[j,1]*Sigma[i,2]*Sigma[k,3]
        + Sigma[k,1]*Sigma[j,2]*Sigma[i,3])*
       (g[2,3] - g[3,1]*g[2,1]/g[1,1]) +
        
       (Sigma[i,j]*(Sigma[k,1]*g[1,1] + Sigma[k,2]*g[2,1] + Sigma[k,3]*g[3,1])
        + Sigma[j,k]*(Sigma[i,1]*g[1,1] + Sigma[i,2]*g[2,1] + Sigma[i,3]*g[3,1])
        + Sigma[k,i]*(Sigma[j,1]*g[1,1] + Sigma[j,2]*g[2,1] + Sigma[j,3]*g[3,1]))/g[1,1] -
        
        (Sigma[i,2]*Sigma[j,2]*Sigma[k,2])*
        3*(g[2,2]-g[2,1]^2/g[1,1])*g[2,1]/g[1,1] -
        
        (Sigma[i,2]*Sigma[j,2]*Sigma[k,3] +
           Sigma[i,2]*Sigma[j,3]*Sigma[k,2] +
           Sigma[i,3]*Sigma[j,2]*Sigma[k,2])*
        ((g[2,2]-g[2,1]^2/g[1,1])*g[3,1]+2*(g[2,3]-g[2,1]*g[3,1]/g[1,1])*g[2,1])/g[1,1] -
        
        (Sigma[i,2]*Sigma[j,3]*Sigma[k,3] +
           Sigma[i,3]*Sigma[j,2]*Sigma[k,3] +
           Sigma[i,3]*Sigma[j,3]*Sigma[k,2])*
        ((g[3,3]-g[3,1]^2/g[1,1])*g[2,1]+2*(g[2,3]-g[2,1]*g[3,1]/g[1,1])*g[3,1])/g[1,1] -
        
        (Sigma[i,3]*Sigma[j,3]*Sigma[k,3])*
        3*(g[3,3]-g[3,1]^2/g[1,1])*g[3,1]/g[1,1]
      
      
      hermite_mu_3[i,j,k,1] <- -a1
      hermite_mu_3[i,j,k,2] <- -a3
    }
  }
}

hermite_mu_1 <- array(rep(0, 3), dim = 3)

for(i in 1:3){
  hermite_mu_1[i] <- (Sigma[i,1]*g[1,1] + Sigma[i,2]*g[2,1] + Sigma[i,3]*g[3,1])/g[1,1] 
}

asymp_Mu <- function(x){
  density <- 0
  
  for(i in 1:3){
    sndMuterm <- 0
    
    for(j in 1:3){
      for(k in 1:3){
        fstMuterm <- 0
        for(l1 in 1:3){
          for(l2 in 1:3){
            fstMuterm <- fstMuterm + Mu[k,l1,l2]*g[l1,i]*g[l2,j]
          }
        }
        
        sndMuterm <- sndMuterm + Mu[i,j,k]*g[j,k]
        trdHermite <- hermite_mu_3[i,j,k,1]*x + hermite_mu_3[i,j,k,2]*x^3
        density <- density + (tilde_kappa[i,j,k]/6 + fstMuterm)*trdHermite
      }
    }
    
    fstHermite <- hermite_mu_1[i]*x
    density <- density + sndMuterm*fstHermite
  }
  
  density <- (density/sqrt(t_max) + 1)*dnorm(x,0,sqrt(g[1,1]))
  density <- max(density, 0)
  
  return(density)
}

curve(Vectorize(asymp_Mu)(x),-10,10,col='red')
curve(dnorm(x,0,sqrt(g[1,1])),-10,10,add=T)


hermite_alpha_3 <- array(rep(0, 27*2), dim = c(3, 3, 3, 2))

for(i in 1:3){
  for(j in 1:3){
    for(k in 1:3){
      
      b3 <- -Sigma[i,2]*Sigma[j,2]*Sigma[k,2] -
        
        (Sigma[i,2]*Sigma[j,2]*(Sigma[k,3]*g[3,2] + Sigma[k,1]*g[1,2])
         +Sigma[j,2]*Sigma[k,2]*(Sigma[i,3]*g[3,2] + Sigma[i,1]*g[1,2])
         +Sigma[k,2]*Sigma[i,2]*(Sigma[j,3]*g[3,2] + Sigma[j,1]*g[1,2]))/g[2,2] -
        
        (Sigma[i,2]*Sigma[j,3]*Sigma[k,3] 
         + Sigma[j,2]*Sigma[k,3]*Sigma[i,3]
         + Sigma[k,2]*Sigma[i,3]*Sigma[j,3])*g[3,2]^2/g[2,2]^2 -
        
        (Sigma[i,2]*Sigma[j,1]*Sigma[k,1] 
         + Sigma[j,2]*Sigma[k,1]*Sigma[i,1]
         + Sigma[k,2]*Sigma[i,1]*Sigma[j,1])*g[1,2]^2/g[2,2]^2 -
        
        (Sigma[i,2]*Sigma[j,3]*Sigma[k,1] 
         + Sigma[j,2]*Sigma[k,3]*Sigma[i,1]
         + Sigma[k,2]*Sigma[i,3]*Sigma[j,1]
         + Sigma[i,2]*Sigma[k,3]*Sigma[j,1] 
         + Sigma[j,2]*Sigma[i,3]*Sigma[k,1]
         + Sigma[k,2]*Sigma[j,3]*Sigma[i,1])*g[1,2]*g[3,2]/g[2,2]^2 -
        
        (Sigma[i,3]*Sigma[j,3]*Sigma[k,3])*g[3,2]^3/g[2,2]^3 -
        
        (Sigma[i,3]*Sigma[j,3]*Sigma[k,1] +
           Sigma[i,3]*Sigma[j,1]*Sigma[k,3] +
           Sigma[i,1]*Sigma[j,3]*Sigma[k,3])*g[1,2]*g[3,2]^2/g[2,2]^3 -
        
        (Sigma[i,3]*Sigma[j,1]*Sigma[k,1] +
           Sigma[i,1]*Sigma[j,3]*Sigma[k,1] +
           Sigma[i,1]*Sigma[j,1]*Sigma[k,3])*g[1,2]^2*g[3,2]/g[2,2]^3 -
        
        (Sigma[i,1]*Sigma[j,1]*Sigma[k,1])*g[1,2]^3/g[2,2]^3
      
      
      b1 <- -(Sigma[i,2]*Sigma[j,3]*Sigma[k,3] 
              + Sigma[j,2]*Sigma[k,3]*Sigma[i,3]
              + Sigma[k,2]*Sigma[i,3]*Sigma[j,3])*
        (g[3,3] - g[3,2]^2/g[2,2]) -
        
        (Sigma[i,2]*Sigma[j,1]*Sigma[k,1] 
         + Sigma[j,2]*Sigma[k,1]*Sigma[i,1]
         + Sigma[k,2]*Sigma[i,1]*Sigma[j,1])*
        (g[1,1] - g[1,2]^2/g[2,2]) -
        
        (Sigma[i,2]*Sigma[j,3]*Sigma[k,1] 
         + Sigma[j,2]*Sigma[k,3]*Sigma[i,1]
         + Sigma[k,2]*Sigma[i,3]*Sigma[j,1]
         + Sigma[i,2]*Sigma[k,3]*Sigma[j,1] 
         + Sigma[j,2]*Sigma[i,3]*Sigma[k,1]
         + Sigma[k,2]*Sigma[j,3]*Sigma[i,1])*
        (g[3,1] - g[1,2]*g[3,2]/g[2,2]) +
        
        (Sigma[i,j]*(Sigma[k,2]*g[2,2] + Sigma[k,3]*g[3,2] + Sigma[k,1]*g[1,2])
         + Sigma[j,k]*(Sigma[i,2]*g[2,2] + Sigma[i,3]*g[3,2] + Sigma[i,1]*g[1,2])
         + Sigma[k,i]*(Sigma[j,2]*g[2,2] + Sigma[j,3]*g[3,2] + Sigma[j,1]*g[1,2]))/g[2,2] -
        
        (Sigma[i,3]*Sigma[j,3]*Sigma[k,3])*
        3*(g[3,3]-g[3,2]^2/g[2,2])*g[3,2]/g[2,2] -
        
        (Sigma[i,3]*Sigma[j,3]*Sigma[k,1] +
           Sigma[i,3]*Sigma[j,1]*Sigma[k,3] +
           Sigma[i,1]*Sigma[j,3]*Sigma[k,3])*
        ((g[3,3]-g[3,2]^2/g[2,2])*g[1,2]+2*(g[3,1]-g[3,2]*g[1,2]/g[2,2])*g[3,2])/g[2,2] -
        
        (Sigma[i,3]*Sigma[j,1]*Sigma[k,1] +
           Sigma[i,1]*Sigma[j,3]*Sigma[k,1] +
           Sigma[i,1]*Sigma[j,1]*Sigma[k,3])*
        ((g[1,1]-g[1,2]^2/g[2,2])*g[3,2]+2*(g[3,1]-g[3,2]*g[1,2]/g[2,2])*g[1,2])/g[2,2] -
        
        (Sigma[i,1]*Sigma[j,1]*Sigma[k,1])*
        3*(g[1,1]-g[1,2]^2/g[2,2])*g[1,2]/g[2,2]
      
      
      hermite_alpha_3[i,j,k,1] <- -b1
      hermite_alpha_3[i,j,k,2] <- -b3
    }
  }
}

hermite_alpha_1 <- array(rep(0, 3), dim = 3)

for(i in 1:3){
  hermite_alpha_1[i] <- (Sigma[i,2]*g[2,2] + Sigma[i,3]*g[3,2] + Sigma[i,1]*g[1,2])/g[2,2] 
}

asymp_Alpha <- function(x){
  density <- 0
  
  for(i in 1:3){
    sndMuterm <- 0
    
    for(j in 1:3){
      for(k in 1:3){
        fstMuterm <- 0
        for(l1 in 1:3){
          for(l2 in 1:3){
            fstMuterm <- fstMuterm + Mu[k,l1,l2]*g[l1,i]*g[l2,j]
          }
        }
        
        sndMuterm <- sndMuterm + Mu[i,j,k]*g[j,k]
        trdHermite <- hermite_alpha_3[i,j,k,1]*x + hermite_alpha_3[i,j,k,2]*x^3
        density <- density + (tilde_kappa[i,j,k]/6 + fstMuterm)*trdHermite
      }
    }
    
    fstHermite <- hermite_alpha_1[i]*x
    density <- density + sndMuterm*fstHermite
  }
  
  density <- (density/sqrt(t_max) + 1)*dnorm(x,0,sqrt(g[2,2]))
  density <- max(density, 0)
  
  return(density)
}

curve(Vectorize(asymp_Alpha)(x),-10,10,col='red')
curve(dnorm(x,0,sqrt(g[2,2])),-10,10,add=T)


hermite_beta_3 <- array(rep(0, 27*2), dim = c(3, 3, 3, 2))

for(i in 1:3){
  for(j in 1:3){
    for(k in 1:3){
      
      c3 <- -Sigma[i,3]*Sigma[j,3]*Sigma[k,3] -
        
        (Sigma[i,3]*Sigma[j,3]*(Sigma[k,1]*g[1,3] + Sigma[k,2]*g[2,3])
         +Sigma[j,3]*Sigma[k,3]*(Sigma[i,1]*g[1,3] + Sigma[i,2]*g[2,3])
         +Sigma[k,3]*Sigma[i,3]*(Sigma[j,1]*g[1,3] + Sigma[j,2]*g[2,3]))/g[3,3] -
        
        (Sigma[i,3]*Sigma[j,1]*Sigma[k,1] 
         + Sigma[j,3]*Sigma[k,1]*Sigma[i,1]
         + Sigma[k,3]*Sigma[i,1]*Sigma[j,1])*g[1,3]^2/g[3,3]^2 -
        
        (Sigma[i,3]*Sigma[j,2]*Sigma[k,2] 
         + Sigma[j,3]*Sigma[k,2]*Sigma[i,2]
         + Sigma[k,3]*Sigma[i,2]*Sigma[j,2])*g[2,3]^2/g[3,3]^2 -
        
        (Sigma[i,2]*Sigma[j,3]*Sigma[k,1] 
         + Sigma[j,2]*Sigma[k,3]*Sigma[i,1]
         + Sigma[k,2]*Sigma[i,3]*Sigma[j,1]
         + Sigma[i,2]*Sigma[k,3]*Sigma[j,1] 
         + Sigma[j,2]*Sigma[i,3]*Sigma[k,1]
         + Sigma[k,2]*Sigma[j,3]*Sigma[i,1])*g[2,3]*g[1,3]/g[3,3]^2 -
        
        (Sigma[i,1]*Sigma[j,1]*Sigma[k,1])*g[1,3]^3/g[3,3]^3 -
        
        (Sigma[i,1]*Sigma[j,1]*Sigma[k,2] +
           Sigma[i,1]*Sigma[j,2]*Sigma[k,1] +
           Sigma[i,2]*Sigma[j,1]*Sigma[k,1])*g[2,3]*g[1,3]^2/g[3,3]^3 -
        
        (Sigma[i,1]*Sigma[j,2]*Sigma[k,2] +
           Sigma[i,2]*Sigma[j,1]*Sigma[k,2] +
           Sigma[i,2]*Sigma[j,2]*Sigma[k,1])*g[2,3]^2*g[1,3]/g[3,3]^3 -
        
        (Sigma[i,2]*Sigma[j,2]*Sigma[k,2])*g[2,3]^3/g[3,3]^3
      
      
      c1 <- -(Sigma[i,3]*Sigma[j,1]*Sigma[k,1] 
              + Sigma[j,3]*Sigma[k,1]*Sigma[i,1]
              + Sigma[k,3]*Sigma[i,1]*Sigma[j,1])*
        (g[1,1] - g[1,3]^2/g[3,3]) -
        
        (Sigma[i,3]*Sigma[j,2]*Sigma[k,2] 
         + Sigma[j,3]*Sigma[k,2]*Sigma[i,2]
         + Sigma[k,3]*Sigma[i,2]*Sigma[j,2])*
        (g[2,2] - g[2,3]^2/g[3,3]) -
        
        (Sigma[i,2]*Sigma[j,3]*Sigma[k,1] 
         + Sigma[j,2]*Sigma[k,3]*Sigma[i,1]
         + Sigma[k,2]*Sigma[i,3]*Sigma[j,1]
         + Sigma[i,2]*Sigma[k,3]*Sigma[j,1] 
         + Sigma[j,2]*Sigma[i,3]*Sigma[k,1]
         + Sigma[k,2]*Sigma[j,3]*Sigma[i,1])*
        (g[1,2] - g[2,3]*g[1,3]/g[3,3]) +
        
        (Sigma[i,j]*(Sigma[k,3]*g[3,3] + Sigma[k,1]*g[1,3] + Sigma[k,2]*g[2,3])
         + Sigma[j,k]*(Sigma[i,3]*g[3,3] + Sigma[i,1]*g[1,3] + Sigma[i,2]*g[2,3])
         + Sigma[k,i]*(Sigma[j,3]*g[3,3] + Sigma[j,1]*g[1,3] + Sigma[j,2]*g[2,3]))/g[3,3] -
        
        (Sigma[i,1]*Sigma[j,1]*Sigma[k,1])*
        3*(g[1,1]-g[1,3]^2/g[3,3])*g[1,3]/g[3,3] -
        
        (Sigma[i,1]*Sigma[j,1]*Sigma[k,2] +
           Sigma[i,1]*Sigma[j,2]*Sigma[k,1] +
           Sigma[i,2]*Sigma[j,1]*Sigma[k,1])*
        ((g[1,1]-g[1,3]^2/g[3,3])*g[2,3]+2*(g[1,2]-g[1,3]*g[2,3]/g[3,3])*g[1,3])/g[3,3] -
        
        (Sigma[i,1]*Sigma[j,2]*Sigma[k,2] +
           Sigma[i,2]*Sigma[j,1]*Sigma[k,2] +
           Sigma[i,2]*Sigma[j,2]*Sigma[k,1])*
        ((g[2,2]-g[2,3]^2/g[3,3])*g[1,3]+2*(g[1,2]-g[1,3]*g[2,3]/g[3,3])*g[2,3])/g[3,3] -
        
        (Sigma[i,2]*Sigma[j,2]*Sigma[k,2])*
        3*(g[2,2]-g[2,3]^2/g[3,3])*g[2,3]/g[3,3]
      
      
      hermite_beta_3[i,j,k,1] <- -c1
      hermite_beta_3[i,j,k,2] <- -c3
    }
  }
}

hermite_beta_1 <- array(rep(0, 3), dim = 3)

for(i in 1:3){
  hermite_beta_1[i] <- (Sigma[i,3]*g[3,3] + Sigma[i,1]*g[1,3] + Sigma[i,2]*g[2,3])/g[3,3] 
}

asymp_Beta <- function(x){
  density <- 0
  
  for(i in 1:3){
    sndMuterm <- 0
    
    for(j in 1:3){
      for(k in 1:3){
        fstMuterm <- 0
        for(l1 in 1:3){
          for(l2 in 1:3){
            fstMuterm <- fstMuterm + Mu[k,l1,l2]*g[l1,i]*g[l2,j]
          }
        }
        
        sndMuterm <- sndMuterm + Mu[i,j,k]*g[j,k]
        trdHermite <- hermite_beta_3[i,j,k,1]*x + hermite_beta_3[i,j,k,2]*x^3
        density <- density + (tilde_kappa[i,j,k]/6 + fstMuterm)*trdHermite
      }
    }
    
    fstHermite <- hermite_beta_1[i]*x
    density <- density + sndMuterm*fstHermite
  }
  
  density <- (density/sqrt(t_max) + 1)*dnorm(x,0,sqrt(g[3,3]))
  density <- max(density, 0)
  
  return(density)
}

curve(Vectorize(asymp_Beta)(x),-10,10,col='red')
curve(dnorm(x,0,sqrt(g[3,3])),-10,10,add=T)


