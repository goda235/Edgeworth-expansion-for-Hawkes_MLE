###### Moment (Fisher Infomation) ###### 

############################################################################
############################################################################

# Compute the expectation for derivatives of l_T(\theta_0)/T, inparticular the Fisher Infomation matrix by the Monte Carlo method.

set.seed(0)

params <-  c(0.5, 1, 1.3, 0) # (mu, alpha, beta, x=0)
t_max <- 300
MC <- 5000 # num of MC

########################################################################################################################
########################################################################################################################
########################################################################################################################

# Hawkes core process X
kernel_X1 <- function(params, s, t){
  alpha  <- params[2]
  beta   <- params[3]
  return(alpha*exp(-beta*(t-s)))
}
kernel_X2 <- function(params, s, t){
  alpha  <- params[2]
  beta   <- params[3]
  return(alpha*(t-s)*exp(-beta*(t-s)))
}
kernel_X3 <- function(params, s, t){
  alpha  <- params[2]
  beta   <- params[3]
  return(alpha*(t-s)^2*exp(-beta*(t-s)))
}

X <- function(params, arrivals){
  total_jump <- length(arrivals) #最後はTでの値．これも含める
  X1 <-rep(0, total_jump)
  X2 <-rep(0, total_jump)
  X3 <-rep(0, total_jump)
  for(i in 2:(total_jump-1)){
    X1 <- X1 + kernel_X1(params, arrivals[i], append(rep(Inf,i), arrivals[(i+1):total_jump]))
    X2 <- X2 + kernel_X2(params, arrivals[i], append(rep(arrivals[i],i), arrivals[(i+1):total_jump]))
    X3 <- X3 + kernel_X3(params, arrivals[i], append(rep(arrivals[i],i), arrivals[(i+1):total_jump]))
  }
  return(list(X1,X2,X3))
}

########################################################################################################################
########################################################################################################################
########################################################################################################################

Fisher_MC <- function(params, t_max, MC){
  mu     <- params[1]
  alpha  <- params[2]
  beta   <- params[3]
  x <- params[4]
  
  g11 <- 0
  g12 <- 0
  g13 <- 0
  
  g22 <- 0
  g23 <- 0
  
  g33 <- 0
  
  for(i in 1:MC){
    t <- simulate_uni_hawkes(params, t_max)
    X <-  X(params, t)
    
    X1 <- X[[1]]
    X2 <- X[[2]]
    
    g11 <- g11 + sum(1/(mu+X1)^2)
    g12 <- g12 + sum(X1/(alpha*(mu+X1)^2))
    g13 <- g13 + sum(-X2/(mu+X1)^2)
    
    g22 <- g22 + sum(X1^2/(alpha*(mu+X1))^2)
    g23 <- g23 + sum(-X1*X2/(alpha*(mu+X1)^2))
    
    g33 <- g33 + sum(X2^2/(mu+X1)^2)
  }
  
  FI <- matrix(c(g11,g12,g13, g12,g22,g23, g13,g23,g33), 3, 3)
  
  return(FI/(MC*t_max))
}

Third_Tensor_MC <- function(params, t_max, MC){
  mu     <- params[1]
  alpha  <- params[2]
  beta   <- params[3]
  x <- params[4]
  
  g111 <- 0
  g112 <- 0
  g113 <- 0
  
  g221 <- 0
  g222 <- 0
  g223 <- 0
  
  g331 <- 0
  g332 <- 0
  g333 <- 0
  
  g123 <- 0
  
  for(i in 1:MC){
    t <- simulate_uni_hawkes(params, t_max)
    X <-  X(params, t)
    
    X1 <- X[[1]]
    X2 <- X[[2]]
    X3 <- X[[3]]
    
    g111 <- g111 + sum(2/(mu+X1)^3)
    g112 <- g112 + sum(2*X1/(alpha*(mu+X1)^3))
    g113 <- g113 + sum(-2*X2/(mu+X1)^3)
    
    g221 <- g221 + sum(2*X1^2/(alpha^2*(mu+X1)^3))
    g222 <- g222 + sum(2*X1^3/(alpha*(mu+X1))^3)
    g223 <- g223 + sum(2*mu*X1*X2/(alpha^2*(mu+X1)^3))
    
    g331 <- g331 + sum((2*X2^2-X3*(mu+X1))/(mu+X1)^3)
    g332 <- g332 + sum((-X1*X3*(mu+X1)-2*mu*X2^2)/(alpha*(mu+X1)^3))
    g333 <- g333 + sum((3*X2*X3*(mu+X1)-2*X2^3)/(mu+X1)^3)
    
    g123 <- g123 + sum((mu*X2-X1*X2)/(alpha*(mu+X1)^3))
    
  }
  
 Nu <- array(rep(0, 27), dim = c(3, 3, 3))
 
 Nu[1,1,1] = g111 / (MC*t_max)
 Nu[1,1,2] = Nu[1,2,1] = Nu[2,1,1] = g112 / (MC*t_max)
 Nu[1,1,3] = Nu[1,3,1] = Nu[3,1,1] = g113 / (MC*t_max)
 Nu[2,2,1] = Nu[2,1,2] = Nu[1,2,2] = g221 / (MC*t_max)
 Nu[2,2,2] = g222 / (MC*t_max)
 Nu[2,2,3] = Nu[2,3,2] = Nu[3,2,2] = g223 / (MC*t_max)
 Nu[3,3,1] = Nu[3,1,3] = Nu[1,3,3] = g331 / (MC*t_max)
 Nu[3,3,2] = Nu[3,2,3] = Nu[2,3,3] = g332 / (MC*t_max)
 Nu[3,3,3] = g333 / (MC*t_max)
 Nu[1,2,3] = Nu[1,3,2] = Nu[2,3,1] = Nu[2,1,3] = Nu[3,1,2] = Nu[3,2,1] = g123 / (MC*t_max)
 
 return(Nu)
}

########################################################################################################################
########################################################################################################################
########################################################################################################################

# Fisher Infomation
(FI <- Fisher_MC(params, t_max, MC))
# inverse
(g <- solve(FI)) 

# 3rd derivatives
(Nu <- Third_Tensor_MC(params, t_max, MC))
