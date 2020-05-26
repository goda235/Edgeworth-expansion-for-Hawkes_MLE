M <- length(res_mu0)
n <- 300 #分位点の分割幅

# muについてQ-Q Plot
mini <- max(min(res_alpha0),-5*sqrt(g[1,1]))
maxi <- min(max(res_alpha0),5*sqrt(g[1,1]))
delta <- (maxi-mini)/n
points <- seq(mini, maxi, delta)

F_mu <- c()
i <- 0
for(x in points){
  F_mu <- append(F_mu, integrate(Vectorize(asymp_Mu),-Inf,x)$value)
  print(i/n*100)
  i <- i+1
}

Edgeworth_F_Mu <- c()
for(i in 1:M){
  Edgeworth_F_Mu <- append(Edgeworth_F_Mu, points[sum(F_mu<i/M)+1])
}

Normal_F_Mu <- c()
for(i in 1:M){
  Normal_F_Mu <- append(Normal_F_Mu, pnorm(i/M, 0, sqrt(g[1,1])))
}

# alphaについてQ-Q Plot
mini <- max(min(res_alpha0),-5*sqrt(g[2,2]))
maxi <- min(max(res_alpha0),5*sqrt(g[2,2]))
delta <- (maxi-mini)/n
points <- seq(mini, maxi, delta)

F_alpha <- c()
i <- 0
for(x in points){
  F_alpha<- append(F_alpha, integrate(Vectorize(asymp_Alpha),-Inf,x)$value)
  print(i/n*100)
  i <- i+1
}

Edgeworth_F_Alpha <- c()
for(i in 1:M){
  Edgeworth_F_Alpha <- append(Edgeworth_F_Alpha, points[sum(F_alpha<i/M)+1])
}

Normal_F_Alpha <- c()
for(i in 1:M){
  Normal_F_Alpha <- append(Normal_F_Alpha, pnorm(i/M, 0, sqrt(g[2,2])))
}

# betaについてQ-Q Plot
mini <- max(min(res_beta0),-5*sqrt(g[3,3]))
maxi <- min(max(res_beta0),5*sqrt(g[3,3]))
delta <- (maxi-mini)/n
points <- seq(mini, maxi, delta)

F_beta <- c()
i <- 0
for(x in points){
  F_beta<- append(F_beta, integrate(Vectorize(asymp_Beta),-Inf,x)$value)
  print(i/n*100)
  i <- i+1
}

Edgeworth_F_Beta <- c()
for(i in 1:M){
  Edgeworth_F_Beta <- append(Edgeworth_F_Beta, points[sum(F_beta<i/M)+1])
}

Normal_F_Beta <- c()
for(i in 1:M){
  Normal_F_Beta <- append(Normal_F_Beta, pnorm(i/M, 0, sqrt(g[3,3])))
}

# 結果をPlot

par(oma = c(0, 0, 2.5, 0))
par(mfrow=c(2,3))

qqplot(x=res_mu0[res_mu0 < 5*sqrt(g[1,1])], y=qnorm(seq(0,1,1/M), 0, sqrt(g[1,1])), col="red", cex = 0.8, ylim=c(-4*sqrt(g[1,1]),4*sqrt(g[1,1])), xlab = "mle of mu", ylab = "Normal Distribution", main="Normal Q-Q Plot: mu")
abline(0,1, col="blue",lwd = 2)
grid()
qqplot(x=res_alpha0[res_alpha0 < 5*sqrt(g[2,2])], y=qnorm(seq(0,1,1/M), 0, sqrt(g[2,2])), col="red", cex = 0.8, ylim=c(-4*sqrt(g[2,2]),4*sqrt(g[2,2])), xlab = "mle of alpha", ylab = "Normal Distribution", main="Normal Q-Q Plot: alpha")
abline(0,1, col="blue",lwd = 2)
grid()
qqplot(x=res_beta0[res_beta0 < 5*sqrt(g[3,3])], y=qnorm(seq(0,1,1/M), 0, sqrt(g[3,3])), col="red", cex = 0.8, ylim=c(-4*sqrt(g[3,3]),4*sqrt(g[3,3])), xlab = "mle of beta", ylab = "Normal Distribution", main="Normal Q-Q Plot: beta")
abline(0,1, col="blue",lwd = 2)
grid()

qqplot(x=sort(res_mu0[res_mu0 < 5*sqrt(g[1,1])]), y=Edgeworth_F_Mu,col="red", cex = 0.8, ylim=c(-4*sqrt(g[1,1]),4*sqrt(g[1,1])), xlab = "mle of mu", ylab = "Edgeworth Expantion", main="Edgeworth Q-Q Plot: mu")
abline(0,1, col="blue",lwd = 2)
grid()
qqplot(x=sort(res_alpha0[res_alpha0 < 5*sqrt(g[2,2])]), y=Edgeworth_F_Alpha,col="red", cex = 0.8, ylim=c(-4*sqrt(g[2,2]),4*sqrt(g[2,2])), xlab = "mle of alpha", ylab = "Edgeworth Expantion", main="Edgeworth Q-Q Plot: alpha")
abline(0,1, col="blue",lwd = 2)
grid()
qqplot(x=sort(res_beta0[res_beta0 < 5*sqrt(g[3,3])]), y=Edgeworth_F_Beta,col="red", cex = 0.8, ylim=c(-4*sqrt(g[3,3]),4*sqrt(g[3,3])), xlab = "mle of beta", ylab = "Edgeworth Expantion", main="Edgeworth Q-Q Plot: beta")
abline(0,1, col="blue",lwd = 2)
grid()

mtext(side=3, line=0, outer=T, text = paste("Q-Q Plot, T=",t_max), cex=1.2)

