##annealing
library(GenSA)


############## beta from data
n = nrow(shuffled.x)
beta.g <- lm(shuffled.y~cbind(rep(1,n),shuffled.x[,c(3,6)])+0)$coef[c(2,3)]
beta.t <- beta.g 
# beta <- round(beta.g,3)
beta.1.g <- lm(shuffled.y~shuffled.x[,3])$coef[2]
# beta.1 <- round(beta.1.g,3)
beta.1.t <- beta.1.g
beta.2.g <- lm(shuffled.y~shuffled.x[,6])$coef[2]
# beta.2 <- round(beta.2.g,3)
beta.2.t <- beta.2.g


theta.ft<-function(vec,beta=beta.t,beta.1=beta.1.t,beta.2=beta.2.t){
p = length(beta)
#create covariance matrix by Sigma = A %*% t(A)
chol.A = matrix(vec,nrow=p)
sigma = chol.A%*%t(chol.A)
X <- mvrnorm(p*20,mu=rep(0,p),Sigma=sigma,empirical=T)
X <- scale(X,scale=F)
X.1 = X[,1]
X.2 = X[,2]
A = rbind(t(X),t(X))
b.1 = t(X)%*%X%*%beta
b.2 = t(X.1)%*%X.1%*%beta.1
b.3 = t(X.2)%*%X.2%*%beta.2
b = as.vector(rbind(b.1,b.2,b.3))
b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
D_X = t(b-b_p)%*%(b-b_p)
theta = rad2deg(acos(cosine(b,b_p)))
return(theta)
}



start.time <- Sys.time()
out <- GenSA(par=c(1,0,0,1),lower=rep(-100,4),upper=rep(100,4),fn=theta.ft,control=list(threshold.stop = 1e-8,max.time=1200))
end.time <- Sys.time()
running.time <- end.time - start.time
running.time

tmp <- matrix(out$par,nrow=2)
res.cov <- tmp%*%t(tmp)


value<- function(sigma){
X <- mvrnorm(p*20,mu=c(0,0),Sigma=sigma,empirical=T)
X <- scale(X,scale=F)
X.1 = X[,1]
X.2 = X[,2]
A = rbind(t(X),t(X))
b.1 = t(X)%*%X%*%beta.t
b.2 = t(X.1)%*%X.1%*%beta.1.t
b.3 = t(X.2)%*%X.2%*%beta.2.t
b = as.vector(rbind(b.1,b.2,b.3))
# qr(cbind(A,b))$rank
# qr(A)$rank
b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
D_X = t(b-b_p)%*%(b-b_p)
theta = rad2deg(acos(cosine(b,b_p)))
return(theta)
}
value(res.cov)
out$value

out1 <- out
res.cov1 <- res.cov
running.time1 <- running.time

#######use Lt(L)
theta.ft<-function(vec,beta=beta.t,beta.1=beta.1.t,beta.2=beta.2.t){
  p = length(beta)
  #create covariance matrix by Sigma = A %*% t(A)
  chol.A = low.mat(vec)
  sigma = chol.A%*%t(chol.A)
  X <- mvrnorm(p*20,mu=rep(0,p),Sigma=sigma,empirical=T)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  return(theta)
}



start.time <- Sys.time()
out <- GenSA(par=c(1,0,1),lower=rep(-100,3),upper=rep(100,3),fn=theta.ft,control=list(threshold.stop = 1e-8,max.time=1200))
end.time <- Sys.time()
running.time <- end.time - start.time
running.time

tmp <- low.mat(out$par)
res.cov <- tmp%*%t(tmp)


value<- function(sigma){
  X <- mvrnorm(p*20,mu=c(0,0),Sigma=sigma,empirical=T)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta.t
  b.2 = t(X.1)%*%X.1%*%beta.1.t
  b.3 = t(X.2)%*%X.2%*%beta.2.t
  b = as.vector(rbind(b.1,b.2,b.3))
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  return(theta)
}
value(res.cov)
out$value

out11 <- out
res.cov11 <- res.cov
running.time11 <- running.time
running.time1
running.time11
################### beta from data and round

beta.g <- lm(shuffled.y~cbind(rep(1,n),shuffled.x[,c(3,6)])+0)$coef[c(2,3)]
beta.t <- round(beta.g,3)
beta.1.g <- lm(shuffled.y~shuffled.x[,3])$coef[2]
beta.1.t <- round(beta.1.g,3)
beta.2.g <- lm(shuffled.y~shuffled.x[,6])$coef[2]
beta.2.t <- round(beta.2.g,3)


theta.ft<-function(vec,beta=beta.t,beta.1=beta.1.t,beta.2=beta.2.t){
  p = length(beta)
  #create covariance matrix by Sigma = A %*% t(A)
  chol.A = matrix(vec,nrow=p)
  sigma = chol.A%*%t(chol.A)
  X <- mvrnorm(p*20,mu=rep(0,p),Sigma=sigma,empirical=T)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  return(theta)
}



out <- GenSA(par=c(1,0,1,0),lower=rep(-100,4),upper=rep(100,4),fn=theta.ft,control=list(threshold.stop = 1e-8,max.time=1200))
tmp <- matrix(out$par,nrow=2)
res.cov <- tmp%*%t(tmp)


value<- function(sigma){
  X <- mvrnorm(p*20,mu=c(0,0),Sigma=sigma,empirical=T)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta.t
  b.2 = t(X.1)%*%X.1%*%beta.1.t
  b.3 = t(X.2)%*%X.2%*%beta.2.t
  b = as.vector(rbind(b.1,b.2,b.3))
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  return(theta)
}
value(res.cov)
out$value

out2 <- out
res.cov2 <- res.cov



## generate random betas
mu1=runif(1,-10,10)
beta.t <- rnorm(2,mean=mu,sd=1)
beta.r <- beta.t
mu2=runif(1,-10,10)
beta.1.t <- rnorm(1,mean=mu2,1)
beta.r.1 <- beta.1.t
mu3=runif(1,-10,10)
beta.2.t <- rnorm(1,mean=mu3,1)
beta.r.2 <- beta.2.t


theta.ft<-function(vec,beta=beta.t,beta.1=beta.1.t,beta.2=beta.2.t){
  p = length(beta)
  #create covariance matrix by Sigma = A %*% t(A)
  chol.A = matrix(vec,nrow=p)
  sigma = chol.A%*%t(chol.A)
  X <- mvrnorm(p*20,mu=rep(0,p),Sigma=sigma,empirical=T)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  return(theta)
}



out <- GenSA(par=c(1,0,1,0),lower=rep(-100,4),upper=rep(100,4),fn=theta.ft,control=list(threshold.stop = 1e-8,max.time=1200))
tmp <- matrix(out$par,nrow=2)
res.cov <- tmp%*%t(tmp)



value<- function(sigma){
  X <- mvrnorm(p*20,mu=c(0,0),Sigma=sigma,empirical=T)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta.t
  b.2 = t(X.1)%*%X.1%*%beta.1.t
  b.3 = t(X.2)%*%X.2%*%beta.2.t
  b = as.vector(rbind(b.1,b.2,b.3))
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  return(theta)
}
value(res.cov)
out$value

out3 <- out
res.cov3 <- res.cov
