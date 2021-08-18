##annealing
library(GenSA)

#Generate beta hat
air <- read_excel("air.xlsx")

air <- as.data.frame(air)
air <- air[1:20,]
Y.air = as.vector(air[,3])
X.air = as.matrix(air[,-c(1,2,3)])
n = length(Y.air)


################### beta from data and round
beta.g <- lm(Y.air~cbind(rep(1,n),X.air[,c(3,4)])+0)$coef[c(2,3)]
beta.t <- round(beta.g,3)
beta.1.g <- lm(Y.air~cbind(rep(1,n),X.air[,3])+0)$coef[2]
beta.1.t <- round(beta.1.g,3)
beta.2.g <- lm(Y.air~X.air[,4])$coef[2]
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



out <- GenSA(par=c(1,0,1,0),lower=rep(-100,4),upper=rep(100,4),fn=theta.ft,control=list(max.time = 3))

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

