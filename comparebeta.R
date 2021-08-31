angles= c()
values=c()

for(i in 1:30){
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
tmp.beta <- c(beta.t,beta.1.t,beta.2.t)
############## beta from data
###generate data from d-sphere
u = c(mvrnorm(1,beta.t,diag(rep(1,2))),rnorm(1,beta.1.t,runif(1,1,50)),rnorm(1,beta.2.t,runif(1,1,50)))  # an array of d normally distributed random variables
u = u/norm(u,'2')
tmp.beta2 = u * norm(tmp.beta,'2')
beta.t = tmp.beta2[c(1,2)]
beta.1.t = tmp.beta2[3]
beta.2.t = tmp.beta2[4]
rad2deg(acos(cosine(tmp.beta,tmp.beta2)))
angles[i] = rad2deg(acos(cosine(tmp.beta,tmp.beta2)))


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



out <- GenSA(par=c(1,0,0,1),lower=rep(-100,4),upper=rep(100,4),fn=theta.ft,control=list(threshold.stop = 1e-8,max.time=1200))

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
values[i] = out$value
}

cor(values,angles,method = 'spearman')
