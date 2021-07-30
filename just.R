library(MASS)
#######Simulation 1#############
#Create beta hat
n=20
p=3
# X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
X.g <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
cov(X.g[,-1])
Y <- rnorm(n,mean=runif(1,-5,5),sd=10)
# coef.Y <- Y
beta.g <- lm(Y~X.g+0)$coef
beta <- round(beta.g,3)
beta.1.g <- lm(Y~X.g[,c(1,2)]+0)$coef
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y~X.g[,c(3,4)]+0)$coef
beta.2 <- round(beta.2.g,3)

##Case1   Use X ~ N(0,1) (no structure)

theta_list_0_c1=c()
d_list_0_c1=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),rnorm(n),rnorm(n),rnorm(n))
  
  A = rbind(t(X[,c(1,2)]),t(X[,c(3,4)]),t(X[,c(1,2)]),t(X[,c(3,4)]))
  b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
  b.2 = t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta.2
  b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  b.4 = t(X[,c(3,4)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_c1[i] = theta
  d_list_0_c1[i] = D_X
}

hist(theta_list_0_c1,breaks = 8)


##Case2 Give X structure


n=20
p=3

theta_list_0_c2=c()
d_list_0_c2=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
  # X <- cbind(rep(1,n),rmvnorm(n*p,sigma=cov(X.g[,-1])))
  
  A = rbind(t(X[,c(1,2)]),t(X[,c(3,4)]),t(X[,c(1,2)]),t(X[,c(3,4)]))
  b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
  b.2 = t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta.2
  b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  b.4 = t(X[,c(3,4)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_c2[i] = theta
  d_list_0_c2[i] = D_X
}

hist(theta_list_0_c2,breaks = 8)

#Case3 use exact covariance structure
theta_list_0_c3=c()
d_list_0_c3=c()
mat_list_c3=list()

n=20
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  # X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
  X <- cbind(rep(1,n),rmvnorm(n,mean=apply(X.g[,-1],2,mean),sigma=cov(X.g[,-1])))
  mat_list_c3[[i]] = X
  A = rbind(t(X[,c(1,2)]),t(X[,c(3,4)]),t(X[,c(1,2)]),t(X[,c(3,4)]))
  b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
  b.2 = t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta.2
  b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  b.4 = t(X[,c(3,4)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_c3[i] = theta
  d_list_0_c3[i] = D_X
}
summary(theta_list_0_c3)

hist(theta_list_0_c3,breaks=8)
par(mfrow=c(1,3))
hist(theta_list_0_c1,main = "case1", xlab="theta", breaks=8)
hist(theta_list_0_c2,main = "case2", xlab="theta", breaks=8)
hist(theta_list_0_c3,main = "case3", xlab="theta", breaks=8)
cov(X.g[,-1])
cor(X.g[,-1])



theta_list=c()
d_list=c()
mat_list=list()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  # X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
  X <- cbind(rep(1,n),mvrnorm(n,mu = c(0,0,0), Sigma=cov(tmp[,-1]),empirical = F))
  # X <- cbind(rep(1,n),mvrnorm(n,mu = c(0,0,0), Sigma=cov(X.g[,-1]),empirical = F))
  mat_list[[i]] = X
  A = rbind(t(X[,c(1,2)]),t(X[,c(3,4)]),t(X[,c(1,2)]),t(X[,c(3,4)]))
  b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
  b.2 = t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta.2
  b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  b.4 = t(X[,c(3,4)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list[i] = theta
  d_list[i] = D_X
}
summary(theta_list)

# X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
# X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
# X <- cbind(rep(1,n),mvrnorm(n*p,mu = c(0,0,0), Sigma=cov(mat_list[[536]][,-1]),empirical = F))
X <- cbind(rep(1,n),mvrnorm(n*p,mu = apply(tmp,2,mean)[-1], Sigma=cov(tmp[,-1]),empirical = T))

A = rbind(t(X[,c(1,2)]),t(X[,c(3,4)]),t(X[,c(1,2)]),t(X[,c(3,4)]))
b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
b.2 = t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta.2
b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4)]%*%beta[c(3,4)]
b.4 = t(X[,c(3,4)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta[c(3,4)]

b = as.vector(rbind(b.1,b.2,b.3,b.4))

# qr(cbind(A,b))$rank
# qr(A)$rank
b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
D_X = t(b-b_p)%*%(b-b_p)
theta = rad2deg(acos(cosine(b,b_p)))
theta


theta_list_0_c3=c()
d_list_0_c3=c()
mat_list_c3=list()

n=20
for(i in 1:100000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  # X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
  X <- cbind(rep(1,n),rmvnorm(n,mean=apply(X.g[,-1],2,mean),sigma=cov(X.g[,-1])))
  mat_list_c3[[i]] = X
  A = rbind(t(X[,c(1,2)]),t(X[,c(3,4)]),t(X[,c(1,2)]),t(X[,c(3,4)]))
  b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
  b.2 = t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta.2
  b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  b.4 = t(X[,c(3,4)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_c3[i] = theta
  d_list_0_c3[i] = D_X
}
summary(theta_list_0_c3)

tmp <- mat_list_c3[[which.min(theta_list_0_c3)]]
X = tmp
A = rbind(t(X[,c(1,2)]),t(X[,c(3,4)]),t(X[,c(1,2)]),t(X[,c(3,4)]))
b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
b.2 = t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta.2
b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4)]%*%beta[c(3,4)]
b.4 = t(X[,c(3,4)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta[c(3,4)]

b = as.vector(rbind(b.1,b.2,b.3,b.4))
rand.y <- pinv(A)%*%b + (diag(rep(1,n)) - pinv(A)%*%A)%*%rnorm(n)
summary(theta_list_0_c3)
beta.g
lm(rand.y ~ X+0)$coef
