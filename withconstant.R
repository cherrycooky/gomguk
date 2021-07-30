#######Simulation 1#############
#Create beta hat
n=20
p=3
# X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
X.g <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)))
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

for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  # X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=apply(X.g[,-1],2,mean),Sigma=cov(X.g[,-1])))
  
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

hist(theta_list_0_c3,breaks=8)
par(mfrow=c(1,3))
hist(theta_list_0_c1,main = "case1", xlab="theta", breaks=8)
hist(theta_list_0_c2,main = "case2", xlab="theta", breaks=8)
hist(theta_list_0_c3,main = "case3", xlab="theta", breaks=8)
text(10,250,labels="simulation 1")
cov(X.g[,-1])
cor(X.g[,-1])


############Simulation 2 ##########
#####exactly same as above, but with no round 

#Create beta hat
n=20
p=3
# X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
# X.g <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
# cov(X.g[,-1])
# Y <- rnorm(n,mean=runif(1,-5,5),sd=10)
# coef.Y <- Y
beta.g <- lm(Y~X.g+0)$coef
beta <- beta.g
beta.1.g <- lm(Y~X.g[,c(1,2)]+0)$coef
beta.1 <- beta.1.g
beta.2.g <- lm(Y~X.g[,c(3,4)]+0)$coef
beta.2 <- beta.2.g

##Case1   Use X ~ N(0,1) (no structure)

theta_list_0_c4=c()
d_list_0_c4=c()
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
  theta_list_0_c4[i] = theta
  d_list_0_c4[i] = D_X
}

hist(theta_list_0_c4,breaks = 8)


#Case2 Give X structure


n=20
p=3

theta_list_0_c5=c()
d_list_0_c5=c()
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
  theta_list_0_c5[i] = theta
  d_list_0_c5[i] = D_X
}

hist(theta_list_0_c5,breaks = 8)

#Case3 use exact covariance structure
theta_list_0_c6=c()
d_list_0_c6=c()

for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  # X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=apply(X.g[,-1],2,mean),Sigma=cov(X.g[,-1])))
  
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
  theta_list_0_c6[i] = theta
  d_list_0_c6[i] = D_X
}

hist(theta_list_0_c6,breaks=8)

par(mfrow=c(1,3))
hist(theta_list_0_c4,main = "case1", xlab="theta", breaks=8)
hist(theta_list_0_c5,main = "case2", xlab="theta", breaks=8)
hist(theta_list_0_c6,main = "case3", xlab="theta", breaks=8)
text(10,400,labels="simulation 2")




#################Simulation 3########
##########

#Generate beta hat

p=3
air <- read_excel("air.xlsx")

air <- as.data.frame(air)
air <- air[1:20,]
Y.air = as.vector(air[,3])
X.air = as.matrix(air[,-c(1,2,3)])
n = length(Y.air)
# Y.air = Y.air[1:500]
# X.air = X.air[1:500,]

corrplot(cor(X.air),method="number")
#highly correlated
corrplot(cor(X.air[,c(1,3,4)]),method="number")
#low correlated
corrplot(cor(X.air[,c(5,6,11)]),method="number")

#get beta hat from data
###If b(X) \in col(A(X)) and rounded. & X,Y from air data (highly correlated)
beta.g <- lm(Y.air~cbind(rep(1,n),X.air[,c(1,3,4)])+0)$coef
beta <- round(beta.g,3)
beta.1.g <- lm(Y.air~cbind(rep(1,n),X.air[,1])+0)$coef
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y.air~X.air[,c(3,4)]+0)$coef
beta.2 <- round(beta.2.g,3)


##Case1   Use X ~ N(0,1) (no structure)

theta_list_0_cc1=c()
d_list_0_cc1=c()
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
  theta_list_0_cc1[i] = theta
  d_list_0_cc1[i] = D_X
}

hist(theta_list_0_cc1,breaks = 8)


##Case2   Give X covariance structure

theta_list_0_cc2=c()
d_list_0_cc2=c()
cov.mat = cov(X.air[,c(1,3,4)])
cov.mat
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=c(0,0,0),Sigma=cov.mat))
  
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
  theta_list_0_cc2[i] = theta
  d_list_0_cc2[i] = D_X
}

hist(theta_list_0_cc2,breaks = 8)

##Case3   Give X covariance structure + mean structure

theta_list_0_cc3=c()
d_list_0_cc3=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=apply(X.air[,c(1,3,4)],2,mean),Sigma=cov.mat))
  
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
  theta_list_0_cc3[i] = theta
  d_list_0_cc3[i] = D_X
}

hist(theta_list_0_cc3,breaks = 8)

par(mfrow=c(1,3))
hist(theta_list_0_cc1,main = "case1", xlab="theta", breaks=8)
hist(theta_list_0_cc2,main = "case2", xlab="theta", breaks=8)
hist(theta_list_0_cc3,main = "case3", xlab="theta", breaks=8)
text(40,500,labels="simulation 3")


##############Simulation 4 ############
#get beta hat from data (realtively low correlation)
###If b(X) \in col(A(X)) and rounded. & X,Y from air data (highly correlated)
beta.g <- lm(Y.air~cbind(rep(1,n),X.air[,c(5,6,11)])+0)$coef
beta <- round(beta.g,3)
beta.1.g <- lm(Y.air~cbind(rep(1,n),X.air[,5])+0)$coef
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y.air~X.air[,c(6,11)]+0)$coef
beta.2 <- round(beta.2.g,3)


##Case1   Use X ~ N(0,1) (no structure)

theta_list_0_ccc1=c()
d_list_0_ccc1=c()
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
  theta_list_0_ccc1[i] = theta
  d_list_0_ccc1[i] = D_X
}

hist(theta_list_0_ccc1,breaks = 8)


##Case2   Give X covariance structure

theta_list_0_ccc2=c()
d_list_0_ccc2=c()
cov.mat = cov(X.air[,c(5,6,11)])
cov.mat
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=c(0,0,0),Sigma=cov.mat))
  
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
  theta_list_0_ccc2[i] = theta
  d_list_0_ccc2[i] = D_X
}

hist(theta_list_0_ccc2,breaks = 8)


##Case3   Give X covariance structure and mean structure

theta_list_0_ccc3=c()
d_list_0_ccc3=c()

for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=apply(X.air[,c(5,6,11)],2,mean),Sigma=cov.mat))
  
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
  theta_list_0_ccc3[i] = theta
  d_list_0_ccc3[i] = D_X
}

hist(theta_list_0_ccc3,breaks = 8)

par(mfrow=c(1,3))
hist(theta_list_0_ccc1,main = "case1", xlab="theta", breaks=8)
hist(theta_list_0_ccc2,main = "case2", xlab="theta", breaks=8)
hist(theta_list_0_ccc3,main = "case3", xlab="theta", breaks=8)
text(20,600,labels="simulation 4")


#######Simulation 5 #############
#random fake betas

###If b(X) \notin col(A(X)) beta ~ N
n=20
p=3
beta <- rnorm(4,mean=runif(1,-5,5),sd=5)
beta.1 <- rnorm(2,mean=runif(1,-5,5),sd=5)
beta.2 <- rnorm(2,mean=runif(1,-5,5),sd=5)

d_list_cccc1=c()
theta_list_cccc1=c()
for(i in 1:300){
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
  theta_list_cccc1[i] = theta
  d_list_cccc1[i] = D_X
}

par(mfrow=c(1,1))
hist(theta_list_cccc1, breaks=8, main = "Simulation 5", xlab="theta")
text(70,60,labels="simulation 5")


