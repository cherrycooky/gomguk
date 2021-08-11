library(pracma)
library(expm)
library(readxl)
library(corrplot)
library(mvtnorm)
library(tictoc)
library(lsa)
#######Simulation 1############# (n=20, p=2)
#Create beta hat
n=20
p=2
# X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
X.g <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)))
cov(X.g[,-1])
Y <- rnorm(n,mean=runif(1,-5,5),sd=10)
# coef.Y <- Y
beta.g <- lm(Y~X.g+0)$coef[c(2,3)]
beta <- round(beta.g,3)
beta.1.g <- lm(Y~X.g[,2])$coef[2]
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y~X.g[,3])$coef[2]
beta.2 <- round(beta.2.g,3)

##Case1   Use X ~ N(0,1) (no structure)

theta_list_a1=c()
d_list_a1=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rnorm(n),rnorm(n))
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_a1[i] = theta
  d_list_a1[i] = D_X
}

hist(theta_list_a1,breaks = 8)


##Case2 Give X structure

theta_list_a2=c()
d_list_a2=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)))
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_a2[i] = theta
  d_list_a2[i] = D_X
}

hist(theta_list_a2,breaks = 8)

#Case3 use exact mean vector and covariance structure
theta_list_a3=c()
d_list_a3=c()

for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  # X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
  X <- cbind(mvrnorm(n*p,mu=apply(X.g[,-1],2,mean),Sigma=cov(X.g[,-1])))
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_a3[i] = theta
  d_list_a3[i] = D_X
}

hist(theta_list_a3,breaks=8)
par(mfrow=c(1,3))
hist(theta_list_a1,main = "case1", xlab="theta", breaks=8)
hist(theta_list_a2,main = "case2", xlab="theta", breaks=8)
hist(theta_list_a3,main = "case3", xlab="theta", breaks=8)
text(4,250,labels="simulation 1")
cov(X.g[,-1])



#################Simulation 2########
##########

#Generate beta hat
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
corrplot(cor(X.air[,c(3,4)]),method="number")
# cov(X.air[,c(1,3,4)])
# det(cov(X.air[,c(1,3,4)]))
#low correlated
corrplot(cor(X.air[,c(6,11)]),method="number")
# cov(X.air[,c(5,6,11)])
# det(cov(X.air[,c(5,6,11)]))
#get beta hat from data
###If b(X) \in col(A(X)) and rounded. & X,Y from air data (highly correlated)
beta.g <- lm(Y.air~cbind(rep(1,n),X.air[,c(3,4)])+0)$coef[c(2,3)]
beta <- round(beta.g,3)
beta.1.g <- lm(Y.air~cbind(rep(1,n),X.air[,3])+0)$coef[2]
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y.air~X.air[,4])$coef[2]
beta.2 <- round(beta.2.g,3)


##Case1   Use X ~ N(0,1) (no structure)

theta_list_aa1=c()
d_list_aa1=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rnorm(n),rnorm(n))
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_aa1[i] = theta
  d_list_aa1[i] = D_X
}

hist(theta_list_aa1,breaks = 8)


##Case2   Give X covariance structure

theta_list_aa2=c()
d_list_aa2=c()
cov.mat = cov(X.air[,c(3,4)])
cov.mat
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=c(0,0),Sigma=cov.mat)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_aa2[i] = theta
  d_list_aa2[i] = D_X
}

hist(theta_list_aa2,breaks = 8)

##Case3   Give X covariance structure + mean structure

theta_list_aa3=c()
d_list_aa3=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=apply(X.air[,c(3,4)],2,mean),Sigma=cov.mat)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_aa3[i] = theta
  d_list_aa3[i] = D_X
}

hist(theta_list_aa3,breaks = 8)

##Case4   Give X + mean structure

theta_list_aa4=c()
d_list_aa4=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=apply(X.air[,c(3,4)],2,mean),Sigma=diag(rep(1,2)))
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_aa4[i] = theta
  d_list_aa4[i] = D_X
}

hist(theta_list_aa4,breaks = 8)

par(mfrow=c(1,4))
hist(theta_list_aa1,main = "case1", xlab="theta", breaks=8)
hist(theta_list_aa2,main = "case2", xlab="theta", breaks=8)
hist(theta_list_aa3,main = "case3", xlab="theta", breaks=8)
hist(theta_list_aa4,main = "case4", xlab="theta", breaks=8)
text(40,350,labels="simulation 2")

##############Simulation 3 ############
#get beta hat from data (realtively low correlation)
###If b(X) \in col(A(X)) and rounded. & X,Y from air data (highly correlated)
beta.g <- lm(Y.air~cbind(rep(1,n),X.air[,c(6,11)])+0)$coef[c(2,3)]
beta <- round(beta.g,3)
beta.1.g <- lm(Y.air~cbind(rep(1,n),X.air[,6])+0)$coef[2]
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y.air~X.air[,11])$coef[2]
beta.2 <- round(beta.2.g,3)


##Case1   Use X ~ N(0,1) (no structure)

theta_list_aaa1=c()
d_list_aaa1=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rnorm(n),rnorm(n))
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_aaa1[i] = theta
  d_list_aaa1[i] = D_X
}

hist(theta_list_aaa1,breaks = 8)


##Case2   Give X covariance structure

theta_list_aaa2=c()
d_list_aaa2=c()
cov.mat = cov(X.air[,c(6,11)])
cov.mat
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=c(0,0),Sigma=cov.mat)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_aaa2[i] = theta
  d_list_aaa2[i] = D_X
}

hist(theta_list_aaa2,breaks = 8)


##Case3   Give X covariance structure and mean structure

theta_list_aaa3=c()
d_list_aaa3=c()

for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=apply(X.air[,c(6,11)],2,mean),Sigma=cov.mat)
  
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_aaa3[i] = theta
  d_list_aaa3[i] = D_X
}

hist(theta_list_aaa3,breaks = 8)

##Case4   Give X mean structure

theta_list_aaa4=c()
d_list_aaa4=c()

for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=apply(X.air[,c(6,11)],2,mean),Sigma=diag(rep(1,2)))
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_aaa4[i] = theta
  d_list_aaa4[i] = D_X
}

hist(theta_list_aaa4,breaks = 8)

par(mfrow=c(1,4))
hist(theta_list_aaa1,main = "case1", xlab="theta", breaks=8)
hist(theta_list_aaa2,main = "case2", xlab="theta", breaks=8)
hist(theta_list_aaa3,main = "case3", xlab="theta", breaks=8)
hist(theta_list_aaa4,main = "case4", xlab="theta", breaks=8)
text(4,270,labels="simulation 3")


#######Simulation 4 #############
#random fake betas

###If b(X) \notin col(A(X)) beta ~ N
n=20
p=2
beta <- rnorm(2,mean=runif(1,-5,5),sd=5)
beta.1 <- rnorm(1,mean=runif(1,-5,5),sd=5)
beta.2 <- rnorm(1,mean=runif(1,-5,5),sd=5)

d_list_aaaa1=c()
theta_list_aaaa1=c()
for(i in 1:1000){
  X <- mvrnorm(n,mu=c(0,0),Sigma=diag(c(1,1)))
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_aaaa1[i] = theta
  d_list_aaaa1[i] = D_X
}

par(mfrow=c(1,1))
hist(theta_list_aaaa1, breaks=8, main = "Simulation 4", xlab="theta")
text(70,60,labels="simulation 4")


#2

d_list_aaaa2=c()
theta_list_aaaa2=c()
for(i in 1:1000){
  X <- mvrnorm(n,mu=c(0,0),Sigma = matrix(c(10,-2,0,1),ncol=2))
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_aaaa2[i] = theta
  d_list_aaaa2[i] = D_X
}

par(mfrow=c(1,2))
hist(theta_list_aaaa1, breaks=8, main = "Case 1", xlab="theta")
hist(theta_list_aaaa2, breaks=8, main = "Case 2", xlab="theta")

