library(pracma)
library(expm)
library(readxl)
library(corrplot)
library(mvtnorm)
library(tictoc)
library(lsa)

#version3
#change n = 100 ~ 10000

######
n = 100
p = 3

#Q1. known : beta , beta*
X = scale(matrix(rnorm(n*p),ncol=p))
Y = scale(rnorm(n))

###If b(X) \in col(A(X)) and rounded. & X ~ N(0,1)
#round( . , 1)
beta.g <- lm(Y~X)$coef[-1]
beta <- round(beta.g,3)
beta.1.g <- lm(Y~X[,1])$coef[-1]
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y~X[,c(2,3)])$coef[-1]
beta.2 <- round(beta.2.g,3)

d_list_0_n1=c()
theta_list_0_n1 = c()
for(i in 1:300){
  X = scale(matrix(rnorm(n*p),ncol=p))
  A = rbind(t(X[,1]),t(X[,c(2,3)]),t(X[,1]),t(X[,c(2,3)]))
  b.1 = t(X[,1])%*%X[,1]%*%beta.1
  b.2 = t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta.2
  b.3 = t(X[,1])%*%X[,1]%*%beta[1] + t(X[,1])%*%X[,c(2,3)]%*%beta[c(2,3)]
  b.4 = t(X[,c(2,3)])%*%X[,1]%*%beta[1] + t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta[c(2,3)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_n1[i] = theta
  d_list_0_n1[i] = D_X
}



######
n = 1000
p = 3
X = scale(matrix(rnorm(n*p),ncol=p))
Y = scale(rnorm(n))


#Q1. known : beta , beta*
###If b(X) \in col(A(X)) and rounded. & X ~ N(0,1)
#round( . , 1)
beta.g <- lm(Y~X)$coef[-1]
beta <- round(beta.g,3)
beta.1.g <- lm(Y~X[,1])$coef[-1]
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y~X[,c(2,3)])$coef[-1]
beta.2 <- round(beta.2.g,3)

d_list_0_n2=c()
theta_list_0_n2 = c()
for(i in 1:300){
  X = scale(matrix(rnorm(n*p),ncol=p))
  A = rbind(t(X[,1]),t(X[,c(2,3)]),t(X[,1]),t(X[,c(2,3)]))
  b.1 = t(X[,1])%*%X[,1]%*%beta.1
  b.2 = t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta.2
  b.3 = t(X[,1])%*%X[,1]%*%beta[1] + t(X[,1])%*%X[,c(2,3)]%*%beta[c(2,3)]
  b.4 = t(X[,c(2,3)])%*%X[,1]%*%beta[1] + t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta[c(2,3)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_n2[i] = theta
  d_list_0_n2[i] = D_X
}




#######
n = 2000
p = 3
X = scale(matrix(rnorm(n*p),ncol=p))
Y = scale(rnorm(n))

#Q1. known : beta , beta*
###If b(X) \in col(A(X)) and rounded. & X ~ N(0,1)
#round( . , 1)
beta.g <- lm(Y~X)$coef[-1]
beta <- round(beta.g,3)
beta.1.g <- lm(Y~X[,1])$coef[-1]
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y~X[,c(2,3)])$coef[-1]
beta.2 <- round(beta.2.g,3)

d_list_0_n3=c()
theta_list_0_n3 = c()
for(i in 1:300){
  if(i==1){tic('first looop')}
  X = scale(matrix(rnorm(n*p),ncol=p))
  A = rbind(t(X[,1]),t(X[,c(2,3)]),t(X[,1]),t(X[,c(2,3)]))
  b.1 = t(X[,1])%*%X[,1]%*%beta.1
  b.2 = t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta.2
  b.3 = t(X[,1])%*%X[,1]%*%beta[1] + t(X[,1])%*%X[,c(2,3)]%*%beta[c(2,3)]
  b.4 = t(X[,c(2,3)])%*%X[,1]%*%beta[1] + t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta[c(2,3)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_n3[i] = theta
  d_list_0_n3[i] = D_X
  if(i==1){toc()}
}



#######
n = 1000
p = 3
X = scale(matrix(rnorm(n*p),ncol=p))
Y = scale(rnorm(n))

#Q1. known : beta , beta*
###If b(X) \in col(A(X)) and rounded. & X ~ N(0,1)
#round( . , 1)
beta.g <- lm(Y~X)$coef[-1]
beta <- round(beta.g,2)
beta.1.g <- lm(Y~X[,1])$coef[-1]
beta.1 <- round(beta.1.g,2)
beta.2.g <- lm(Y~X[,c(2,3)])$coef[-1]
beta.2 <- round(beta.2.g,2)

d_list_0_n4=c()
theta_list_0_n4 = c()
for(i in 1:300){
  if(i==1){tic('first looop')}
  X = scale(matrix(rnorm(n*p),ncol=p))
  A = rbind(t(X[,1]),t(X[,c(2,3)]),t(X[,1]),t(X[,c(2,3)]))
  b.1 = t(X[,1])%*%X[,1]%*%beta.1
  b.2 = t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta.2
  b.3 = t(X[,1])%*%X[,1]%*%beta[1] + t(X[,1])%*%X[,c(2,3)]%*%beta[c(2,3)]
  b.4 = t(X[,c(2,3)])%*%X[,1]%*%beta[1] + t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta[c(2,3)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_n4[i] = theta
  d_list_0_n4[i] = D_X
  if(i==1){toc()}
}


change_n <- function(n){
######

p = 3
X = scale(matrix(rnorm(n*p),ncol=p))
Y = scale(rnorm(n))


#Q1. known : beta , beta*
###If b(X) \in col(A(X)) and rounded. & X ~ N(0,1)
#round( . , 1)
beta.g <- lm(Y~X)$coef[-1]
beta <- round(beta.g,3)
beta.1.g <- lm(Y~X[,1])$coef[-1]
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y~X[,c(2,3)])$coef[-1]
beta.2 <- round(beta.2.g,3)

d_list_0_n=c()
theta_list_0_n = c()
for(i in 1:300){
  X = scale(matrix(rnorm(n*p),ncol=p))
  A = rbind(t(X[,1]),t(X[,c(2,3)]),t(X[,1]),t(X[,c(2,3)]))
  b.1 = t(X[,1])%*%X[,1]%*%beta.1
  b.2 = t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta.2
  b.3 = t(X[,1])%*%X[,1]%*%beta[1] + t(X[,1])%*%X[,c(2,3)]%*%beta[c(2,3)]
  b.4 = t(X[,c(2,3)])%*%X[,1]%*%beta[1] + t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta[c(2,3)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_n[i] = theta
  d_list_0_n[i] = D_X
}
return(list(beta = beta, beta.1 = beta.1 , beta.2 = beta.2, theta_list =theta_list_0_n, d_list = d_list_0_n))
}


k<-change_n(1000)
