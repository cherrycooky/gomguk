library(pracma)
library(expm)
library(readxl)
library(corrplot)
library(mvtnorm)
library(tictoc)
library(lsa)

n = 100
p = 3

#Q1. known : beta , beta*
X = scale(matrix(rnorm(n*p),ncol=p))
Y = scale(rnorm(n))

###If b(X) \in col(A(X)) and rounded. & X ~ N(0,1)
beta.g <- lm(Y~X)$coef[-1]
beta <- round(beta.g,3)
beta.1.g <- lm(Y~X[,1])$coef[-1]
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y~X[,c(2,3)])$coef[-1]
beta.2 <- round(beta.2.g,3)
d_list_0=c()
theta_list_0 = c()
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
  theta_list_0[i] = theta
  d_list_0[i] = D_X
}

###If b(X) \in col(A(X)) and rounded. & X ~ U(0,1)

d_list_0_u=c()
theta_list_0_u = c()
for(i in 1:300){
  X = scale(matrix(runif(n*p),ncol=p))
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
  theta_list_0_u[i] = theta
  d_list_0_u[i] = D_X
}

###If b(X) \in col(A(X)) and rounded. & X ~ N(-1,1)

d_list_0_2=c()
theta_list_0_2=c()
for(i in 1:300){
  X = scale(matrix(rnorm(n*p,-1,1),ncol=p))
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
  theta_list_0_2[i] = theta
  d_list_0_2[i] = D_X
}

###If b(X) \in col(A(X)) and rounded. & X ~ N(0,10)

d_list_0_3=c()
theta_list_0_3=c()
for(i in 1:300){
  X = scale(matrix(rnorm(n*p,0,10),ncol=p))
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
  theta_list_0_3[i] = theta
  d_list_0_3[i] = D_X
}

#https://archive.ics.uci.edu/ml/machine-learning-databases/00360/AirQualityUCI.zip
#Q1_2. known : beta , beta*
air <- read_excel("air.xlsx")

air <- as.data.frame(air)
air <- air[1:500,]
Y.air = as.vector(air[,3])
X.air = as.matrix(air[,-c(1,2,3)])
n = length(Y.air)
# Y.air = Y.air[1:500]
# X.air = X.air[1:500,]

corrplot(cor(X.air),method="number")
corrplot(cor(X.air[,c(1,3,4)]),method="number")
corrplot(cor(X.air[,c(2,10,12)]),method="number")

###If b(X) \in col(A(X)) and rounded. & X,Y from air data (highly correlated)
beta.g <- lm(Y.air~X.air[,c(1,3,4)])$coef[-1]
beta <- round(beta.g,3)
beta.1.g <- lm(Y.air~X.air[,1])$coef[-1]
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y.air~X.air[,c(3,4)])$coef[-1]
beta.2 <- round(beta.2.g,3)
d_list_0_a=c()
theta_list_0_a=c()
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
  theta_list_0_a[i] = theta
  d_list_0_a[i] = D_X
}

###If b(X) \in col(A(X)) and rounded. & X,Y from air data (almost independent)
beta.g <- lm(Y.air~X.air[,c(2,10,12)])$coef[-1]
beta <- round(beta.g,3)
beta.1.g <- lm(Y.air~X.air[,2])$coef[-1]
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y.air~X.air[,c(10,12)])$coef[-1]
beta.2 <- round(beta.2.g,3)
d_list_0_a2=c()
theta_list_0_a2=c()
for(i in 1:300){
  if(i==1){tic("looop")}
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
  theta_list_0_a2[i] = theta
  d_list_0_a2[i] = D_X
  if(i==300){toc()}
}

###If b(X) \in col(A(X)) and rounded. & X,Y from air data (almost independent)
cor.mat = cor(X.air[,c(1,3,4)])
beta.g <- lm(Y.air~X.air[,c(1,3,4)])$coef[-1]
beta <- round(beta.g,3)
beta.1.g <- lm(Y.air~X.air[,1])$coef[-1]
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y.air~X.air[,c(3,4)])$coef[-1]
beta.2 <- round(beta.2.g,3)
d_list_0_a3=c()
theta_list_0_a3=c()
for(i in 1:300){
  if(i==1){tic("looop")}
  X = scale(matrix(rmvnorm(n*p,mean=c(0,0,0),sigma=cor.mat),ncol=p))
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
  theta_list_0_a3[i] = theta
  d_list_0_a3[i] = D_X
  if(i==300){toc()}
}


###If b(X) \notin col(A(X)) beta ~ N(0,0.1^2)
beta <- rnorm(3,sd=0.1)
beta.1 <- rnorm(1,sd=0.1)
beta.2 <- rnorm(2,sd=0.1)

d_list_1=c()
theta_list_1=c()
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
  theta_list_1[i] = theta
  d_list_1[i] = D_X
}

###If b(X) \notin col(A(X)), beta ~ N(0,0.01^2)
beta <- rnorm(3,sd=0.01)
beta.1 <- rnorm(1,sd=0.01)
beta.2 <- rnorm(2,sd=0.01)

d_list_2=c()
theta_list_2=c()
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
  theta_list_2[i] = theta
  d_list_2[i] = D_X
}

###If b(X) \notin col(A(X)), beta ~ N(0,0.001^2)
beta <- rnorm(3,sd=0.001)
beta.1 <- rnorm(1,sd=0.001)
beta.2 <- rnorm(2,sd=0.001)

d_list_3=c()
theta_list_3=c()
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
  theta_list_3[i] = theta
  d_list_3[i] = D_X
}


#Q2.

