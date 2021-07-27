n=20
p=3
# X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
Y <- rnorm(n)
# coef.Y <- Y
beta.g <- lm(Y~X+0)$coef
beta <- round(beta.g,3)
beta.1.g <- lm(Y~X[,c(1,2)]+0)$coef
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y~X[,c(3,4)]+0)$coef
beta.2 <- round(beta.2.g,3)

theta_list_0_c1=c()
d_list_0_c1=c()
for(i in 1:1000){
# X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
X <- cbind(rep(1,n),rnorm(n,),rnorm(n),rnorm(n))

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

hist(theta_list_0_c1)



change.round <- function(d){
  n=20
  p=3
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),rnorm(n),rnorm(n),rnorm(n))
  Y <- rnorm(n)
  # coef.Y <- Y
  beta.g <- lm(Y~X+0)$coef
  beta <- round(beta.g,d)
  beta.1.g <- lm(Y~X[,c(1,2)]+0)$coef
  beta.1 <- round(beta.1.g,d)
  beta.2.g <- lm(Y~X[,c(3,4)]+0)$coef
  beta.2 <- round(beta.2.g,d)
  
  theta_list_0_c1=c()
  d_list_0_c1=c()
  for(i in 1:1000){
    # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
    X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
    
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
  
  hist(theta_list_0_c1)
}




###no round
n=20
p=3
# X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
Y <- rnorm(n)
# coef.Y <- Y
beta.g <- lm(Y~X+0)$coef
beta <- beta.g
beta.1.g <- lm(Y~X[,c(1,2)]+0)$coef
beta.1 <- beta.1.g
beta.2.g <- lm(Y~X[,c(3,4)]+0)$coef
beta.2 <- beta.2.g

theta_list_0_c2=c()
d_list_0_c2=c()
for(i in 1:100){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
  
  A = rbind(t(X[,c(1,2)]),t(X[,c(3,4)]),t(X[,c(1,2)]),t(X[,c(3,4)]))
  b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
  b.2 = t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta.2
  b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  b.4 = t(X[,c(3,4)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  qr(cbind(A,b))$rank
  qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_c2[i] = theta
  d_list_0_c2[i] = D_X
}

#round 2
beta <- round(beta.g,2)
beta.1 <- round(beta.1.g,2)
beta.2 <- round(beta.2.g,2)

theta_list_0_c3=c()
d_list_0_c3=c()
for(i in 1:100){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
  
  A = rbind(t(X[,c(1,2)]),t(X[,c(3,4)]),t(X[,c(1,2)]),t(X[,c(3,4)]))
  b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
  b.2 = t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta.2
  b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  b.4 = t(X[,c(3,4)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  qr(cbind(A,b))$rank
  qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_c3[i] = theta
  d_list_0_c3[i] = D_X
}


beta <- round(beta.g,2)
beta.1 <- round(beta.1.g,2)
beta.2 <- round(beta.2.g,2)


  X <- cbind(rep(1,n),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-5,5)),rnorm(n,mean=runif(1,-3,3)))
  
  A = rbind(t(X[,c(1,2)]),t(X[,c(3,4)]),t(X[,c(1,2)]),t(X[,c(3,4)]))
  b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
  b.2 = t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta.2
  b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  b.4 = t(X[,c(3,4)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  qr(cbind(A,b))$rank
  qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
PA = A%*%pinv(t(A)%*%A)%*%t(A)
I.n = diag(rep(1,n))
h=rnorm(n,mean=runif(1,-100,100))
(A%*%pinv(A)%*%PA%*%b + A%*%(I.n - pinv(A)%*%A)%*%h - PA%*%b)
2*t(PA%*%b-b) %*%(A%*%pinv(A)%*%PA%*%b + A%*%(I.n - pinv(A)%*%A)%*%h - PA%*%b)
D_X