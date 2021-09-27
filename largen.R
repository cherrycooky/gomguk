#Large n 

x1 <- rnorm(109,mean=10,sd=100)
x2 <- rnorm(109,mean=-3,sd=100)
y.o2 <- scale(rnorm(109,mean=0,sd=100),scale=F)
X.o2 <- cbind(x1,x2)
X.o2 <- scale(X.o2,scale=F)
beta.t=lm(y.o2~X.o2+0)$coef
beta.1.t=lm(y.o2~X.o2[,1]+0)$coef
beta.2.t=lm(y.o2~X.o2[,2]+0)$coef
summary(lm(y.o2~X.o2+0))
summary(lm(y.o2~X.o2[,1]+0))
summary(lm(y.o2~X.o2[,2]+0))


out33 <- GenSA(par=c(1,0,1),lower=rep(-100,3),upper=rep(100,3),fn=theta.ft,control=list(threshold.stop = 1e-8,max.time=1200))
tmp <- low.mat(out33$par)
res.cov <- tmp%*%t(tmp)

vals3=c()
for(i in 21:30){
  n=i
  X <- mvrnorm(n,mu=rep(0,p),Sigma=res.cov,empirical=T)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta.t
  b.2 = t(X.1)%*%X.1%*%beta.1.t
  b.3 = t(X.2)%*%X.2%*%beta.2.t
  b = as.vector(rbind(b.1,b.2,b.3))
  find.h <- function(h){
    h=as.vector(h)
    Y = pinv(A)%*%b + (diag(rep(1,n))-pinv(A)%*%A)%*%h
    # var((diag(rep(1,n))-p.A%*%A)%*%h)
    # sqrt(t(Y - X%*%solve(t(X)%*%X)%*%t(X)%*%Y) %*% (Y - X%*%solve(t(X)%*%X)%*%t(X)%*%Y)/(n-p-1) * solve(t(X)%*%X)[1,1])
    
    lm1 <- summary(lm(Y~X+0))
    lm2 <- summary(lm(Y~X[,1]+0))
    lm3 <- summary(lm(Y~X[,2]+0))
    lm1.std.error <- lm1$coefficients[,2]
    lm2.std.error <- lm2$coefficients[,2]
    lm3.std.error <- lm3$coefficients[,2]
    v1 <- (lm1.std.error[1]-0.0907)^2 + (lm1.std.error[2]-0.09446)^2
    v2 <- (lm2.std.error-0.09843)^2
    v3 <- (lm3.std.error-0.09573)^2
    error = v1 + v2 + v3
    return(error = error)
  }
  out3 <- GenSA(par=rep(0,n),lower=rep(-10000,n),upper=rep(10000,n),fn=find.h,control=list(max.time=10))
  vals3[i-2]=out3$value
}

vals3

n=25
X <- mvrnorm(n,mu=rep(0,p),Sigma=res.cov,empirical=T)
X <- scale(X,scale=F)
X.1 = X[,1]
X.2 = X[,2]
A = rbind(t(X),t(X))
b.1 = t(X)%*%X%*%beta.t
b.2 = t(X.1)%*%X.1%*%beta.1.t
b.3 = t(X.2)%*%X.2%*%beta.2.t
b = as.vector(rbind(b.1,b.2,b.3))
find.h <- function(h){
  h=as.vector(h)
  Y = pinv(A)%*%b + (diag(rep(1,n))-pinv(A)%*%A)%*%h
  # var((diag(rep(1,n))-p.A%*%A)%*%h)
  # sqrt(t(Y - X%*%solve(t(X)%*%X)%*%t(X)%*%Y) %*% (Y - X%*%solve(t(X)%*%X)%*%t(X)%*%Y)/(n-p-1) * solve(t(X)%*%X)[1,1])
  
  lm1 <- summary(lm(Y~X+0))
  lm2 <- summary(lm(Y~X[,1]+0))
  lm3 <- summary(lm(Y~X[,2]+0))
  lm1.std.error <- lm1$coefficients[,2]
  lm2.std.error <- lm2$coefficients[,2]
  lm3.std.error <- lm3$coefficients[,2]
  v1 <- (lm1.std.error[1]-0.0907)^2 + (lm1.std.error[2]-0.09446)^2
  v2 <- (lm2.std.error-0.09843)^2
  v3 <- (lm3.std.error-0.09573)^2
  error = v1 + v2 + v3
  return(error = error)
}

out44 <- GenSA(par=rep(0,n),lower=rep(-10000,n),upper=rep(10000,n),fn=find.h,control=list(max.time=10))
out44$value

X3 <- X
y3 = pinv(A)%*%b + (diag(rep(1,n)) - pinv(A)%*%A)%*%out44$par
summary(lm(y3~X3+0))
summary(lm(y3~X3[,1]+0))
summary(lm(y3~X3[,2]+0))


