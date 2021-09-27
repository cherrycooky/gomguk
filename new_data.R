n=10
x1 <- rnorm(10,mean=5,sd=10)
x2 <- rnorm(10,mean=-3,sd=100)
y.o <- scale(rnorm(10,mean=0,sd=10),scale=F)
X.o <- cbind(x1,x2)
X.o <- scale(X.o,scale=F)
beta.t=lm(y.o~X.o+0)$coef
beta.1.t=lm(y.o~X.o[,1]+0)$coef
beta.2.t=lm(y.o~X.o[,2]+0)$coef
summary(lm(y.o~X.o+0))
summary(lm(y.o~X.o[,1]+0))
summary(lm(y.o~X.o[,2]+0))


out22 <- GenSA(par=c(1,0,1),lower=rep(-100,3),upper=rep(100,3),fn=theta.ft,control=list(threshold.stop = 1e-8,max.time=1200))
tmp <- low.mat(out22$par)
res.cov <- tmp%*%t(tmp)

vals2=c()
for(i in 3:15){
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
  v1 <- (lm1.std.error[1]-0.34806)^2 + (lm1.std.error[2]-0.03807)^2
  v2 <- (lm2.std.error-0.3340)^2
  v3 <- (lm3.std.error-0.04588)^2
  error = v1 + v2 + v3
  return(error = error)
}
out3 <- GenSA(par=rep(0,n),lower=rep(-10000,n),upper=rep(10000,n),fn=find.h,control=list(max.time=10))
vals2[i-2]=out3$value
}

vals2

n=10
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
  v1 <- (lm1.std.error[1]-0.34806)^2 + (lm1.std.error[2]-0.03807)^2
  v2 <- (lm2.std.error-0.3340)^2
  v3 <- (lm3.std.error-0.04588)^2
  error = v1 + v2 + v3
  return(error = error)
}

out4 <- GenSA(par=rep(0,n),lower=rep(-10000,n),upper=rep(10000,n),fn=find.h,control=list(max.time=10))
out4$value

X2 <- X
y2 = pinv(A)%*%b + (diag(rep(1,n)) - pinv(A)%*%A)%*%out4$par
summary(lm(y2~X2+0))
summary(lm(y2~X2[,1]+0))
summary(lm(y2~X2[,2]+0))


