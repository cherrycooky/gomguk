####KYD prof. data


vals=c()
for(i in 3:20){
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
    v1 <- (lm1.std.error[1]-0.029)^2 + (lm1.std.error[2]-0.470)^2
    v2 <- (lm2.std.error-0.035)^2
    v3 <- (lm3.std.error-1.182)^2
    error = v1 + v2 + v3
    return(error = error)
  }
  
  
  out2 <- GenSA(par=rep(0,n),lower=rep(-10000,n),upper=rep(10000,n),fn=find.h,control=list(max.time=10))
  vals[i-2]=out2$value
}
vals

######Create X

n=9
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
  v1 <- (lm1.std.error[1]-0.029)^2 + (lm1.std.error[2]-0.470)^2
  v2 <- (lm2.std.error-0.035)^2
  v3 <- (lm3.std.error-1.182)^2
  error = v1 + v2 + v3
  return(error = error)
}


out2 <- GenSA(par=rep(0,n),lower=rep(-10000,n),upper=rep(10000,n),fn=find.h,control=list(max.time=10))
out2$value


y = pinv(A)%*%b + (diag(rep(1,n)) - pinv(A)%*%A)%*%out2$par
summary(lm(y~X+0))
summary(lm(y~X[,1]+0))
summary(lm(y~X[,2]+0))

