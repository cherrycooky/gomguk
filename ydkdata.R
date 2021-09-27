####KYD prof. data
####define beta
beta.t <- c(0.301,-2.091)
beta.1.t <- 0.205
beta.2.t <- 1.625

#function for create lower matrix
# low.mat <- function(vec){
#   x = length(vec)
#   n = (-1 + sqrt(8*x+1))/2
#   mat = matrix(0, ncol = n, nrow = n)
#   mat[upper.tri(mat, diag = TRUE)] <- vec #Upper triangle
#   return(t(mat))
# }

###using Lt(L) for sigma
theta.ft<-function(vec,beta=beta.t,beta.1=beta.1.t,beta.2=beta.2.t){
  p = length(beta)
  #create covariance matrix by Sigma = A %*% t(A)
  chol.A = low.mat(vec)
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


start.time <- Sys.time()
out <- GenSA(par=c(1,0,1),lower=rep(-100,3),upper=rep(100,3),fn=theta.ft,control=list(threshold.stop = 1e-8,max.time=1200))
end.time <- Sys.time()
running.time <- end.time - start.time
running.time
tmp <- low.mat(out$par)
res.cov <- tmp%*%t(tmp)

vals1=c()
for(i in 3:11){
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
  vals1[i-2]=out2$value
}
vals1

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

X1 <- X
y1 = pinv(A)%*%b + (diag(rep(1,n)) - pinv(A)%*%A)%*%out2$par
summary(lm(y1~X1+0))
summary(lm(y1~X1[,1]+0))
summary(lm(y1~X1[,2]+0))

