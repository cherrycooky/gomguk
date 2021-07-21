library(expm)
library(randcorr)
library(pracma)

ran.val<-function(){
  #case n=1922 p=3 + with intercept term
  n=1922
  p=3
  X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  Y <- rnorm(n)
  # coef.Y <- Y
  beta <<- lm(Y~X+0)$coef
  beta.1 <<- lm(Y~X[,c(1,2)]+0)$coef
  beta.2 <<- lm(Y~X[,c(3,4)]+0)$coef
  
  A = rbind(t(X[,c(1,2)]),t(X[,c(3,4)]),t(X[,c(1,2)]),t(X[,c(3,4)]))
  b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
  b.2 = t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta.2
  b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  b.4 = t(X[,c(3,4)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  
  b = rbind(b.1,b.2,b.3,b.4)
  
  qr(cbind(A,b))$rank
  qr(A)$rank
  
  yhat<-pinv(A)%*%b
  
  #check - okay
  # lm(yhat~X+0)$coef
  # beta
  # 
  # lm(yhat~X[,c(1,2)]+0)$coef
  # beta.1
  # 
  # lm(yhat~X[,c(3,4)]+0)$coef
  # beta.2
  
  
  #make arbitrary y & check
  rand.y <- pinv(A)%*%b + (diag(rep(1,n)) - pinv(A)%*%A)%*%rnorm(n)
  
  coef.Y <- yhat
  
  #####
  
  
  pic<-find.border(pikachu$x.list,pikachu$y.list)
  Yhat<-pic$x
  R0<-pic$y

  
  p = 4 #can be any number
  n = length(Yhat)
  Y = Yhat + R0 
  
  PYhat = Yhat%*%solve(t(Yhat)%*%Yhat)%*%t(Yhat)
  PR0 = R0%*%solve(t(R0)%*%R0)%*%t(R0)
  
  #M is arbitrary n * p-1 matrix
  M = matrix(rnorm(n*(p-1)),nrow=n,ncol=p-1)
  
  I.n = diag(rep(1,n))
  
  #X.M is n*p
  part=(I.n-PR0)%*%(I.n-PYhat)%*%M
  X.M = as.matrix(cbind(Yhat,part))
  
  X.M = scale(X.M)
  
  #Q.1 is arbitrary nonsingular p * p matrix .
  Q.1 = diag(rep(1,p))
  X.1 = X.M %*% Q.1
  
  t.C <- randcorr(p)
  t.beta <- beta
  
  #to get negative power of matrix
  
  Q.2 = sqrt(n-1) * solve(sqrtm(t(X.1)%*%X.1)$B) %*% sqrtm(t.C)$B
  X.2 = X.1 %*% Q.2
  
  Betahat.X2 = (lm(Y~X.2)$coef[-1])
  
  Q.3 = diag(Betahat.X2)%*%solve(diag(t.beta))
  
  X.M.Q = X.2 %*% Q.3
  
  res<-lm(Y~X.M.Q)$res
  fitted<-lm(Y~X.M.Q)$fitted
  
  image_Y <- fitted
  image_X <- X.M.Q
  
  HX = image_X%*%solve(t(image_X)%*%image_X)%*%t(image_X)

  #############
  id.m=diag(rep(1,n))
  first.term = t(id.m-pinv(A)%*%A)%*%(id.m-pinv(A)%*%A)
  second.term = t(id.m-pinv(A)%*%A)%*%(Y-pinv(A)%*%b)
  
  h = pinv(first.term)%*%second.term
  
  new.yy <- pinv(A)%*%b + (diag(rep(1,n)) - pinv(A)%*%A)%*%h
  
  #rand.y
  lm(new.yy~X+0)$coef
  beta
  
  lm(new.yy~X[,c(1,2)]+0)$coef
  beta.1
  
  lm(new.yy~X[,c(3,4)]+0)$coef
  beta.2
  
  
  par(mfrow=c(1,1))
  res.rand<-lm(new.yy~image_X)$res
  fitted.rand<-lm(new.yy~image_X)$fitted
  plot(fitted.rand,res.rand,pch=".")
  
  Y <- new.yy
  X <- X
  return(
    list(abs_diff=sum(abs(image_Y-HX%*%(pinv(A)%*%b+(id.m-pinv(A)%*%A)%*%h))),
    Y = new.yy,
    X = X)
  )
}



repeat{
  k <- ran.val()$abs_diff
  print(k)
  print(beta)
  if (k < 0.65) {
    break
  }
}
lm(Y~X+0)$coef
beta
lm(Y~X[,c(1,2)]+0)$coef
beta.1
lm(Y~X[,c(3,4)]+0)$coef
beta.2




