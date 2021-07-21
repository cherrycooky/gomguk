# ####### approximation make funciton!!!!
# #case3 cut digits by round(a,l)
# 
X <- scale(matrix(rnorm(30),ncol=3))
Y <- scale(rnorm(10))


cut.digit <- function(X,Y,l){
  beta <- round(lm(Y~X)$coef[-1],l)
  beta.1 <- round(lm(Y~X[,1])$coef[-1],l)
  beta.2 <- round(lm(Y~X[,c(2,3)])$coef[-1],l)
  
  beta_r <- lm(Y~X)$coef[-1]
  beta.1_r <- lm(Y~X[,1])$coef[-1]
  beta.2_r <- lm(Y~X[,c(2,3)])$coef[-1]
  
  #round 3
  A = rbind(t(X[,1]),t(X[,c(2,3)]),t(X[,1]),t(X[,c(2,3)]))
  b.1 = t(X[,1])%*%X[,1]%*%beta.1
  b.2 = t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta.2
  b.3 = t(X[,1])%*%X[,1]%*%beta[1] + t(X[,1])%*%X[,c(2,3)]%*%beta[c(2,3)]
  b.4 = t(X[,c(2,3)])%*%X[,1]%*%beta[1] + t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta[c(2,3)]
  
  b = rbind(b.1,b.2,b.3,b.4)
  
  qr(cbind(A,b))$rank
  qr(A)$rank
  
  yhat<-pinv(A)%*%b
  
  #round x
  A_r = rbind(t(X[,1]),t(X[,c(2,3)]),t(X[,1]),t(X[,c(2,3)]))
  b.1_r = t(X[,1])%*%X[,1]%*%beta.1_r
  b.2_r = t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta.2_r
  b.3_r = t(X[,1])%*%X[,1]%*%beta_r[1] + t(X[,1])%*%X[,c(2,3)]%*%beta_r[c(2,3)]
  b.4_r = t(X[,c(2,3)])%*%X[,1]%*%beta_r[1] + t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta_r[c(2,3)]
  
  b_r = rbind(b.1_r,b.2_r,b.3_r,b.4_r)
  
  qr(cbind(A_r,b_r))$rank
  qr(A_r)$rank
  
  yhat_r<-pinv(A_r)%*%b_r
  
  #check - okay
  lm(yhat~X)$coef[-1]
  lm(yhat_r~X)$coef[-1]
  beta_r
  
  lm(yhat~X[,1])$coef[-1]
  lm(yhat_r~X[,1])$coef[-1]
  beta.1_r
  
  lm(yhat~X[,c(2,3)])$coef[-1]
  lm(yhat_r~X[,c(2,3)])$coef[-1]
  beta.2_r
  
  #calculate ||b - b'||
  b = b
  b_p = A%*%pinv(t(A)%*%A)%*%t(A)%*%b
  
  return(t(b-b_p)%*%(b-b_p))
}
cut.digit(X,Y,1)
cut.digit(X,Y,2)
cut.digit(X,Y,3)
