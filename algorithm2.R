library(expm)
library(randcorr)
pic<-find.border(pikachu$x.list,pikachu$y.list)
Yhat<-pic$x
R0<-pic$y
plot(Yhat,R0,pch=".")


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

#t.C : target correlation , t.beta = target coefficient vector
t.C <- randcorr(p)
t.beta <- as.vector(rnorm(p))

#to get negative power of matrix
Q.2 = sqrt(n-1) *eigen(t(X.1)%*%X.1)$vectors %*% diag((sqrt(eigen(t(X.1)%*%X.1)$values))^(-1)) %*% t(eigen(t(X.1)%*%X.1)$vectors) %*% eigen(t.C)$vectors %*% diag(sqrt(eigen(t.C)$values)) %*% t(eigen(t.C)$vectors)

X.2 = X.1 %*% Q.2

Betahat.X2 = (lm(Y~X.2)$coef[-1])

Q.3 = diag(Betahat.X2)%*%solve(diag(t.beta))

X.M.Q = X.2 %*% Q.3

cor(X.M.Q)
lm(Y~X.M.Q)$coef[-1]
t.C  #error at calculating negative power of matrix in Q.2
t.beta

res<-lm(Y~X.M.Q)$res
fitted<-lm(Y~X.M.Q)$fitted

plot(fitted,res,pch=".")




