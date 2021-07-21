library(expm)
library(randcorr)
library(pracma)

#case n=1922 p=3 + with intercept term
n=1922
p=3
X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
Y <- rnorm(n)
# coef.Y <- Y
beta <- lm(Y~X+0)$coef
beta.1 <- lm(Y~X[,c(1,2)]+0)$coef
beta.2 <- lm(Y~X[,c(3,4)]+0)$coef

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
lm(yhat~X+0)$coef
beta

lm(yhat~X[,c(1,2)]+0)$coef
beta.1

lm(yhat~X[,c(3,4)]+0)$coef
beta.2


#make arbitrary y & check
rand.y <- pinv(A)%*%b + (diag(rep(1,n)) - pinv(A)%*%A)%*%rnorm(n)
#rand.y
lm(rand.y~X+0)$coef
beta

lm(rand.y~X[,c(1,2)]+0)$coef
beta.1

lm(rand.y~X[,c(3,4)]+0)$coef
beta.2

coef.Y <- yhat

#####


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
#u.1 <- matrix(runif(p*p,-1,1),nrow=p,ncol=p)
#tmp <- (u.1+t(u.1))/2
#t.C <- tmp - diag(diag(tmp)) + diag(rep(1,p))
t.C <- randcorr(p)
t.beta <- beta

#to get negative power of matrix

Q.2 = sqrt(n-1) * solve(sqrtm(t(X.1)%*%X.1)$B) %*% sqrtm(t.C)$B
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

image_Y <- fitted
image_X <- X.M.Q

# HX = image_X%*%solve(t(image_X)%*%image_X)%*%t(image_X)
# first.term = t((diag(rep(1,n))-pinv(A)%*%A))%*%(diag(rep(1,n))-pinv(A)%*%A) + t((diag(rep(1,n))-pinv(A)%*%A))%*%HX%*%(diag(rep(1,n))-pinv(A)%*%A)
# second.term = t((diag(rep(1,n))-pinv(A)%*%A))%*%(coef.Y - pinv(A)%*%b + HX%*%image_Y - HX%*%pinv(A)%*%b)
# 
# h = pinv(first.term)%*%second.term
# 
# rand.yy <- pinv(A)%*%b + (diag(rep(1,n)) - pinv(A)%*%A)%*%h
# 
# #rand.y
# lm(rand.yy~X+0)$coef
# beta
# 
# lm(rand.yy~X[,c(1,2)]+0)$coef
# beta.1
# 
# lm(rand.yy~X[,c(3,4)]+0)$coef
# beta.2
# 
# 
# par(mfrow=c(1,2))
# res.rand<-lm(rand.yy~image_X)$res
# fitted.rand<-lm(rand.yy~image_X)$fitted
# plot(fitted.rand,res.rand,pch=".")
# lm(res.rand~fitted.rand)
# res<-lm(Y~X.M.Q)$res
# fitted<-lm(Y~X.M.Q)$fitted
# plot(fitted,res,pch=".")

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


par(mfrow=c(1,2))
res.rand<-lm(new.yy~image_X)$res
fitted.rand<-lm(new.yy~image_X)$fitted
plot(fitted.rand,res.rand,pch=".")
lm(res.rand~fitted.rand)
res<-lm(Y~X.M.Q)$res
fitted<-lm(Y~X.M.Q)$fitted
plot(fitted,res,pch=".")


sum(abs(image_Y-HX%*%(pinv(A)%*%b+(id.m-pinv(A)%*%A)%*%h)))

