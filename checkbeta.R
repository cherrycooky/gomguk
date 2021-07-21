library(pracma)
##### no 
#case1 n=10 p=3 
X <- matrix(rnorm(30),ncol=3)
Y <- rnorm(10)

beta <- lm(Y~X)$coef[-1]
beta.1 <- lm(Y~X[,1])$coef[-1]
beta.2 <- lm(Y~X[,c(2,3)])$coef[-1]

A = rbind(t(X[,1]),t(X[,c(2,3)]),t(X[,1]),t(X[,c(2,3)]))
b.1 = t(X[,1])%*%X[,1]%*%beta.1
b.2 = t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta.2
b.3 = t(X[,1])%*%X[,1]%*%beta[1] + t(X[,1])%*%X[,c(2,3)]%*%beta[c(2,3)]
b.4 = t(X[,c(2,3)])%*%X[,1]%*%beta[1] + t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta[c(2,3)]

b = rbind(b.1,b.2,b.3,b.4)

qr(cbind(A,b))$rank
qr(A)$rank
yhat<-pinv(A)%*%b

lm(yhat~X)$coef[-1]
beta

lm(yhat~X[,1])$coef[-1]
beta.1

lm(yhat~X[,c(2,3)])$coef[-1]
beta.2


###### okay
#case2 n=10 p=3 + scale  (no constant term!!)
X <- scale(matrix(rnorm(30),ncol=3))
Y <- scale(rnorm(10))

beta <- lm(Y~X)$coef[-1]
beta.1 <- lm(Y~X[,1])$coef[-1]
beta.2 <- lm(Y~X[,c(2,3)])$coef[-1]

A = rbind(t(X[,1]),t(X[,c(2,3)]),t(X[,1]),t(X[,c(2,3)]))
b.1 = t(X[,1])%*%X[,1]%*%beta.1
b.2 = t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta.2
b.3 = t(X[,1])%*%X[,1]%*%beta[1] + t(X[,1])%*%X[,c(2,3)]%*%beta[c(2,3)]
b.4 = t(X[,c(2,3)])%*%X[,1]%*%beta[1] + t(X[,c(2,3)])%*%X[,c(2,3)]%*%beta[c(2,3)]

b = rbind(b.1,b.2,b.3,b.4)

qr(cbind(A,b))$rank
qr(A)$rank

yhat<-pinv(A)%*%b

#check - okay
lm(yhat~X)$coef[-1]
beta

lm(yhat~X[,1])$coef[-1]
beta.1

lm(yhat~X[,c(2,3)])$coef[-1]
beta.2


#make arbitrary y & check
rand.y <- pinv(A)%*%b + (diag(rep(1,10)) - pinv(A)%*%A)%*%rnorm(10)
rand.y
lm(rand.y~X)$coef[-1]
beta

lm(rand.y~X[,1])$coef[-1]
beta.1

lm(rand.y~X[,c(2,3)])$coef[-1]
beta.2


####### approximation make funciton!!!!
#case3 cut digits by round(a,l)

X <- scale(matrix(rnorm(30),ncol=3))
Y <- scale(rnorm(10))
cut.digit(X,Y,1)
cut.digit(X,Y,2)
cut.digit(X,Y,3)


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


######## okay
#case4 n=100 p=5 + scale  (no constant term!!)
X <- scale(matrix(rnorm(500),ncol=5))
Y <- scale(rnorm(100))

beta <- lm(Y~X)$coef[-1]
beta.1 <- lm(Y~X[,c(1,2)])$coef[-1]
beta.2 <- lm(Y~X[,c(3,4,5)])$coef[-1]

A = rbind(t(X[,c(1,2)]),t(X[,c(3,4,5)]),t(X[,c(1,2)]),t(X[,c(3,4,5)]))
b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
b.2 = t(X[,c(3,4,5)])%*%X[,c(3,4,5)]%*%beta.2
b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4,5)]%*%beta[c(3,4,5)]
b.4 = t(X[,c(3,4,5)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4,5)])%*%X[,c(3,4,5)]%*%beta[c(3,4,5)]

b = rbind(b.1,b.2,b.3,b.4)

qr(cbind(A,b))$rank
qr(A)$rank

yhat<-pinv(A)%*%b

#check - okay
lm(yhat~X)$coef[-1]
beta

lm(yhat~X[,c(1,2)])$coef[-1]
beta.1

lm(yhat~X[,c(3,4,5)])$coef[-1]
beta.2


#make arbitrary y & check
rand.y <- pinv(A)%*%b + (diag(rep(1,100)) - pinv(A)%*%A)%*%rnorm(100)

lm(rand.y~X)$coef[-1]
beta

lm(rand.y~X[,c(1,2)])$coef[-1]
beta.1

lm(rand.y~X[,c(3,4,5)])$coef[-1]
beta.2


###### okay
#case5 n=10 p=6 + scale  (no constant term!!) (n < 2p case)
X <- scale(matrix(rnorm(60),ncol=6))
Y <- scale(rnorm(10))

beta <- lm(Y~X)$coef[-1]
beta.1 <- lm(Y~X[,c(1,2)])$coef[-1]
beta.2 <- lm(Y~X[,c(3,4,5,6)])$coef[-1]

A = rbind(t(X[,c(1,2)]),t(X[,c(3,4,5,6)]),t(X[,c(1,2)]),t(X[,c(3,4,5,6)]))
b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
b.2 = t(X[,c(3,4,5,6)])%*%X[,c(3,4,5,6)]%*%beta.2
b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4,5,6)]%*%beta[c(3,4,5,6)]
b.4 = t(X[,c(3,4,5,6)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4,5,6)])%*%X[,c(3,4,5,6)]%*%beta[c(3,4,5,6)]

b = rbind(b.1,b.2,b.3,b.4)

qr(cbind(A,b))$rank
qr(A)$rank

yhat<-pinv(A)%*%b

#check - okay
lm(yhat~X)$coef[-1]
beta

lm(yhat~X[,c(1,2)])$coef[-1]
beta.1

lm(yhat~X[,c(3,4,5,6)])$coef[-1]
beta.2


#make arbitrary y & check
rand.y <- pinv(A)%*%b + (diag(rep(1,10)) - pinv(A)%*%A)%*%rnorm(10)

lm(rand.y~X)$coef[-1]
beta

lm(rand.y~X[,c(1,2)])$coef[-1]
beta.1

lm(rand.y~X[,c(3,4,5,6)])$coef[-1]
beta.2


###### of course, no
#case6 n=10 p=10 + scale  (no constant term!!) (n << 2p case)
X <- scale(matrix(rnorm(100),ncol=10))
Y <- scale(rnorm(10))

beta <- lm(Y~X)$coef[-1]
beta.1 <- lm(Y~X[,c(1,2)])$coef[-1]
beta.2 <- lm(Y~X[,c(3,4,5,6,7,8,9,10)])$coef[-1]

A = rbind(t(X[,c(1,2)]),t(X[,c(3,4,5,6,7,8,9,10)]),t(X[,c(1,2)]),t(X[,c(3,4,5,6,7,8,9,10)]))
b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
b.2 = t(X[,c(3,4,5,6,7,8,9,10)])%*%X[,c(3,4,5,6,7,8,9,10)]%*%beta.2
b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4,5,6,7,8,9,10)]%*%beta[c(3,4,5,6,7,8,9,10)]
b.4 = t(X[,c(3,4,5,6,7,8,9,10)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4,5,6,7,8,9,10)])%*%X[,c(3,4,5,6,7,8,9,10)]%*%beta[c(3,4,5,6,7,8,9,10)]


b = rbind(b.1,b.2,b.3,b.4)

qr(cbind(A,b))$rank
qr(A)$rank

yhat<-pinv(A)%*%b

#check - okay
lm(yhat~X)$coef[-1]
beta

lm(yhat~X[,c(1,2)])$coef[-1]
beta.1

lm(yhat~X[,c(3,4,5,6)])$coef[-1]
beta.2


#make arbitrary y & check
rand.y <- pinv(A)%*%b + (diag(rep(1,10)) - pinv(A)%*%A)%*%rnorm(10)

lm(rand.y~X)$coef[-1]
beta

lm(rand.y~X[,c(1,2)])$coef[-1]
beta.1

lm(rand.y~X[,c(3,4,5,6)])$coef[-1]
beta.2




###### okay
#case7 n=10 p=9 + scale  (no constant term!!) (n << 2p case)
X <- scale(matrix(rnorm(90),ncol=9))
Y <- scale(rnorm(10))

beta <- lm(Y~X)$coef[-1]
beta.1 <- lm(Y~X[,c(1,2)])$coef[-1]
beta.2 <- lm(Y~X[,c(3,4,5,6,7,8,9)])$coef[-1]

A = rbind(t(X[,c(1,2)]),t(X[,c(3,4,5,6,7,8,9)]),t(X[,c(1,2)]),t(X[,c(3,4,5,6,7,8,9)]))
b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
b.2 = t(X[,c(3,4,5,6,7,8,9)])%*%X[,c(3,4,5,6,7,8,9)]%*%beta.2
b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4,5,6,7,8,9)]%*%beta[c(3,4,5,6,7,8,9)]
b.4 = t(X[,c(3,4,5,6,7,8,9)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4,5,6,7,8,9)])%*%X[,c(3,4,5,6,7,8,9)]%*%beta[c(3,4,5,6,7,8,9)]


b = rbind(b.1,b.2,b.3,b.4)

qr(cbind(A,b))$rank
qr(A)$rank

yhat<-pinv(A)%*%b

#check - okay
lm(yhat~X)$coef[-1]
beta

lm(yhat~X[,c(1,2)])$coef[-1]
beta.1

lm(yhat~X[,c(3,4,5,6)])$coef[-1]
beta.2


#make arbitrary y & check
rand.y <- pinv(A)%*%b + (diag(rep(1,10)) - pinv(A)%*%A)%*%rnorm(10,100)
rand.y

lm(rand.y~X)$coef[-1]
beta

lm(rand.y~X[,c(1,2)])$coef[-1]
beta.1

lm(rand.y~X[,c(3,4,5,6,7,8,9)])$coef[-1]
beta.2


rand.y <- pinv(A)%*%b + (diag(rep(1,10)) - pinv(A)%*%A)%*%rnorm(10,100,100)
rand.y
random.vector <- rep(c(-10000,1000),5)
rand.y <- pinv(A)%*%b + (diag(rep(1,10)) - pinv(A)%*%A)%*%random.vector
rand.y




###### fake betas 
#case8 n=100 p=2 + scale  (no constant term!!) && application
X <- scale(matrix(rnorm(200),ncol=2))
Y <- scale(rnorm(100))

beta <- c(0.301,-2.091)
beta.1 <- 0.205
beta.2 <- 1.625

A = rbind(t(X[,1]),t(X[,2]),t(X[,1]),t(X[,2]))
b.1 = t(X[,1])%*%X[,1]%*%beta.1
b.2 = t(X[,2])%*%X[,2]%*%beta.2
b.3 = t(X[,1])%*%X[,1]%*%beta[1] + t(X[,1])%*%X[,2]%*%beta[2]
b.4 = t(X[,2])%*%X[,1]%*%beta[1] + t(X[,2])%*%X[,2]%*%beta[2]


b = rbind(b.1,b.2,b.3,b.4)

qr(cbind(A,b))$rank
qr(A)$rank

yhat<-pinv(A)%*%b

#check - okay
lm(yhat~X)$coef[-1]
beta

lm(yhat~X[,1])$coef[-1]
beta.1

lm(yhat~X[,2])$coef[-1]
beta.2


#make arbitrary y & check
rand.y <- pinv(A)%*%b + (diag(rep(1,100)) - pinv(A)%*%A)%*%rnorm(100)
rand.y

lm(rand.y~X)$coef[-1]
beta

lm(rand.y~X[,1])$coef[-1]
beta.1

lm(rand.y~X[,2])$coef[-1]
beta.2

mean(rand.y)



