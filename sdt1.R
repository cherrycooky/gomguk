###test for sd
###11/25

####KYD prof. data
####define beta
beta.t <- c(0.301,-2.091)
beta.1.t <- 0.205
beta.2.t <- 1.625

beta.kyd = list(beta.t,beta.1.t,beta.2.t)
index.kyd = list(c(1,2),c(1),c(2))

sd.t <- c(0.029,0.470)
sd.1.t <- 0.035
sd.2.t <- 1.182
sd.kyd = list(sd.t,sd.1.t,sd.2.t)

kyd.gensa <- get.gensa.res(beta.kyd,index.kyd)
kyd.else <- betastotheta_else(kyd.gensa$par,beta.kyd,index.kyd)

kyd.sig = low.mat(kyd.gensa$par) %*% t(low.mat(kyd.gensa$par))
kyd.res$p
######Create X
n=9
X <- mvrnorm(n,mu=rep(0,kyd.res$p),Sigma=kyd.sig,empirical=T)
X <- scale(X,scale=F)
X.1 = X[,1]
X.2 = X[,2]
A = rbind(t(X),t(X))
b.1 = t(X)%*%X%*%beta.t
b.2 = t(X.1)%*%X.1%*%beta.1.t
b.3 = t(X.2)%*%X.2%*%beta.2.t
b = as.vector(rbind(b.1,b.2,b.3))


###choose h?? rnorm(n) or rep(1,n) -> same result
set.seed(2021)
h.t = rnorm(n)
X1.t <- X
tmp.y1.t = pinv(A)%*%b + (diag(rep(1,n)) - pinv(A)%*%A)%*%(h.t)
tmp.res1 <- summary(lm(tmp.y1.t~X1.t+0))
tmp.res2 <- summary(lm(y1.t~X1.t[,1]+0))
tmp.res3 <- summary(lm(y1.t~X1.t[,2]+0))
const = 0.029/tmp.res1$coefficients[,2][1]
y1.t = pinv(A)%*%b + (diag(rep(1,n)) - pinv(A)%*%A)%*%(const*h.t)
summary(lm(y1.t~X1.t+0))


h.1 = rep(1,n)
X1.t <- X
tmp.y1.t2 = pinv(A)%*%b + (diag(rep(1,n)) - pinv(A)%*%A)%*%(h.1)
k=(t(tmp.y1.t2)%*%(diag(rep(1,n))-X1.t%*%solve(t(X1.t)%*%X1.t)%*%t(X1.t))%*%tmp.y1.t2)/(n-2)
k=as.numeric(k)
solve(t(X1.t)%*%X1.t)*k

tmp.res2 <- summary(lm(tmp.y1.t2~X1.t+0))
const2 = 0.029/tmp.res2$coefficients[,2][1]
y1.t2 = pinv(A)%*%b + (diag(rep(1,n)) - pinv(A)%*%A)%*%(const2*h.1)
summary(lm(y1.t2~X1.t+0))
summary(lm(y1.t2~X1.t[,1]+0))
summary(lm(y1.t2~X1.t[,2]+0))



#by GenSA algorithm.... fix n>p & calculate l2norm.

tmp.func <- function(h,n,beta,index,sigma.vec,target.var.vec){
  n=as.integer(n)
  else.tmp<-betastotheta_else(n=n,beta=beta,index=index,vec=var.vec)
  A <- else.tmp$A
  X <- else.tmp$X
  b <- else.tmp$b
  pinv.A <- pinv(else.tmp$A)
  Y <- pinv.A%*%b + (diag(rep(1,n))-pinv.A%*%A)%*%h
  hat.var.vec <- summary(lm(Y~X))$coefficients[,2][-1]
  l2norm = norm(target.var.vec - hat.var.vec , type="2")
  return(l2norm)
}

tmp.func2 <- function(nh,beta,index,sigma.vec,target.var.vec){
  n = as.integer(nh[1])
  if(n<0){n=1}
  h = nh[2:(n+1)]
  else.tmp<-betastotheta_else(n=n,beta=beta,index=index,vec=sigma.vec)
  A <- else.tmp$A
  X <- else.tmp$X
  b <- else.tmp$b
  pinv.A <- pinv(else.tmp$A)
  Y <- pinv.A%*%b + (diag(rep(1,n))-pinv.A%*%A)%*%h
  hat.var.vec <- summary(lm(Y~X))$coefficients[,2][-1]
  l2norm = norm(target.var.vec - hat.var.vec , type="2")
  return(l2norm)
}


gensa.l2norm <- function(n,beta,index=list(0),target.var.vec,sigma.vec){
  formals(tmp.func)$n <- n
  formals(tmp.func)$beta <- beta
  formals(tmp.func)$index <- index
  formals(tmp.func)$target.var.vec <- target.var.vec
  formals(tmp.func)$sigma.vec <- sigma.vec
  init.vec <- rnorm(n)
  out.ft <- GenSA(par=init.vec,lower=rep(-1000,n),upper=rep(1000,n),fn=tmp.func,control=list(threshold.stop = 1e-8,max.time=120))
  return(out.ft)
}



tmp.values=c()
for(i in 6:30){
  t.res = paste("tmp.result.",i,sep="")
  assign.res = gensa.l2norm(n=i,beta=beta.kyd,index=index.kyd,target.var.vec = c(0.029,0.470),sigma.vec=kyd.gensa$par)
  assign(t.res,assign.res)
  tmp.values[i-5]=assign.res$value
}






#optimize n and h at once.
optim.l2norm <- function(beta,index,target.var.vec,sigma.vec){
  formals(tmp.func2)$beta <- beta
  formals(tmp.func2)$index <- index
  formals(tmp.func2)$target.var.vec <- target.var.vec
  formals(tmp.func2)$sigma.vec <- sigma.vec
  p <- (-1+sqrt(1+8*length(sigma.vec)))/2
  init.n <- 3*p
  init.h <- rnorm(init.n)
  out <- optim(rep(c(init.n,init.h),5),tmp.func2)
  return(out)
}


out.tmp243 <- optim.l2norm(beta=beta.kyd,index=index.kyd,target.var.vec=c(0.029,0.470),sigma.vec=kyd.gensa$par)

######### data iris

iris.values=c()
for(i in 6:30){
  t.res = paste("iris.result.",i,sep="")
  assign.res = gensa.l2norm(n=i,beta=beta.iris,target.var.vec = target.iris,sigma.vec=tmp3.gen$par)
  assign(t.res,assign.res)
  iris.values[i-5]=assign.res$value
}

iris.150 <- gensa.l2norm(n=150,beta=beta.iris,target.var.vec=target.iris,sigma.vec=tmp3.gen$par)
iris.150$value

iris24.else <- betastotheta_else(n=24,iris.gen$par,beta.iris)
iris.Y = pinv(iris24.else$A)%*%iris24.else$b + (diag(rep(1,24))-pinv(iris24.else$A)%*%iris24.else$A)%*%iris.result.24$par

summary(lm(iris.Y ~ iris24.else$X))


ydk.values=c()
####data kyd
for(i in 7:11){
  t.res = paste("ydk.result.",i,sep="")
  assign.res = gensa.l2norm(n=i,beta=beta.kyd,index=index.kyd,target.var.vec = c(0.029,0.470),sigma.vec=kyd.gensa$par)
  assign(t.res,assign.res)
  ydk.values[i-6]=assign.res$value
}

ydk8.else <- betastotheta_else(n=8,kyd.gensa$par,index=index.kyd,beta.kyd)
ydk.Y = pinv(ydk8.else$A)%*%ydk8.else$b + (diag(rep(1,8))-pinv(ydk8.else$A)%*%ydk8.else$A)%*%ydk.result.8$par

summary(lm(ydk.Y ~ ydk8.else$X))


