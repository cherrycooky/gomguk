

border <- function(x,y,m,size,alpha,blur){
  #x = Yhat , y = R0
  #m is number of points
  #size is 0~1 the distance from original
  #blur is sd of blur for border
  #alpha = 1 : points on the frame are equally spaced (0~100)
  range.x<-diff(range(x))
  range.y<-diff(range(y))
  
  x.min<-min(x)-size*range.x
  x.max<-max(x)+size*range.x
  y.min<-min(y)-size*range.y
  y.max<-max(y)+size*range.y
  
  range.x <- x.max - x.min
  range.y <- y.max - y.min
  
  u.alpha <- seq(from = 0 , to = 1, length = m)^alpha
  
  left.x <- rnorm(m, mean=x.min , sd=blur)
  left.y <- y.min +(range.y)*(u.alpha)
  
  right.x<-rnorm(m, mean=x.max , sd=blur)
  right.y<-y.min + (range.y)*(1-u.alpha)
  
  bottom.x<-x.min+(range.x)*(u.alpha)
  bottom.y<-rnorm(m, mean=y.min, sd=blur)
  
  top.x<-x.min+(range.x)*(1-u.alpha)
  top.y<-rnorm(m, mean=y.max, sd=blur)
  
  border.x <- c(bottom.x, top.x, left.x, right.x)
  border.y <- c(bottom.y, top.y, left.y, right.y)
  
  x <- c(x, border.x)
  y <- c(y, border.y)
  
  list(x=as.vector(scale(x)), y=as.vector(scale(y)))
}
#

#
slop <- function(x, y, m, size, alpha, blur){
  tmp <- border(x, y, m, size, alpha, blur)
  yy<-tmp$y
  xx<-tmp$x
  as.vector(lm(yy~xx)$coef[2])
}
#


#Find best alpha value between 0 and 100
#function 1
find.border<-function(x, y, m=100, size=0.05, alpha.min=0.000001, alpha.max=100, blur = 0){
  alpha <- uniroot(slop, c(alpha.min, alpha.max),
                   tol = 10e-14, x=x, y=y, m=m, size=size, blur=blur)
  if(abs(alpha$f.root)>0.001) print("error at size m or alpha range")
#  print(paste("alpha is ",alpha$root))
  border(x, y, m, size, alpha$root,blur)
}
#
#function 2
find.border2<-function(x, y, m=100, size=0.05, alpha.min=0.0001, alpha.max=100, blur=0){
  alpha.domain<-seq(from=alpha.min,to=alpha.max,by=alpha.min)
  index=which(alpha.domain==min(alpha.domain))[1]
  alpha <- alpha.domain[index]
  border(x, y, m, size, alpha,blur)
}



#
alg<-function(Yhat, R0, coef, p){
  #Yhat and R0 are orthogonal
  #p = # of variables
  #j=1~p
  
  j = p
  sd.Yhat<-sd(Yhat)
  sd.R0<-sd(R0)
  tau <- sd.R0
  gam <- sd.Yhat
  n<-length(Yhat)
  beta0 <- 1
  #beta can be any values with abs(beta > 0)
  betas <- runif(p,0,2)
  
  Yhat <- (sd.R0/sd.Yhat) * sqrt(coef/(1-coef))*Yhat
  
  Z<-rnorm(n, mean=0, sd=tau)
  M <- matrix(rnorm(n*p, mean=0, sd=gam), n , p)
  PR0 <- outer(R0,R0,"*")/sum(R0^2)
  
  val = 1
  k=0
  while ( val > 10e-13){

    W <- cbind(rep(1,n),(diag(rep(1,n))-PR0)%*%M)
    A <- W%*%solve(t(W)%*%W)%*%t(W)
    
    t1 <- A%*%Z
    t2 <- PR0 %*% M %*% betas
    t3 <- apply(t(t(M)*betas)[,-j], 1, sum)
    
    M.j <- (1/betas[j]) * (Yhat - beta0 - t1 +t2 -t3)
    new.M <- M
    new.M[,j] <- M.j
    
    delta_k <- (diag(rep(1,n))-PR0)%*%((new.M-M))
    val = max(abs(delta_k))
    M <- new.M
    k=k+1
  }
  W <- cbind(rep(1,n),(diag(rep(1,n))-PR0)%*%M)
  A <- W%*%solve(t(W)%*%W)%*%t(W)
  
  eps <- as.vector(R0 + A%*%Z)
  X <-(diag(rep(1,n)) - PR0)%*%M
  Y <- beta0 + X%*%betas + eps
  
  return(list(X=X, Y=Y, k=k, betas=betas))
}

#make border
pic<-find.border(pikachu$x.list,pikachu$y.list)
Yhat<-pic$x
R0<-pic$y
plot(Yhat,R0,main="Original (orthogonal) Data",type="n")
points(Yhat,R0,pch=20)

data <- alg(Yhat,R0,coef=.05,p=5)
str(data)
reg<-lm(data$Y~data$X)
summary(reg)
plot(reg$fitted,reg$resid,type="n",main="Residual plot from data")
points(reg$fitted,reg$resid,pch=20,col="lightpink3")
abline(lm(reg$resid~reg$fitted),col="pink")




