###


colnames(car)
diag(cov(car))
some.val
reg1.car <- lm(price ~ wheelbase + carheight + carlength + stroke + curbweight,data=car)
reg2.car <- lm(price ~ carlength + carheight + stroke,data=car)
reg3.car <- lm(price ~ wheelbase + carlength,data=car)
beta.car = list(reg1.car,reg2.car,reg3.car)
beta = list(reg1.car$coef[-1],reg2.car$coef[-1],reg3.car$coef[-1])
car1130.X <- car[,c(1,4,2,8,5)]
car1130.sigma <- cov(car1130.X)
car1130.p = 5

car1130.n = 3 * car1130.p

set.seed(1234)
X <- mvrnorm(car1130.n,mu=rep(0,car1130.p),Sigma=car1130.sigma,empirical=T)
X <- scale(X,scale=F)
#get variables X1 - Xp.
# for(i in 1:ncol(X)){
#   columns = paste("X",i,sep="")
#   assign(columns,X[,i])
# }
index = list(c(1,2,3,4,5),c(3,2,4),c(1,3))
#create A
for(i in 1:length(index)){
  if (i==1){
    A = X[,index[[i]]]
    A = t(A)
    As = paste("A",i,sep="")
    assign(As,A)
  }else{
    As = paste("A",i,sep="")
    assign.mat = X[,index[[i]]]
    assign.mat = t(assign.mat)
    assign(As,assign.mat)
    A = rbind(A,assign.mat)
  }
}

#create b 
for(i in 1:length(index)){
  if(i==1){
    if(length(index[[i]])==1){
      b = t(A[index[[i]],])%*%A[index[[i]],]%*%beta[[i]]
      bs = paste("b",i,sep="")
      assign(bs,b)
    }else{
      b = A[index[[i]],]%*%t(A[index[[i]],])%*%beta[[i]]
      bs = paste("b",i,sep="")
      assign(bs,b)
    }
  }else{
    if(length(index[[i]])==1){
      bs = paste("b",i,sep="")
      assign.b = t(A[index[[i]],])%*%A[index[[i]],]%*%beta[[i]]
      b = rbind(b,assign.b)
      assign(bs,assign.b)
    }else{
      bs = paste("b",i,sep="")
      assign.b = A[index[[i]],]%*%t(A[index[[i]],])%*%beta[[i]]
      b = rbind(b,assign.b)
      assign(bs,assign.b)
    }
  }
}
b <- as.vector(b)
A <- as.matrix(A)
#calculate theta
b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
theta = rad2deg(acos(cosine(b,b_p)))
theta



set.seed(1234)
car1130.Y = pinv(A)%*%b + (diag(rep(1,car1130.n))-pinv(A)%*%A)%*%rnorm(car1130.n)

beta.car[[1]]$coef
lm(car1130.Y ~ X[,index[[1]]])$coef
