#change cos function
#11/25 get beta coefficients...
#if result is weird, then tune the parameters of GenSA function
library(GenSA)
library(MASS)
library(pracma)
library(lsa)
library(readxl)

#application
result.gsa <- get.gensa.res(betas,index)
vec.gsa <- result.gsa$par
lowmat.gsa <- low.mat(vec.gsa)
sigma.gsa <- lowmat.gsa%*%t(lowmat.gsa)
sigma.gsa

result11 <- betastotheta_else(vec.gsa,betas,index)
Y.gsa = pinv(result11$A)%*%result11$b + (diag(rep(1,ncol(result11$A)))-pinv(result11$A)%*%result11$A)%*%rnorm(ncol(result11$A))

lm(Y.gsa ~ result11$X + 0)



###
tmp2.gen <- get.gensa.res(betas,index=index)
vec.tmp2 <- tmp2.gen$par
lowmat.tmp2 <- low.mat(vec.tmp2)
sigma.tmp2 <- lowmat.tmp2%*%t(lowmat.tmp2)

n=10000
X.tmp = mvrnorm(n,mu=c(0,0),Sigma=sigma.tmp2,empirical=T)
X.tmp = scale(X.tmp,scale=F)
A.tmp = rbind(t(X.tmp),t(X.tmp))
b.1 = A.tmp[1:2,]%*%t(A.tmp[1:2,])%*%betas[[1]]
b.2 = t(A.tmp[3,])%*%A.tmp[3,]%*%betas[[2]]
b.3 = t(A.tmp[4,])%*%A.tmp[4,]%*%betas[[3]]
b.tmp = c(b.1,b.2,b.3)

Y.tmp = pinv(A.tmp)%*%b.tmp + (diag(rep(1,n))-pinv(A.tmp)%*%A.tmp)%*%rnorm(n)

lm(Y.tmp ~ X.tmp)$coef
lm(Y.tmp ~ X.tmp[,1])$coef
lm(Y.tmp ~ X.tmp[,2])$coef

iris.lm1 = lm(Petal.Width ~ Sepal.Length + Sepal.Width + Petal.Length, data=iris)
iris.lm2 = lm(Petal.Width ~ Sepal.Width + Petal.Length, data=iris)
iris.lm3 = lm(Petal.Width ~ Sepal.Length + Petal.Length, data=iris)
beta.iris = list(iris.lm1,iris.lm2,iris.lm3)

iris.gen <- get.gensa.res(beta.iris)
vec.tmp3<- tmp3.gen$par
lowmat.tmp3 <- low.mat(vec.tmp3)
sigma.tmp3 <- lowmat.tmp3%*%t(lowmat.tmp3)

n=100
X.iris = mvrnorm(n,mu=c(0,0,0),Sigma=sigma.tmp3,empirical=T)
X.iris = scale(X.iris,scale=F)
A.iris = rbind(t(X.iris),t(X.iris[,c(2,3)]))
A.iris = rbind(A.iris,t(X.iris[,c(1,3)]))
b.1.ir = A.iris[1:3,]%*%t(A.iris[1:3,])%*%beta.iris[[1]]$coef[-1]
b.2.ir = (A.iris[c(2,3),])%*%t(A.iris[c(2,3),])%*%beta.iris[[2]]$coef[-1]
b.3.ir = (A.iris[c(1,3),])%*%t(A.iris[c(1,3),])%*%beta.iris[[3]]$coef[-1]
b.iris = c(b.1.ir,b.2.ir,b.3.ir)

Y.iris = pinv(A.iris)%*%b.iris + (diag(rep(1,n))-pinv(A.iris)%*%A.iris)%*%rnorm(n)

beta.iris[[1]]$coef
lm(Y.iris ~ X.iris)$coef

beta.iris[[2]]$coefficients
lm(Y.iris~ X.iris[,c(2,3)])$coef

beta.iris[[3]]$coefficients
lm(Y.iris ~ X.iris[,c(1,3)])$coef

##energy
energy <- read_excel("ENB2012_data.xlsx")

energy <- as.data.frame(energy)

reg1.eng2 <- lm(Y2 ~ X1+X2+X3,data=energy)
reg2.eng2 <- lm(Y2 ~ X1,data=energy)
reg3.eng2 <- lm(Y2 ~ X2,data=energy)
beta.eng2 = list(reg1.eng2,reg2.eng2,reg3.eng2)
beta = beta.eng2

tmp4.gen <- get.gensa.res(beta)


##fish
fish <- read.csv("fish.csv")
reg1.fish <- lm(Weight ~ Height + Width,data=fish)
reg2.fish <- lm(Weight ~ Length1 + Length2,data=fish)
reg3.fish <- lm(Weight ~ Length2 + Length3,data=fish)
beta.fish = list(reg1.fish,reg2.fish,reg3.fish)

fish.gen <- get.gensa.res(beta.fish)
fish.gen$value
vec.fish<- fish.gen$par
lowmat.fish <- low.mat(vec.fish)
sigma.fish <- lowmat.fish%*%t(lowmat.fish)
fish.index
n=30
X.fish = mvrnorm(n,mu=rep(0,nrow(sigma.fish)),Sigma=sigma.fish,empirical=T)
X.fish = scale(X.fish,scale=F)
A.fish = rbind(t(X.fish[,c(1,2)]),t(X.fish[,c(3,4)]))
A.fish = rbind(A.fish,t(X.fish[,c(4,5)]))
b.1.fi = A.fish[c(1,2),]%*%t(A.fish[c(1,2),])%*%beta.fish[[1]]$coef[-1]
b.2.fi = (A.fish[c(3,4),])%*%t(A.fish[c(3,4),])%*%beta.fish[[2]]$coef[-1]
b.3.fi = (A.fish[c(5,6),])%*%t(A.fish[c(5,6),])%*%beta.fish[[3]]$coef[-1]
b.fish = c(b.1.fi,b.2.fi,b.3.fi)

Y.fish = pinv(A.fish)%*%b.fish + (diag(rep(1,n))-pinv(A.fish)%*%A.fish)%*%rnorm(n)

beta.fish[[1]]
lm(Y.fish~X.fish[,1:2])$coef

beta.fish[[2]]
lm(Y.fish~X.fish[,c(3,4)])$coef

beta.fish[[3]]
lm(Y.fish~X.fish[,c(4,5)])$coef



###car
####problematic
head(car)
car <- read.csv("car.csv")
car <- car[,c(10,11,12,13,14,17,19,20,21,22,24,25,26)]

# some.val = c()
# for(i in 1:ncol(car)){
#   some.val[i] = cov(car)[i,i] * cor(car[,ncol(car)],car[,i])
# }
colnames(car)
diag(cov(car))
some.val
reg1.car <- lm(price ~ wheelbase + carheight + carlength + stroke + curbweight,data=car)
reg2.car <- lm(price ~ carheight + carlength + stroke,data=car)
reg3.car <- lm(price ~ wheelbase + carlength,data=car)
beta.car = list(reg1.car,reg2.car,reg3.car)


car.gen <- get.gensa.res(beta.car)
car.gen$value
vec.car<- car.gen$par
  
car.else = betastotheta_else(n=20,vec.car,beta.car)
  
X.car = car.else$X
A.car = car.else$A
b.car = car.else$b
# norm(b.car,type="2")
# norm(A.car%*%pinv(t(A.car)%*%A.car)%*%t(A.car)%*%b.car,type="2")
n = ncol(A.car)
Y.car = pinv(A.car)%*%b.car + (diag(rep(1,n))-pinv(A.car)%*%A.car)%*%rnorm(n)
  
car.index = car.else$index
beta.car[[1]]$coefficients
lm(Y.car ~ X.car[,car.index[[1]]])$coef

beta.car[[2]]
lm(Y.car ~ X.car[,car.index[[2]]])$coef

beta.car[[3]]
lm(Y.car ~ X.car[,car.index[[3]]])$coef




