##p=2, 3, 4 full combination simulation
##by data fish


# fish <- read.csv("fish.csv")

head(fish)


#fish p = 2, Length1, Length2
reg1.fish <- lm(Weight ~ Length1 + Length2, data= fish)
reg2.fish <- lm(Weight ~ Length1, data= fish)
reg3.fish <- lm(Weight ~ Length2, data= fish)

beta.fish2 = list(reg1.fish,reg2.fish,reg3.fish)

fish2.gen <- get.gensa.res(beta.fish2)
fish2.else <- betastotheta_else(n=6,vec=fish2.gen$par,beta=beta.fish2)

set.seed(2023)
fish2.Y <- pinv(fish2.else$A)%*%fish2.else$b + (diag(rep(1,6))-pinv(fish2.else$A)%*%fish2.else$A)%*%rnorm(6)
summary(lm(fish2.Y ~ fish2.else$X))
summary(reg1.fish)


#what if n=3 & p=2
fish22.else <- betastotheta_else(n=3,vec=fish2.gen$par,beta=beta.fish2)
set.seed(2023)
fish22.Y <- pinv(fish22.else$A)%*%fish22.else$b + (diag(rep(1,3))-pinv(fish22.else$A)%*%fish22.else$A)%*%rnorm(3)
summary(lm(fish22.Y ~ fish22.else$X))$coef
summary(reg1.fish)$coef


# fish p = 3  ,Length1 - Length3
reg1.fish <- lm(Weight ~ Length1 + Length2 + Length3,data=fish)
reg2.fish <- lm(Weight ~ Length1 + Length2,data=fish)
reg3.fish <- lm(Weight ~ Length2 + Length3,data=fish)
reg4.fish <- lm(Weight ~ Length1 + Length3,data=fish)
reg5.fish <- lm(Weight ~ Length1,data=fish)
reg6.fish <- lm(Weight ~ Length2,data=fish)
reg7.fish <- lm(Weight ~ Length3,data=fish)


beta.fish3 = list(reg1.fish,reg2.fish,reg3.fish,reg4.fish,reg5.fish,reg6.fish,reg7.fish)

fish3.gen <- get.gensa.res(beta.fish3)
fish3.else <- betastotheta_else(n=15,vec=fish3.gen$par,beta=beta.fish3)

set.seed(2023)
fish3.Y <- pinv(fish3.else$A)%*%fish3.else$b + (diag(rep(1,15))-pinv(fish3.else$A)%*%fish3.else$A)%*%(100*rnorm(15))
summary(lm(fish3.Y ~ fish3.else$X))
summary(reg1.fish)
summary(lm(fish3.Y ~ fish3.else$X[,c(1,2)]))
summary(reg2.fish)

#what if n=4 & p=3
fish33.else <- betastotheta_else(n=4,vec=fish3.gen$par,beta=beta.fish3)

set.seed(2023)
fish33.Y <- pinv(fish33.else$A)%*%fish33.else$b + (diag(rep(1,4))-pinv(fish33.else$A)%*%fish33.else$A)%*%(rnorm(4))
summary(lm(fish33.Y ~ fish33.else$X))
summary(reg1.fish)


# fish p=4
head(fish)

reg1.fish = lm(Weight ~ Length1 + Length2 + Length3 + Height ,data=fish)
reg2.fish = lm(Weight ~ Length1 + Length2 + Length3 ,data=fish)
reg3.fish = lm(Weight ~ Length1 + Length2 + Height ,data=fish)
reg4.fish = lm(Weight ~ Length1 + Length3 + Height ,data=fish)
reg5.fish = lm(Weight ~ Length2 + Length3 + Height ,data=fish)
reg6.fish = lm(Weight ~ Length1 + Length2 ,data=fish)
reg7.fish = lm(Weight ~ Length1 + Length3 ,data=fish)
reg8.fish = lm(Weight ~ Length1 + Height ,data=fish)
reg9.fish = lm(Weight ~ Length2 + Length3 ,data=fish)
reg10.fish = lm(Weight ~ Length2 + Height ,data=fish)
reg11.fish = lm(Weight ~ Length3 + Height ,data=fish)
reg12.fish = lm(Weight ~ Length1 ,data=fish)
reg13.fish = lm(Weight ~ Length2 ,data=fish)
reg14.fish = lm(Weight ~ Length3 ,data=fish)
reg15.fish = lm(Weight ~ Height ,data=fish)

beta.fish4 = list(reg1.fish,reg2.fish,reg3.fish,reg4.fish,reg5.fish,reg6.fish,reg7.fish,reg8.fish,
                  reg9.fish,reg10.fish,reg11.fish,reg12.fish,reg13.fish,reg14.fish,reg15.fish)


fish4.gen <- get.gensa.res(beta.fish4)
###n=15
fish4.else <- betastotheta_else(n=15,vec=fish4.gen$par,beta=beta.fish4)
first.part.fish4 <- pinv(fish4.else$A)%*%fish4.else$b
second.part.fish4 <- (diag(rep(1,15))-pinv(fish4.else$A)%*%fish4.else$A)%*%(rnorm(15))
fish4.Y <- first.part.fish4 + second.part.fish4
summary(reg1.fish)
summary(lm(fish4.Y ~ fish4.else$X))


###control standard deviation
const.fish4 = 2.6905
fish4.const.Y = first.part.fish4 + const.fish4 * second.part.fish4
summary(reg1.fish)
summary(lm(fish4.const.Y ~ fish4.else$X))


###n=159, control standard deviation
fish4.else.2 <- betastotheta_else(n=159,vec=fish4.gen$par,beta=beta.fish4)
first.part.fish4.2 <- pinv(fish4.else.2$A)%*%fish4.else.2$b
second.part.fish4.2 <- (diag(rep(1,159))-pinv(fish4.else.2$A)%*%fish4.else.2$A)%*%(rnorm(159))
fish4.Y.2 <- first.part.fish4.2 + second.part.fish4.2
summary(reg1.fish)
summary(lm(fish4.Y.2 ~ fish4.else.2$X))
const.fish4.2 = 9.28395
fish4.const.Y.2 = first.part.fish4.2 + const.fish4.2 * second.part.fish4.2
summary(reg1.fish)
summary(lm(fish4.const.Y.2 ~ fish4.else.2$X))


