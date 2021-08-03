
#################Simulation 2########
##########

#Generate beta hat

p=3
air <- read_excel("air.xlsx")

air <- as.data.frame(air)
air <- air[1:20,]
Y.air = as.vector(air[,3])
X.air = as.matrix(air[,-c(1,2,3)])
n = length(Y.air)
# Y.air = Y.air[1:500]
# X.air = X.air[1:500,]

corrplot(cor(X.air),method="number")
#highly correlated
corrplot(cor(X.air[,c(1,3,4)]),method="number")
#low correlated
corrplot(cor(X.air[,c(5,6,11)]),method="number")

#get beta hat from data
###If b(X) \in col(A(X)) and rounded. & X,Y from air data (highly correlated)
beta.g <- lm(Y.air~cbind(rep(1,n),X.air[,c(1,3,4)])+0)$coef
beta <- round(beta.g,3)
beta.1.g <- lm(Y.air~cbind(rep(1,n),X.air[,1])+0)$coef
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(Y.air~X.air[,c(3,4)]+0)$coef
beta.2 <- round(beta.2.g,3)


##Case1   Use X ~ N(0,1) (no structure)

t_t1=c()
t_d1=c()
for(i in 1:1000000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),rnorm(n),rnorm(n),rnorm(n))
  
  A = rbind(t(X[,c(1,2)]),t(X[,c(3,4)]),t(X[,c(1,2)]),t(X[,c(3,4)]))
  b.1 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta.1
  b.2 = t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta.2
  b.3 = t(X[,c(1,2)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(1,2)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  b.4 = t(X[,c(3,4)])%*%X[,c(1,2)]%*%beta[c(1,2)] + t(X[,c(3,4)])%*%X[,c(3,4)]%*%beta[c(3,4)]
  
  b = as.vector(rbind(b.1,b.2,b.3,b.4))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  t_t1[i] = theta
  t_d1[i] = D_X
}

hist(t_t1,breaks = 8)
tmp = t_t1 - min(t_t1)
hist(tmp,freq = F,ylim=c(0,0.4))
chi.data = rchisq(length(tmp),df=mean(tmp))
lines(x=density(x=chi.data),col='red')

