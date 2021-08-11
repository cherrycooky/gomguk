#Real estate data
#https://archive.ics.uci.edu/ml/datasets/Real+estate+valuation+data+set
###Generate shuffled data
p=2
estate <- read_excel("estate.xlsx")

estate <- as.data.frame(estate)


head(estate)


y.estate <- estate$'Y house price of unit area'
x.estate <- as.matrix(estate[,c(-1,-8)])
corrplot(cor(x.estate),method="number")

shuffled.x <- matrix(nrow=nrow(x.estate), ncol = ncol(x.estate))
colnames(shuffled.x) <- colnames(x.estate)
colnames(x.estate) <- NULL
rownames(x.estate) <- NULL
shuffled.row <- sample(nrow(x.estate),nrow(x.estate))

for(i in 1:nrow(x.estate)){
  shuffled.x[i,] = as.numeric(x.estate[shuffled.row[i],])
}

shuffled.y <- c()
for(i in 1:length(y.estate)){
  shuffled.y[i] = y.estate[shuffled.row[i]]
}

###Simulation5 (higly correlated)
n=nrow(shuffled.x)
beta.g <- lm(shuffled.y~cbind(rep(1,n),shuffled.x[,c(3,6)])+0)$coef[c(2,3)]
beta <- round(beta.g,3)
beta.1.g <- lm(shuffled.y~shuffled.x[,3])$coef[2]
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(shuffled.y~shuffled.x[,6])$coef[2]
beta.2 <- round(beta.2.g,3)


##Case1   Use X ~ N(0,1) (no structure)

theta_list_0_ad1=c()
d_list_0_ad1=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rnorm(n),rnorm(n))
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_ad1[i] = theta
  d_list_0_ad1[i] = D_X
}

hist(theta_list_0_ad1,breaks = 8)


##Case2   Give X covariance structure

theta_list_0_ad2=c()
d_list_0_ad2=c()
cov.mat = cov(head(shuffled.x[,c(3,6)],30))
cov.mat
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=c(0,0),Sigma=cov.mat)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_ad2[i] = theta
  d_list_0_ad2[i] = D_X
}

hist(theta_list_0_ad2,breaks = 8)

##Case3   Give X covariance structure + mean structure

theta_list_0_ad3=c()
d_list_0_ad3=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=apply(head(shuffled.x[,c(3,6)],30),2,mean),Sigma=cov.mat)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_ad3[i] = theta
  d_list_0_ad3[i] = D_X
}

hist(theta_list_0_ad3,breaks = 8)

##Case4   Give X full covariance structure + mean structure
cov.mat = cov(shuffled.x[,c(3,6)])
theta_list_0_ad4=c()
d_list_0_ad4=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=apply(shuffled.x[,c(3,6)],2,mean),Sigma=cov.mat)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_ad4[i] = theta
  d_list_0_ad4[i] = D_X
}

hist(theta_list_0_ad4,breaks = 8)

##Case5   Give X mean structure

theta_list_0_ad5=c()
d_list_0_ad5=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=apply(head(shuffled.x[,c(3,6)],30),2,mean),Sigma=diag(rep(1,2)))
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_ad5[i] = theta
  d_list_0_ad5[i] = D_X
}

hist(theta_list_0_ad5,breaks = 8)

par(mfrow=c(1,5))
hist(theta_list_0_ad1,main = "case1", xlab="theta", breaks=8)
hist(theta_list_0_ad2,main = "case2", xlab="theta", breaks=8)
hist(theta_list_0_ad3,main = "case3", xlab="theta", breaks=8)
hist(theta_list_0_ad4,main = "case4", xlab="theta", breaks=8)
hist(theta_list_0_ad5,main = "case5", xlab="theta", breaks=8)
text(30,500,labels='Simulation 5')



####Simulation6 (not correlated)
n=nrow(shuffled.x)
beta.g <- lm(shuffled.y~cbind(rep(1,n),shuffled.x[,c(2,3)])+0)$coef[c(2,3)]
beta <- round(beta.g,3)
beta.1.g <- lm(shuffled.y~cbind(rep(1,n),shuffled.x[,2])+0)$coef[2]
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(shuffled.y~shuffled.x[,3])$coef[2]
beta.2 <- round(beta.2.g,3)
cov(shuffled.x[,c(2,3)])

##Case1   Use X ~ N(0,1) (no structure)

theta_list_0_adi1=c()
d_list_0_adi1=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rnorm(n),rnorm(n))
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_adi1[i] = theta
  d_list_0_adi1[i] = D_X
}

hist(theta_list_0_adi1,breaks = 8)


##Case2   Give X covariance structure

theta_list_0_adi2=c()
d_list_0_adi2=c()
cov.mat = cov(head(shuffled.x[,c(2,3)],30))
cov.mat
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=c(0,0),Sigma=cov.mat)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_adi2[i] = theta
  d_list_0_adi2[i] = D_X
}

hist(theta_list_0_adi2,breaks = 8)

##Case3   Give X sampled covariance structure + mean structure

theta_list_0_adi3=c()
d_list_0_adi3=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=apply(head(shuffled.x[,c(2,3)],30),2,mean),Sigma=cov.mat)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_adi3[i] = theta
  d_list_0_adi3[i] = D_X
}

hist(theta_list_0_adi3,breaks = 8)

##Case4   Give X full covariance structure + mean structure
cov.mat = cov(shuffled.x[,c(2,3)])
theta_list_0_adi4=c()
d_list_0_adi4=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=apply(shuffled.x[,c(2,3)],2,mean),Sigma=cov.mat)
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_adi4[i] = theta
  d_list_0_adi4[i] = D_X
}

hist(theta_list_0_adi4,breaks = 8)

##Case5   Give X mean structure

theta_list_0_adi5=c()
d_list_0_adi5=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- mvrnorm(n,mu=apply(head(shuffled.x[,c(2,3)],30),2,mean),Sigma=diag(c(1,1)))
  X <- scale(X,scale=F)
  X.1 = X[,1]
  X.2 = X[,2]
  A = rbind(t(X),t(X))
  b.1 = t(X)%*%X%*%beta
  b.2 = t(X.1)%*%X.1%*%beta.1
  b.3 = t(X.2)%*%X.2%*%beta.2
  b = as.vector(rbind(b.1,b.2,b.3))
  
  # qr(cbind(A,b))$rank
  # qr(A)$rank
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  D_X = t(b-b_p)%*%(b-b_p)
  theta = rad2deg(acos(cosine(b,b_p)))
  theta_list_0_adi5[i] = theta
  d_list_0_adi5[i] = D_X
}

hist(theta_list_0_adi5,breaks = 8)

par(mfrow=c(1,5))
hist(theta_list_0_adi1,main = "case1", xlab="theta", breaks=8)
hist(theta_list_0_adi2,main = "case2", xlab="theta", breaks=8)
hist(theta_list_0_adi3,main = "case3", xlab="theta", breaks=8)
hist(theta_list_0_adi4,main = "case4", xlab="theta", breaks=8)
hist(theta_list_0_adi5,main = "case5", xlab="theta", breaks=8)
text(5,300,labels="Simulation 6")
