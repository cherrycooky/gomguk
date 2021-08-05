#Real estate data
#https://archive.ics.uci.edu/ml/datasets/Real+estate+valuation+data+set
###Generate shuffled data
p=3
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
beta.g <- lm(shuffled.y~cbind(rep(1,n),shuffled.x[,c(4,5,6)])+0)$coef
beta <- round(beta.g,3)
beta.1.g <- lm(shuffled.y~cbind(rep(1,n),shuffled.x[,4])+0)$coef
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(shuffled.y~shuffled.x[,c(5,6)]+0)$coef
beta.2 <- round(beta.2.g,3)


##Case1   Use X ~ N(0,1) (no structure)

theta_list_0_s1=c()
d_list_0_s1=c()
for(i in 1:1000){
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
  theta_list_0_s1[i] = theta
  d_list_0_s1[i] = D_X
}

hist(theta_list_0_s1,breaks = 8)


##Case2   Give X covariance structure

theta_list_0_s2=c()
d_list_0_s2=c()
cov.mat = cov(head(shuffled.x[,c(4,5,6)],30))
cov.mat
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=c(0,0,0),Sigma=cov.mat))
  
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
  theta_list_0_s2[i] = theta
  d_list_0_s2[i] = D_X
}

hist(theta_list_0_s2,breaks = 8)

##Case3   Give X covariance structure + mean structure

theta_list_0_s3=c()
d_list_0_s3=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=apply(head(shuffled.x[,c(4,5,6)],30),2,mean),Sigma=cov.mat))
  
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
  theta_list_0_s3[i] = theta
  d_list_0_s3[i] = D_X
}

hist(theta_list_0_s3,breaks = 8)

##Case4   Give X full covariance structure + mean structure
cov.mat = cov(shuffled.x[,c(4,5,6)])
theta_list_0_s4=c()
d_list_0_s4=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=apply(shuffled.x[,c(4,5,6)],2,mean),Sigma=cov.mat))
  
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
  theta_list_0_s4[i] = theta
  d_list_0_s4[i] = D_X
}

hist(theta_list_0_s4,breaks = 8)

##Case5   Give X mean structure

theta_list_0_s5=c()
d_list_0_s5=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=apply(shuffled.x[,c(4,5,6)],2,mean),Sigma=diag(rep(1,3))))
  
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
  theta_list_0_s5[i] = theta
  d_list_0_s5[i] = D_X
}

hist(theta_list_0_s5,breaks = 8)

par(mfrow=c(1,5))
hist(theta_list_0_s1,main = "case1", xlab="theta", breaks=8)
hist(theta_list_0_s2,main = "case2", xlab="theta", breaks=8)
hist(theta_list_0_s3,main = "case3", xlab="theta", breaks=8)
hist(theta_list_0_s4,main = "case4", xlab="theta", breaks=8)
hist(theta_list_0_s5,main = "case5", xlab="theta", breaks=8)
text(30,500,labels='Simulation 5')



####Simulation6 (not correlated)
n=nrow(shuffled.x)
beta.g <- lm(shuffled.y~cbind(rep(1,n),shuffled.x[,c(1,2,3)])+0)$coef
beta <- round(beta.g,3)
beta.1.g <- lm(shuffled.y~cbind(rep(1,n),shuffled.x[,1])+0)$coef
beta.1 <- round(beta.1.g,3)
beta.2.g <- lm(shuffled.y~shuffled.x[,c(2,3)]+0)$coef
beta.2 <- round(beta.2.g,3)
cov(shuffled.x[,c(1,2,3)])

##Case1   Use X ~ N(0,1) (no structure)

theta_list_0_ss1=c()
d_list_0_ss1=c()
for(i in 1:1000){
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
  theta_list_0_ss1[i] = theta
  d_list_0_ss1[i] = D_X
}

hist(theta_list_0_ss1,breaks = 8)


##Case2   Give X covariance structure

theta_list_0_ss2=c()
d_list_0_ss2=c()
cov.mat = cov(head(shuffled.x[,c(1,2,3)],30))
cov.mat
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=c(0,0,0),Sigma=cov.mat))
  
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
  theta_list_0_ss2[i] = theta
  d_list_0_ss2[i] = D_X
}

hist(theta_list_0_ss2,breaks = 8)

##Case3   Give X sampled covariance structure + mean structure

theta_list_0_ss3=c()
d_list_0_ss3=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=apply(head(shuffled.x[,c(1,2,3)],30),2,mean),Sigma=cov.mat))
  
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
  theta_list_0_ss3[i] = theta
  d_list_0_ss3[i] = D_X
}

hist(theta_list_0_ss3,breaks = 8)

##Case4   Give X full covariance structure + mean structure
cov.mat = cov(shuffled.x[,c(1,2,3)])
theta_list_0_ss4=c()
d_list_0_ss4=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=apply(shuffled.x[,c(1,2,3)],2,mean),Sigma=cov.mat))
  
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
  theta_list_0_ss4[i] = theta
  d_list_0_ss4[i] = D_X
}

hist(theta_list_0_ss4,breaks = 8)

##Case5   Give X mean structure

theta_list_0_ss5=c()
d_list_0_ss5=c()
for(i in 1:1000){
  # X <- cbind(rep(1,n),matrix(rnorm(n*p),ncol=p))
  X <- cbind(rep(1,n),mvrnorm(n*p,mu=apply(shuffled.x[,c(1,2,3)],2,mean),Sigma=diag(rep(1,3))))
  
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
  theta_list_0_ss5[i] = theta
  d_list_0_ss5[i] = D_X
}

hist(theta_list_0_ss5,breaks = 8)

par(mfrow=c(1,5))
hist(theta_list_0_ss1,main = "case1", xlab="theta", breaks=8)
hist(theta_list_0_ss2,main = "case2", xlab="theta", breaks=8)
hist(theta_list_0_ss3,main = "case3", xlab="theta", breaks=8)
hist(theta_list_0_ss4,main = "case4", xlab="theta", breaks=8)
hist(theta_list_0_ss5,main = "case5", xlab="theta", breaks=8)
text(5,300,labels="Simulation 6")
