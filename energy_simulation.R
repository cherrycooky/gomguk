library(readxl)

#Generate beta hat
energy <- read_excel("ENB2012_data.xlsx")

energy <- as.data.frame(energy)

#dim : 768 x 10
dim(energy)


#arbitary regression setting Y2 ~
#index set = (1,2,5), (1,3,4), (3,7,8)

reg1.eng <- lm(Y2 ~ X1+X2+X5,data=energy)
reg2.eng <- lm(Y2 ~ X1+X3+X4,data=energy)
reg3.eng <- lm(Y2 ~ X3+X7+X8,data=energy)
beta.eng = list(reg1.eng,reg2.eng,reg3.eng)


index.eng = list(c(1,2,5),c(1,3,4),c(3,7,8))
p.eng = 8
res.eng <- get.gensa.res(beta.eng)
res.eng$value
vec.eng <- res.eng$par
lowmat.eng <- low.mat(vec.eng)
sigma.eng <- lowmat.eng%*%t(lowmat.eng)
betastotheta(vec.eng,beta.eng)

result22 <- betastotheta_else(vec.eng,beta.eng)
Y.eng = pinv(result22$A)%*%result22$b + (diag(rep(1,ncol(result22$A)))-pinv(result22$A)%*%result22$A)%*%rnorm(ncol(result22$A))

lm(Y.eng ~ result22$X[,c(1,2,3)])
reg1.eng

beta.eng2 = list(reg1.eng$coefficients[-1],reg2.eng$coefficients[-1],reg2.eng$coefficients[-1])
result33 <- betastotheta_else(vec.eng,beta.eng,index=index.eng)
