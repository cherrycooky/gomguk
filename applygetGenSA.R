##application



result.gsa <- get.gensa.res(betas,index)
vec.gsa <- result.gsa$par
lowmat.gsa <- low.mat(vec.gsa)
sigma.gsa <- lowmat.gsa%*%t(lowmat.gsa)
sigma.gsa

result11 <- betastotheta_else(vec.gsa,betas,index)
Y.gsa = pinv(result11$A)%*%result11$b + (diag(rep(1,ncol(result11$A)))-pinv(result11$A)%*%result11$A)%*%rnorm(ncol(result11$A))

lm(Y.gsa ~ result11$X + 0)


