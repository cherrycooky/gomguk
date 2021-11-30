#make Y by using result of betastotheta_else (X, A, b)
#ex. else.ex = betastotheta_else(vec,beta)
#    ex.Y = makeY(rnorm(else.ex$n),else.ex$X,else.ex$A,else.ex$b)

makeY <- function(h,X,A,b){
  Y = pinv(A)%*%b + (diag(rep(1,ncol(A)))-pinv(A)%*%A)%*%h
  return(Y)
}
