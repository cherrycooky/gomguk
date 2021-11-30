##build GenSA function 
##function for algorithm GenSA...  

get.gensa.res <- function(beta,index=list(0)){
  if(is.vector(beta[[1]])==0){
    coef.name = c()
    for(i in 1:length(beta)){
      coef.name = union(coef.name,names(beta[[i]]$coefficients))
    }
    coef.name = setdiff(coef.name,"(Intercept)")
    p = length(coef.name)
  }else{
    ind = c()
    for(i in 1:length(index)){
      ind = union(ind,index[[i]])
    }
    p = length(ind)
  }
  formals(betastotheta)$beta <- beta
  formals(betastotheta)$index <- index
  k = p*(p+1)/2
  vec.init = rep(1e-05,k)
  # l=0
  # ##to make identity matrix
  # for(i in 1:p){
  #   l=l+i
  #   vec.init[l]=1
  # }
  out.ft <- GenSA(par=vec.init,lower=rep(-1,k),upper=rep(1,k),fn=betastotheta,control=list(threshold.stop = 1e-10,max.time=3000))
  return(out.ft)
}
