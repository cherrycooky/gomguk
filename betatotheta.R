#Let's make a function for beta coefficients

#we'll get input beta coefficients and index set and full column rank : p

#first create input data for function by using ydk data.

#beta coefficeints = 0.301, -2.091 for (1,2)
#                    0.205 for (1)
#                    1.625 for (2)
#we will use an synthetic data which generates beta coefficients of ydk data
# reg1 <- lm(y1~X1)
# reg2 <- lm(y1~X1[,1])
# reg3 <- lm(y1~X1[,2])
# betas = list(reg1,reg2,reg3)
beta = list(c(0.301,-2.091),c(0.205),c(1.625))

index = list(c(1,2),c(1),c(2))
p = 2



##################
#function for create lower matrix
low.mat <- function(vec){
  x = length(vec)
  n = (-1 + sqrt(8*x+1))/2
  mat = matrix(0, ncol = n, nrow = n)
  mat[upper.tri(mat, diag = TRUE)] <- vec #Upper triangle
  return(t(mat))
}

betastotheta <- function(vec,beta,index=list(0)){
  #beta : list type variable, get two types. type1 : regression result ex. lm(y~X), type2 : just beta coefficeints ex. c(0.301,-2.091)
  #index : index for beta input, must match order with input beta
  #p : full column rank of unknown X
  #vec : use it to create P.S.D matrix(covariance matrix)
  
  #check whether the input is type1 or type2.
  if(is.vector(beta[[1]])==0){
    tempor.betas = list()
    index.set = list()
    coef.name = c()
    for(i in 1:length(beta)){
      coef.name = union(coef.name,names(beta[[i]]$coefficients))
    }
    coef.name = setdiff(coef.name,"(Intercept)")
    for(i in 1:length(beta)){
      tmp.set=c()
      for(j in 1:length(beta[[i]]$coefficients)){
        tmp.set = c(tmp.set,which(coef.name==names(beta[[i]]$coefficients[j])))
      }
      index.set[[i]] = tmp.set
    }
    index = index.set
    p = length(coef.name)
    for(i in 1:length(beta)){
      tempor.betas[[i]] = as.vector(beta[[i]]$coefficients[-1])
    }
    beta = tempor.betas
  }else{
    if(index[[1]][1]==0){
      print("Missing index set")
    }
    ind = c()
    for(i in 1:length(index)){
      ind = union(ind,index[[i]])
    }
    p = length(ind)
  }
  
  #pick n which satisfies n > 2p
  n = 3 * p
  chol.A = low.mat(vec)
  sigma = chol.A%*%t(chol.A)
  set.seed(1234)
  X <- mvrnorm(n,mu=rep(0,p),Sigma=sigma,empirical=T)
  X <- scale(X,scale=F)
  #get variables X1 - Xp.
  # for(i in 1:ncol(X)){
  #   columns = paste("X",i,sep="")
  #   assign(columns,X[,i])
  # }
  
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
  j=1
  #create b 
  for(i in 1:length(index)){
    if(i==1){
      if(length(index[[i]])==1){
        b = t(A[j,])%*%A[j,]%*%beta[[i]]
        bs = paste("b",i,sep="")
        assign(bs,b)
        j=j+1
      }else{
        b = A[j:(j+length(index[[i]])-1),]%*%t(A[j:(j+length(index[[i]])-1),])%*%beta[[i]]
        bs = paste("b",i,sep="")
        assign(bs,b)
        j=j+length(index[[i]])
      }
    }else{
      if(length(index[[i]])==1){
        bs = paste("b",i,sep="")
        assign.b = t(A[j,])%*%A[j,]%*%beta[[i]]
        b = rbind(b,assign.b)
        assign(bs,assign.b)
        j=j+1
      }else{
        bs = paste("b",i,sep="")
        assign.b = A[(j):(j+length(index[[i]])-1),]%*%t(A[(j):(j+length(index[[i]])-1),])%*%beta[[i]]
        b = rbind(b,assign.b)
        assign(bs,assign.b)
        j=j+length(index[[i]])
      }
    }
  }
  b <- as.vector(b)
  A <- as.matrix(A)
  #calculate theta
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  theta = rad2deg(acos(cosine(b,b_p)))
  return(theta=theta)
}

betastotheta(vec=vec.test,beta=beta,index=index)


#After algorithm
betastotheta_else <- function(n,vec,beta,index=list(0)){
  #beta : list type variable, get two types. type1 : regression result ex. lm(y~X), 
  #                                          type2 : just beta coefficeints ex. c(0.301,-2.091)
  #index : index for beta input, must match order with input beta
  #p : full column rank of unknown X
  #vec : use it to create P.S.D matrix(covariance matrix)
  
  #check whether the input is type1 or type2 and find out what p is
  if(is.vector(beta[[1]])==0){
    tempor.betas = list()
    index.set = list()
    coef.name = c()
    for(i in 1:length(beta)){
      coef.name = union(coef.name,names(beta[[i]]$coefficients))
    }
    coef.name = setdiff(coef.name,"(Intercept)")
    for(i in 1:length(beta)){
      tmp.set=c()
      for(j in 1:length(beta[[i]]$coefficients)){
        tmp.set = c(tmp.set,which(coef.name==names(beta[[i]]$coefficients[j])))
      }
      index.set[[i]] = tmp.set
    }
    index = index.set
    p = length(coef.name)
    for(i in 1:length(beta)){
      tempor.betas[[i]] = as.vector(beta[[i]]$coefficients[-1])
    }
    beta = tempor.betas
  }else{
    if(index[[1]][1]==0){
      print("Missing index set")
    }
    ind = c()
    for(i in 1:length(index)){
      ind = union(ind,index[[i]])
    }
    p = length(ind)
  }
  #pick n which satisfies n > p
  n.X=n
  
  # vec = rnorm(p*(p+1)/2)
  chol.A = low.mat(vec)
  sigma = chol.A%*%t(chol.A)
  set.seed(1234)
  X <- mvrnorm(n.X,mu=rep(0,p),Sigma=sigma,empirical=T)
  X <- scale(X,scale=F)
  #get variables X1 - Xp.
  # for(i in 1:ncol(X)){
  #   columns = paste("X",i,sep="")
  #   assign(columns,X[,i])
  # }
  
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
  
  j=1
  #create b 
  for(i in 1:length(index)){
    if(i==1){
      if(length(index[[i]])==1){
        b = t(A[j,])%*%A[j,]%*%beta[[i]]
        bs = paste("b",i,sep="")
        assign(bs,b)
        j=j+1
      }else{
        b = A[j:(j+length(index[[i]])-1),]%*%t(A[j:(j+length(index[[i]])-1),])%*%beta[[i]]
        bs = paste("b",i,sep="")
        assign(bs,b)
        j=j+length(index[[i]])
      }
    }else{
      if(length(index[[i]])==1){
        bs = paste("b",i,sep="")
        assign.b = t(A[j,])%*%A[j,]%*%beta[[i]]
        b = rbind(b,assign.b)
        assign(bs,assign.b)
        j=j+1
      }else{
        bs = paste("b",i,sep="")
        assign.b = A[(j):(j+length(index[[i]])-1),]%*%t(A[(j):(j+length(index[[i]])-1),])%*%beta[[i]]
        b = rbind(b,assign.b)
        assign(bs,assign.b)
        j=j+length(index[[i]])
      }
    }
  }
  b <- as.vector(b)
  A <- as.matrix(A)
  #calculate theta
  b_p = as.vector(A%*%pinv(t(A)%*%A)%*%t(A)%*%b)
  theta = rad2deg(acos(cosine(b,b_p)))
  return(list(theta=theta,X=X,A=A,b=b,index=index,n=ncol(A),p=p))
}
