## Alignment X_shapes with the template mu, number of Fourier functions 

multi_al=function(X_shapes, mu, n_random=2, K=10){ 
  all=NULL
  n_i=length(X_shapes[[1]])
  for(j in 1:n_i){ 
    message(j)
    all=append(all, list(
      ICP(mget_ind(X_shapes, j), mu, n_random , K)
    ))
  }
  xx=fusion(lapply(1:n_i, function(x)all[[x]]$aligned))
  if(is.null(attr(xx, 'time'))){ 
    attr(xx, 'time')=attr(X_shapes, 'time')
  }
  xx
}

## ICP for one shape Xk to the template MU

ICP=function(Xk, MU, n_random=1, K=30){ 
  ## Merge 
  X_out=NULL
  for(k in 1:length(MU)){ 
    rk=list(MU[[k]][[1]], Xk[[k]][[1]])
    rk=as.shapes(rk)
    attr(rk, 'time')=attr(MU, 'time')
    X_out=append(X_out, list(rk))
  }
  X_out=as.mvshapes(X_out)
  attr(X_out, 'time')=attr(MU, 'time')
  ## Basis expansion
  al=mget_alpha(X_out, K )
  
  A=lapply(1:3, function(x)rbind(al[[x]]$alpha_1[2, ], al[[x]]$alpha_2[2, ]) )
  B=lapply(1:3, function(x)rbind(al[[x]]$alpha_1[1, ], al[[x]]$alpha_2[1, ]) )
  
  DS=NULL
  iniz=NULL
  for(nn in 1:n_random){ 
    ss=sample(seq(0, 1, by=0.001), 1)
    ds=altern(A, B, delta =ss)
    DS=append(DS, list(ds) )
    iniz=c(iniz, ss)
  }
  best=which.min(unlist(lapply(1:n_random, function(x)DS[[x]]$min_val)) )
  DD=DS[[best]]
  Xa=param2(rotate(Xk, -DD$opts[1]), 1-DD$opts[-1])
  
  list(aligned=Xa, res_align=DD$opts, iniz=iniz[best])
  
}

## Using the score A of the multivariate to align to B, determine theta and delta 
## by an alternating minimization  

altern=function(A, B, max_iter=100, theta=0, delta=0, verb=F){ 
  d=length(A)
  if(length(delta)==1) delta=rep(delta, d)
  ## Initialization
  for(k in 1:d){ 
    B[[k]][, -1]=B[[k]][, -1]%*%Beta(delta[k], ncol(B[[k]])-1)
    A[[k]]=Ort(-theta)%*%A[[k]]
  }
  B_a=NULL
  A_a=NULL
  for(k in 1:d){ 
    B_a=cbind(B_a, B[[k]])
    A_a=cbind(A_a, A[[k]])
  }
  res_k=norm(A_a-B_a, type='F')
  #res_k=norm_test(A, B)
  for(k in 1:max_iter){ 
    if(verb) cat(paste0(k, '\n' ) )
    res_1=get_theta(A_a,B_a)
    theta=(theta+res_1$theta)%%(2*pi)
    
    A_a=t(res_1$ort)%*%A_a
    for(j in 1:d){ 
      A[[j]]=t(res_1$ort)%*%A[[j]]
    }
    
    if(round(res_k-norm(A_a-B_a, type='F'), 4)==0){ 
      #if(round(res_k-norm_test(A, B), 4)==0){ 
      break
    }else{ 
      res_k=norm(A_a-B_a, type='F')
      if(verb) message(res_k)
      
      B_a=NULL
      for(j in 1:d){ 
        res_2=get_delta(A[[j]][, -1],B[[j]][, -1])
        delta[j]=(delta[j]+res_2$delta)%%1
        B[[j]][, -1]=B[[j]][, -1]%*%t(res_2$ort)
        B_a=cbind(B_a, B[[j]])
      }
      
      res_k=norm(A_a-B_a, type='F')
      
    }
  }
  list(opts=c(theta, delta), min_val=res_k, max_iter=k)
  
}
## From scores A and B, get the optimal angle (Procrustes)

get_theta=function(A, B){ 
  M=crossprod(t(A), t(B) )
  
  m_s=svd(M)
  Cor=diag(1, nrow = nrow(M) )
  R=(m_s$u)%*%t(m_s$v)
  if(round(det(R))==-1){ 
    Cor[nrow(M),nrow(M) ]=-1
    R=(m_s$u)%*% Cor%*%t(m_s$v)
  }
  theta=atan2(-R[1, 2], R[1, 1]) 
  
  list(theta=theta,ort=R, value=obj(theta/(2*pi), R)  )
}


## From scores A, B, get the optimal reparametrization starting point delta 
get_delta=function(A, B){ 
  M_full=crossprod(A, B)
  M=M_full*0
  index=1:nrow(M)
  index=index[! index%%2==0]
  
  for(k in index){ 
    M[k+c(0, 1),k+c(0, 1) ]=M_full[ k+c(0, 1),k+c(0, 1)]
  }
  ## 
  opts=bisection(function(x)grad(x, M), 0, 1)
  
  opts=c(opts, bisection(function(x)grad(x, M), opts+0.001, 1))
  
  ## If the first one is max, the second one is a min. It suffices to have two optimums !
  
  est_delta=opts[which.min(unlist(lapply(opts, function(x)obj(x, M))))]
  min_val=min(unlist(lapply(opts, function(x)obj(x, M))))
  list(delta=est_delta,ort=Beta(-est_delta, ncol(A)),  value=min_val)
}


### The objective function when looking for \delta (here x) i.e obj(x, R)=\sum_{k}\norm{ R_k- Ort((k+1)*pi*x )}_F^2
obj=function(x, R){ 
  index=1:nrow(R)
  index=index[! index%%2==0]
  res=0
  for(k in index){
    res=res+norm(R[k+c(0, 1), k+c(0, 1)]-Ort((k+1) *pi*x), type='F')^2
  }
  res
} 

## Obj(x) derivative function 
grad=function(x, M){ 
  c_1=0
  c_2=0
  index=1:nrow(M)
  index=index[index%%2==1]
  res=0
  for(k in index ){ 
    c_1=c_1+(k+1)*sum(diag(M[k+c(0, 1), k+c(0, 1)]))
    c_2=c_2+(k+1)*sum(diag(M[k+c(0, 1), k+c(0, 1)]%*%Ort(pi/2)))
    res=res+c_1*sin(2*pi*x)+c_2*cos(2*pi*x)
  }
  
  res
  
}

## Uniform reparametrize of a single multivariate planar curve 

param2=function(Xk, delta){ 
  X_out=NULL
  tt=attr(Xk, 'time')
  for(k in 1:length(Xk)){ 
    X_out=append(X_out, list(list(param(as.shape(Xk[[k]][[1]] ,tt), delta[k])))  ) 
  }
  X_out=as.mvshapes(X_out)
  class(X_out)='mvshapes'
  X_out
} 

