## for (a) shape(s) X, gives X^{(dim)}
get_dim=function(shape, dim){ 
  if(all(class(shape)=='shape') ){ 
    as.matrix(shape[dim, ])
  }else{ 
    matrix(unlist(lapply(1:length(shape), function(x)shape[[x]][dim, ])), nrow=length(shape), byrow = T )  
  }
}


## The Beta(delta) function with M fourier functions 

Beta=function(delta=0, M){ 
  I_M=1:M
  I_M=I_M[(I_M+1)%%2==0 ]
  B=NULL
  for(k in 1:length(I_M)){ 
    B_k=Ort(-pi*(I_M[k]+1)*delta )
    B=append(B, list(B_k))
  }
  as.matrix(Matrix::bdiag(B))
}

## For angle in [0, 2\pi], gives the Orthogonal matrix 
Ort=function(angle){ 
  if(class(angle)=='numeric')class(angle)='rotation'
  
  if(class(angle)=='rotation'){ 
    matrix(c(cos(angle),sin(angle), -sin(angle), cos(angle)), nrow=2)
  }else{ 
    matrix(c(cos(angle),sin(angle), sin(angle), -cos(angle)), nrow=2)
  }
  
}

## rotate a shape of angle theta \in [0, 2\pi]

rot=function(shape, angle){ 
  res=Ort(angle)%*%shape
  class(res)='shape'
  attr(res, 'time')=attr(shape, 'time')
  res
}

## for shapes X, get the observations ind: X_i, i \in ind
get_ind=function(shapes, ind){ 
  if(length(ind)>1){ 
    res=as.shapes(lapply(ind, function(x)shapes[[x]]))
    attr(res, 'time')=attr(shapes, 'time')
    class(res)='shapes'
    res
  }else{ 
    out=shapes[[ind]]
    ti=attr(shapes, 'time')
    as.shape(out, ti)
  }
} 

## This function takes a matrix of 2\times k dimensions: XX and associated ti of k-elements representing t \in [0,1] 
## The output is "shape" object, which has associated functions: plot, etc. 
as.shape=function(XX, ti){ 
  class(XX)='shape'
  attr(XX, 'time')=ti
  XX
}

## Consider the XX as a shapes (list of shapes with one individual) 
as.shape2=function(XX, ti){ 
  OO=list(XX)
  attr(OO, 'time')=ti
  class(OO)='shapes'
  OO
}
## Transform X_t into fda shapes
to_fda=function(X_t, col){
  if(class(X_t)=='shapes'){ 
    ti=attr(X_t, 'time')
    n=length(X_t)
    if(length(col)==1){ 
      funData(ti, t(simplify2array(lapply(1:n, function(x)X_t[[x]][col, ]))))
    }else{
      out=NULL
      for( cc in col){ 
        out=append(out, list(funData(ti,
                                     t(simplify2array(lapply(1:n, function(x)X_t[[x]][cc, ]))))) )
      }
      multiFunData(out)
    }
  }else{ 
    ti=attr(X_t, 'time')
    X_t=append(list(X_t), list(X_t) )
    X_t=as.shapes(X_t)
    n=2
    if(length(col)==1){ 
      funData(ti, t(simplify2array(lapply(1:n, function(x)X_t[[x]][col, ]))))[1]
    }else{
      out=NULL
      for( cc in col){ 
        out=append(out, list(funData(ti,
                                     t(simplify2array(lapply(1:n, function(x)X_t[[x]][cc, ]))))) )
      }
      multiFunData(out)[1]
    }
  }
}

## From shapes, get the translate vectors 
get_trans=function(shapes){ 
  n=length(shapes)
  T_1=unlist(lapply(1:n, function(x)integrate(to_fda(shapes, 1))[x]))
  T_2=unlist(lapply(1:n, function(x)integrate(to_fda(shapes, 2))[x]))
  
  cbind(T_1, T_2)
  
}

## From a shape xx, get the norm 
get_norm=function(xx){ 
  xx=to_fda(xx, 1:2)
  sqrt(scalarProduct(xx, xx))
}

## Parametrize the shape with delta 
param=function(shape, delta){ 
  if(delta<0) delta=h(delta)
  ti=attr(shape, 'time')
  if(is.null(ti)) ti=seq(0, 1, length.out=ncol(shape))
  xi=which.min(abs(delta-ti))
  ind_i=c(xi:(length(ti)-1),1:(xi) ) 
  
  res=shape
  if(!(xi %in%  c(1,length(ti)))  ){ 
    for(i in 1:nrow(res)){ 
      res[i, ]= res[i, ind_i]
    }
  }
  res
}

##  multiply shape by rho (a scalar)
scal=function(shape, rho){ 
  res=rho*shape
  class(res)='shape'
  attr(res, 'time')=attr(shape, 'time')
  res
}
## sum two shapes
add=function(shape_1, shape_2){ 
  res=shape_1
  for(i in 1:length(shape_1)){ 
    
    res[[i]]=shape_1[[i]]+shape_2[[i]]
    attr(res[[i]], 'time')=NULL
    attr(res[[i]], 'class')=NULL
  }
  class(res)=class(shape_1)
  
  res
}
