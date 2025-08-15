## Plot mvshapes Xs, 
## Supplementary argument: arr, if True=> plot the arrows of starting point reparametrization 

plot.mvshapes=function(Xs,arr=T, xlim=NULL, ylim=NULL, ...){ 
  if(is.null(xlim)){ 
    xlim=unlist(lapply(Xs, function(xx)lapply(xx, function(x)x[1, ]) ))
    xlim=c(min(xlim), max(xlim))
  }
  if(is.null(ylim)){ 
    ylim=unlist(lapply(Xs, function(xx)lapply(xx, function(x)x[2, ]) ))
    ylim=c(min(ylim), max(ylim))
  }
  
  plot(Xs[[1]][[1]][1, ], Xs[[1]][[1]][2, ], type='l', xlim=xlim, ylim=ylim, col=1, ...)
  for(k in 2:length(Xs)) points(Xs[[k]][[1]][1, ], Xs[[k]][[1]][2, ], type='l', xlim=xlim, ylim=ylim, col=1, ...)
  
  for(j in 1:length(Xs[[1]])){ 
    for(k in 1:length(Xs)){ 
      points(Xs[[k]][[j]][1, ], Xs[[k]][[j]][2, ], type='l', xlim=xlim, ylim=ylim, col=j, ...)
      if(arr){ 
        xs=Xs[[k]][[j]][1, ]
        ys=Xs[[k]][[j]][2, ]
        arrows(x0=xs[1], y0=ys[1], x1=xs[2], y1=ys[2],col='red', ...)
      }
    }
  }
  
}

## Get the ind individuals from Xs: mshapes 
mget_ind=function(Xs, ind){ 
  res=NULL
  for(k in 1:length(Xs)){ 
    res_k=as.shapes(Xs[[k]][ind])
    attr(res_k, 'time')=attr(Xs, 'time')
    res=append(res, list(res_k))
  }
  res=as.mvshapes(res)
  res
}
## Get the m-shape of Xs: Xs^{(m)}
get_mshape=function(Xs, m){ 
  x_m=Xs[[m]]
  attr(x_m, 'time')=attr(Xs, 'time')
  x_m
}
## Rotate the shapes Xs with angle in [0, 2pi]
rotate=function(Xs, angle){ 
  out=NULL
  for(m in 1:length(Xs)){ 
    x_m=get_mshape(Xs, m)
    if(length(angle)==1) angle=rep(angle, length(x_m))
    for(k in 1:length(x_m)){ 
      x_m[[k]]=rot(get_ind(x_m, k), angle[k])
    }
    x_m=as.shapes(x_m)
    out=append(out, list(x_m))
  }
  as.mvshapes(out)  
}

## Parametrize multivariate shapes Xs with d_opt(matrix)
parametrize=function(Xs, d_opt){ 
  X_out=NULL
  for(k in 1:length(Xs)){ 
    xs=get_mshape(Xs, k)
    x_out=NULL
    for(j in 1:nrow(d_opt))x_out=append(x_out, list(param(get_ind(xs,j), d_opt[j, k] )) ) 
    X_out=append(X_out, list(as.shapes(x_out) ))
  }
  
  as.mvshapes(X_out)
  
}
## From list of mshape (dd) to mvshapes 
fusion=function(dd){ 
  X_out=NULL
  for(k in 1:length(dd[[1]])){ 
    rk=lapply(1:length(dd), function(x)dd[[x]][[k]][[1]])
    rk=as.shapes(rk)
    attr(rk, 'time')=attr(dd[[1]], 'time')
    X_out=append(X_out, list(rk))
  }
  X_out=as.mvshapes(X_out)
  attr(X_out, 'time')=attr(rk, 'time')
  X_out
}


## Get the Fourier coefficients (dependent to get_alpha function)

mget_alpha=function(X_s, ...){ 
  res=NULL
  for(k in 1:length(X_s)){ 
    res=append(res, list({x_k=X_s[[k]] ;
    attr(x_k, 'time')= attr(X_s, 'time'); 
    get_alpha(x_k, ...)} ))
  }
  res
}
## From m-shapes, get the preshapes 
muntr=function(X_1){ 
  tt=mget_trans(X_1)
  X_0=X_1
  for(k in 1:length(X_1)){ 
    for(i in 1:length(X_0[[k]]) ){ 
      X_0[[k]][[i]]=X_0[[k]][[i]]-matrix(tt[i, ], ncol=1)[, rep(1, length(attr(X_0, 'time')) )]
    }
  }
  rho=mget_norm(X_0)
  for(k in 1:length(X_1)){ 
    for(i in 1:length(X_0[[k]]) ){ 
      X_0[[k]][[i]]=1/rho[i]*X_0[[k]][[i]]
    }
  }
  X_0
  
}
## From m-shapes, get the trans 
mget_trans=function(X_1){ 
  ts=lapply(1:length(X_1), function(k){ x_k=X_1[[k]]; 
  attr(x_k, 'time')=attr(X_1, 'time'); get_trans(x_k)} )
  T_1=apply(as.data.frame(ts)[, seq(1, 2*length(ts), by=2)], mean, MARGIN=1)
  T_2=apply(as.data.frame(ts)[, -seq(1, 2*length(ts), by=2)], mean, MARGIN=1)
  cbind(T_1, T_2)
}

## From m-shapes Xs, get the norms
mget_norm=function(Xs){ 
  rhos=lapply(1:length(Xs), function(k){ x_k=Xs[[k]]; 
  attr(x_k, 'time')=attr(Xs, 'time'); get_norm(x_k)^2})
  sqrt(apply(as.data.frame(rhos), sum, MARGIN=1))
}

## Sum to mshapes X1, X2
madd=function(X1, X2){ 
  X_out=NULL
  for(k in 1:length(X1)){ 
    X_out=append(X_out, list(add(get_mshape(X1, k), get_mshape(X2,k)  )))
  }
  as.mvshapes(X_out) 
}
## multiply X1 by rho (a scalar)
mscal=function(X1, rho){ 
  X_out=NULL
  for(k in 1:length(X1)){ 
    X_out=append(X_out, list(list(scal(get_mshape(X1, k)[[1]], rho))) ) 
  }  
  X_out=as.mvshapes(X_out) 
  X_out
}
