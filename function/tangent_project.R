## Get the exponential map to the tangent space for one mvshape X1 

mget_v=function(mu, X1){
  p_s=0
  for(k in 1:length(X1)){ 
    p_s=p_s+scalarProduct(to_fda(get_mshape(mu, k), 1:2), to_fda( get_mshape(X1, k), 1:2))
  }
  p_s<-round(p_s, 2)
  omega=acos(p_s)
  if(omega==0){ 
    v=mscal(X1, rho = 0)
  }else{ 
    v=madd(mscal(X1, rho = omega/sin(omega)), mscal(mu, rho = -cos(omega)*omega/sin(omega)))
  }
  attr(v, 'time')=attr(X1, 'time')
  v
}

## Get the exponential map to the tangent space for several mvshapes

mget_V=function(MU, Xs){ 
  
  V=NULL
  p=length(Xs)
  n=length(Xs[[1]])
  n_va=paste('v', 1:p, sep='_')
  eval(parse(text=paste0(paste0(n_va, '<-NULL'), collapse=';') ))
  for(i in 1:n){ 
    vx=mget_v(MU, mget_ind(Xs, i) ) 
    for(k in 1:p){ 
      e_txt=paste0(n_va[k], '<-append(', n_va[k], ', vx[[', k, ']])')
      eval(parse(text = e_txt))
    }
  }
  e_f=paste0('Xf=as.mvshapes(list(', paste0( paste0('as.shapes(', n_va, ')'), collapse = ','), '))')
  eval(parse(text = e_f))
  attr(Xf, 'time')=attr(Xs, 'time')
  Xf
}
