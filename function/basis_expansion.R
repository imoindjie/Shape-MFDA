### Gives the $M$ Fourier basis coefficients of shapes variable 

get_alpha=function(shapes, M=3){ 
  ti=attr(shapes, 'time')
  basis=create.fourier.basis(rangeval = c(min(ti),max(ti)), M)
  sb1=smooth.basis(ti, t(get_dim(shapes, 1)), basis)
  sb2=smooth.basis(ti, t(get_dim(shapes, 2)), basis)
  alpha_1=t(sb1$fd$coefs)
  alpha_2=t(sb2$fd$coefs)
  
  shapes_sm=NULL
  
  for(i in 1:length(shapes)){ 
    ff_1=eval.fd(ti, smooth.basis(ti, shapes[[i]][1, ], basis)$fd)
    ff_2=eval.fd(ti, smooth.basis(ti, shapes[[i]][2, ], basis)$fd)
    
    X_s_i=t(cbind( ff_1, ff_2))
    attr(X_s_i, 'time')=ti
    class(X_s_i)='shape'
    shapes_sm=append(shapes_sm, list(X_s_i))
  }
  shapes_sm=as.shapes(shapes_sm)
  res=list(alpha_1=alpha_1, alpha_2=alpha_2,sm=shapes_sm, basis=basis, fdobj=list(sb1, sb2 ) )
  class(res)='alpha'
  res
  
}
