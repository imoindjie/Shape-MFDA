mu_fre=function(Xs, K=30, n_random=1){ 
  n=length(Xs[[1]])
  
  MUs=NULL
  res=NULL
  Align=NULL
  for(xs in 1:n_random){ 
    print(paste0(xs, '/', n_random))
    MU=mget_ind(Xs, sample(1:n, 1))
    
    cat('Alignement')
    message('')
    all=NULL
    for(i in 1:n){ 
      message(i)
      all=append(all, list(
        ICP(mget_ind(Xs, i), MU, n_random = 1, K=K)
      ))
    }
    X_tilde=fusion(lapply(1:n, function(x)all[[x]]$aligned))
    if(is.null(attr(X_tilde, 'time'))){ 
      attr(X_tilde, 'time')=attr(Xs, 'time')
    }
    cat('Estimation Frechet Mean')
    MU=frech_mean(X_tilde, K)
    cat('\n')
    cat('Alignement final')
    cat('\n')
    all=NULL
    for(i in 1:n){ 
      message(i)
      all=append(all, list(
        ICP(mget_ind(Xs, i), MU, n_random = 1, K=K)
      ))
    }
    X_tilde=fusion(lapply(1:n, function(x)all[[x]]$aligned))
    if(is.null(attr(X_tilde, 'time'))){ 
      attr(X_tilde, 'time')=attr(Xs, 'time')
    }
    
    
    V=lapply(1:n, function(x)mget_v(MU, mget_ind(X_tilde, x) ))
    f_1=NULL
    f_2=NULL
    for(k in 1:3 ){
      xx=t(as.matrix(as.data.frame(lapply(V,
                                          function(x)x[[k]][[1]][1, ])))) 
      colnames(xx)=NULL
      rownames(xx)=NULL
      
      f_1=append(f_1, funData(attr(Xs, 'time'), xx ))
      yy=t(as.matrix(as.data.frame(lapply(V, function(x)x[[k]][[1]][2, ]))))
      colnames(yy)=NULL
      rownames(yy)=NULL
      
      f_2=append(f_2, funData(attr(Xs, 'time'), yy ))
    }
    
    X_data=multiFunData(append(f_1, f_2) )
    mean_V=meanFunction(X_data)
    s=norm(mean_V)
    MUs=append(MUs, list(MU) )
    res=c(res, s)
    Align=append(Align, list(X_tilde))
  }
  if(n_random==1){ 
    list(Mu=MU, value=s, align=X_tilde)
  }else{ 
    opt=which.min(res)  
    list(Mu=MUs[[opt]], value=res[opt], align=Align[[opt]])
  }
  
}



# This function is the shapes::frech_mean without the centering part (see package shapes) 
ff_m=function (x, mean = "intrinsic") 
{
  if (mean == "intrinsic") {
    option <- 1
  }
  if (mean == "partial.procrustes") {
    option <- 2
  }
  if (mean == "full.procrustes") {
    option <- 3
  }
  if (mean == "mle") {
    option <- 4
  }
  if (is.double(mean)) {
    if (mean > 0) {
      option <- -mean
    }
  }
  n <- dim(x)[3]
  # for (i in 1:n) {
  #   x[, , i] <- x[, , i]/centroid.size(x[, , i])
  # }
  if (option < 4) {
    pm <- procGPA(x, scale = FALSE, tol1 = 10^(-8))$mshape
    m <- dim(x)[2]
    k <- dim(x)[1]
    ans <- list(mshape = 0, var = 0, code = 0, gradient = 0)
    out <- nlm(objfun, hessian = TRUE, c(pm), uu = x, option = option, 
               iterlim = 1000)
    B <- matrix(out$estimate, k, m)
    ans$mshape <- procOPA(pm, B)$Bhat
    ans$var <- out$minimum
    ans$code <- out$code
    ans$gradient <- out$gradient
  }
  if (option == 4) {
    pm <- procGPA(x, scale = FALSE, tol1 = 10^(-8))$mshape
    m <- dim(x)[2]
    k <- dim(x)[1]
    if (m == 2) {
      theta <- c(log(centroid.size(pm)^2/(4 * 0.1^2)), 
                 pm)
      ans <- list(mshape = 0, kappa = 0, code = 0, gradient = 0)
      out <- nlm(objfun4, hessian = TRUE, theta, uu = x, 
                 iterlim = 1000)
      B <- matrix(out$estimate[-1], k, m)
      ans$mshape <- procOPA(pm, B)$Bhat
      ans$kappa <- exp(out$estimate[1])
      ans$loglike <- -out$minimum
      ans$code <- out$code
      ans$gradient <- out$gradient
    }
    if (m != 2) {
      print("MLE is only appropriate for planar shapes")
    }
  }
  ans
}

## Get the Frechet mean from shapes Xs, using K Fourier coefficient 
frech_mean=function(Xs, K){ 
  bas=mget_alpha(Xs , M=K)
  
  XX=cbind(bas[[1]]$alpha_1, bas[[1]]$alpha_2)
  
  for(k in (2):length(bas)) XX=shapes::abind(XX,cbind(bas[[k]]$alpha_1, bas[[k]]$alpha_2)  )
  
  L=aperm(XX, c( 2, 3,1) )
  
  ## Centering and standardizing 
  C=matrix(apply(colMeans(L[, ,]), mean, MARGIN=1), nrow = 1)[rep(1, dim(L)[1]), ]
  for(k in 1:dim(L)[3]){ 
    L[, , k]=L[, , k]-C; 
    L[, , k]=L[, , k]/norm(L[, , k], 'F')
  }
  
  ff=ff_m(L, mean = "intrinsic")
  
  SS=ff$mshape+C ## Add the center in the estimate 
  SS=SS/norm(SS, 'F')
  tt=attr(Xs, 'time')
  Basis=eval.basis(attr(Xs, 'time'), bas[[1]]$basis)
  X_1=as.shape2(t(cbind(Basis%*%SS[1:(K+1), 1], Basis%*%SS[-(1:(K+1)), 1 ])), tt)
  X_2=as.shape2(t(cbind(Basis%*%SS[1:(K+1), 2], Basis%*%SS[-(1:(K+1)), 2 ])), tt)
  X_3=as.shape2(t(cbind(Basis%*%SS[1:(K+1), 3], Basis%*%SS[-(1:(K+1)), 3 ])), tt)
  as.mvshapes(list(X_1, X_2, X_3))
}
