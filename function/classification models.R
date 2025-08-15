## For the group lasso methologies for several types of group and different lambdas
## find the best prediction with regard to lambdas  

GPLasso=function(B_exp, Alpha, Alpha_val, Y_val, Y_apr, lambda=27*(0.96)^seq(0, 148, by=5)){ 
  index_1=rep(1:length(B_exp), each=2*(ncol(B_exp[[1]]$alpha_1)) ) ## Each planar curve is a group 
  index_2=rep(1:(2*length(B_exp)), each=ncol(B_exp[[1]]$alpha_1) ) ## Each functional covariate is a group 
  
  GL1 <- grplasso(x=Alpha, y = as.vector(Y_apr*1), index = index_1, lambda = lambda, model = LogReg(),
                  penscale = sqrt,
                  control = grpl.control(update.hess = "lambda", trace = 0), center=F)
  
  GL2 <- grplasso(x=Alpha, y = as.vector(Y_apr*1), index = index_2, lambda = lambda, model = LogReg(),
                  penscale = sqrt,
                  control = grpl.control(update.hess = "lambda", trace = 0), center=F)
  
  ## Prediction 
  
  GLPred1=predict(GL1, Alpha_val)
  GLPred2=predict(GL2, Alpha_val)
  
  ## Get the scores of all the prediction 
  
  GLRes1=unlist(lapply(1:ncol(GLPred1), function(x)sum(diag(table((GLPred1[, x]>0.5), Y_val)))/length(Y_val) ))
  GLRes2=unlist(lapply(1:ncol(GLPred2), function(x)sum(diag(table((GLPred2[, x]>0.5), Y_val)))/length(Y_val) )) 
  res=c(max(GLRes1), max(GLRes2))## Return the scores of each model 
  names(res)=c('GL1', 'GL2')
  res
  }

## For linear discriminant analysis based (PLS and PCR)
## Test all components and give the best prediction 

LDA=function(Alpha, Alpha_val, Y_val, Y_apr){ 
  ## PLS methodology testing all the number of components  
  PLS=NULL
  for(k in 1:(min(c(nrow(Alpha), ncol(Alpha)))-1) ){ 
    prec_k=NULL
    x1=plsr(ifelse(Y_apr, -1, 1)~Alpha, data=as.data.frame(Alpha), ncomp=k)
    p1=predict(x1, Alpha_val)<0
    prec_k=c(prec_k, mean(p1==Y_val)) ## get the score of prediction for k
    PLS=c(PLS, mean(prec_k) )
  }
  ## PCR methodology testing all the number of components  
  PCR=NULL
  for(k in 1:(min(c(nrow(Alpha), ncol(Alpha)))-1) ){ 
    prec_k=NULL
    x1=pcr(ifelse(Y_apr, -1, 1)~Alpha, data=as.data.frame(Alpha), ncomp=k)
    p1=predict(x1, Alpha_val)<0
    prec_k=c(prec_k, mean(p1==Y_val)) 
    PCR=c(PCR, mean(prec_k) )## get the score of prediction for k
  }
  ## Return the best scores 
  res=c(max(PLS), max(PCR))
  names(res)=c('PLS', 'PCR')
  res  
  }
