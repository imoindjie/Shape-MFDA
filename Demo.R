setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('function/packages.R')
source('function/object_format.R') ## The definition of the format shapes and mvshapes
source('function/mvshapes_functions.R') ## Associate functions for mvshapes
source('function/align_functions.R') ## Functions for alignment (optimization) 
source('function/Frechet_mean.R') ## Estimation of the Frechet mean 
source('function/basis_expansion.R') ## Basis expansion functions (fda)
source('function/uni_functions.R') ## Functions defined for univariate shape (A functional approach etc. )
source('function/tangent_project.R') ## The transformation into the tangent space 
source('function/classification models.R') ## The used methods for classification 

XY=readRDS('cardio/data.RDS') ## Cardiomegaly data with parameterization and rotation

## Select only a subsample of 50 individuals (For fast results) 

X_t=XY[[1]] 
Y=XY[[2]]

n_0=length(Y)

sub_sample=sample(1:n_0, 50)

X_t=mget_ind(XY[[1]], sub_sample)

## plot(X_t) ## show the curves 
## plot(muntr(X_t)) ## show the preshapes 
Y=XY[[2]][sub_sample]

n=length(Y)

M=10 ### The number of Fourier functions to consider


## The learning dataset

X_apr=mget_ind(X_t, 1:25) ## select the first 25th-individuals as training datasets
Y_apr=Y[1:25] ## 

## The validation set 
X_val=mget_ind(X_t, -(1:25) )
Y_val=Y[-(1:25)]

#### Our approach: use the tangent space for the prediction 
##/!\ Might takes some time ... 
mu=mu_fre(muntr(X_apr), K=M, n_random = 2)$Mu

## Align the shapes with mu
X_val2=multi_al(muntr(X_val), mu, K=M ) ## muntr gives the preshape variables by descarding translation and scaling
X_apr2=multi_al(muntr(X_apr), mu, K=M)

## Figures on the alignment (the estimation of the shape data) 
par(mfrow=c(1, 3))
plot(mu, ylab=expression(X^{(2)}), xlab=expression(X^{(1)}), main='The Fréchet mean', xlim=c(-1, 1), ylim=c(-1, 1))
plot(muntr(X_apr), ylab=expression(X^{(2)}), xlab=expression(X^{(1)}), main='The pre-shapes', ylim=c(-1, 1), xlim=c(-1, 1))
plot(X_apr2, ylab=expression(X^{(2)}), xlab=expression(X^{(1)}), main='The shapes',  ylim=c(-1, 1),xlim=c(-1, 1))
dev.off()

## Transformation of shape to the tangent mapping vectors 
X_apr3=mget_V(mu, X_apr2)
X_val3=mget_V(mu, X_val2)
#### Get the Fourier coefficients
B_exp=mget_alpha(X_apr3, M) ## For the training 

Alpha=NULL
for(kk in 1:length(B_exp)) Alpha=cbind(Alpha, cbind(B_exp[[kk]]$alpha_1,B_exp[[kk]]$alpha_2 ) ) ## Extract the coefficient

B_val=mget_alpha(X_val3, M) ## For the validation dataset

Alpha_val=NULL
for(kk in 1:length(B_val)) Alpha_val=cbind(Alpha_val, cbind(B_val[[kk]]$alpha_1,B_val[[kk]]$alpha_2 ) )

## Results of the models 

Our_approach=c(GPLasso(B_exp, Alpha, Alpha_val, Y_val, Y_apr), LDA(Alpha, Alpha_val, Y_val, Y_apr))

#### Classical approach use the functional covariate without transformation 
##/!\ Might takes some time ... 
B0=mget_alpha(X_apr, M) ## For the training 

A0=NULL
for(kk in 1:length(B0)) A0=cbind(A0, cbind(B0[[kk]]$alpha_1,B0[[kk]]$alpha_2 ) )

B0_val=mget_alpha(X_val, M) ## For the validation dataset

A0_val=NULL
for(kk in 1:length(B0_val)) A0_val=cbind(A0_val, cbind(B0_val[[kk]]$alpha_1,B0_val[[kk]]$alpha_2 ) )

## Results of the models 

Classical=c(GPLasso(B0, A0, A0_val, Y_val, Y_apr), LDA(A0, A0_val, Y_val, Y_apr))

## Summary of the performances (Accuracies in %) of the models 

message(paste0('Our approach: ', paste0((round(Our_approach, 2)*100), collapse = ',')))
message(paste0('Classical FDA: ', paste0((round(Classical, 2)*100), collapse = ',')))

