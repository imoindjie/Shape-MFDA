packages=c("mclust", "abind", "shapes", 
           "xtable","MFPCA", 
           "caret", "MASS", 
           "pls", "grplasso", 
           "caret", "fda", "cmna")

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

for(pkg in packages) eval(parse(text=paste0('library(', pkg, ')') ))
