ber<-function(Y,b){
# Y = expression matrix with n subjects (rows) and g genes (columns)
# b = vector indicating the batch of each subject; it must be an object of class factor

library(MASS)

################################################################################

n<-dim(Y)[1]   # number of samples
g<-dim(Y)[2]   # number of genes

m <- nlevels(b)  # number of batches              

################################################################################

X <-cbind(rep(1,n),model.matrix(~-1+b,b))

################################################################################

Xinv<-ginv(X)

B<-Xinv%*%Y

e_1<-c(1,rep(0,m))
e_1_matrix<-matrix(rep(e_1,n),ncol=m+1,byrow=TRUE)
adjY<-Y-X%*%B

################################################################################
# residuals

res<-Y-X%*%B
res_squared<-res^2 

################################################################################

Bdouble<-Xinv%*%res_squared

genes_var<-Bdouble[1,]
genes_var_matrix<-matrix(rep(genes_var,m),byrow=TRUE,nrow=m)
ones<-matrix(rep(1,m*g),nrow=m)
ScaleFactors<-ones/sqrt(ones+(Bdouble[2:(m+1),]/genes_var_matrix))

adjYdouble<-ScaleFactors[b,]*adjY
adjYdouble<-adjYdouble + e_1_matrix%*%B

################################################################################

return(adjYdouble)
}
            