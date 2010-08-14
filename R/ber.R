ber<-function(Y,b){
# Y = expression matrix with n subjects (rows) and g genes (columns)
# b = vector indicating the batch of each subject; they must be indicated with the number 1,2,...,m

library(MASS)

################################################################################

n<-dim(Y)[1]   # number of samples
g<-dim(Y)[2]   # number of genes

m<-length(unique(b))  # number of batches

################################################################################
# removal of batch-location effects

X<-matrix(rep(0,n*(m+1)),nrow=n)
X[,1]<-rep(1,n)
for(i in 1:n){
for(j in 1:m){
if(b[i]==j)X[i,j+1]<-1
}
}

Xinv<-ginv(X)

B<-Xinv%*%Y

e_1<-c(1,rep(0,m))
e_1_matrix<-matrix(rep(e_1,n),ncol=m+1,byrow=TRUE)
adjY<-Y-X%*%B+e_1_matrix%*%B

A<-matrix(rep(1,n*n),ncol=n)
means_matrix<-(1/n)*A%*%adjY

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

adjYdouble<-adjY-means_matrix
adjYdouble<-ScaleFactors[b,]*adjYdouble # this line replaces this :
#for(i in 1:n){
#adjYdouble[i,]<-ScaleFactors[b[i],]*adjYdouble[i,]
#}  
adjYdouble<-adjYdouble+means_matrix

################################################################################

return(adjY)
}
            