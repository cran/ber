berCov<-function (Y, b, covariates,stage2_cov=FALSE){
# Y = expression matrix with n subjects (rows) and g genes (columns)
# b = vector indicating the batch of each subject; it must be an object of class factor
# covariates = object of class data.frame where each column represents a quantitative 
# (object of class numeric) or qualitative variable (object of class factor)

library(MASS)
    
################################################################################
    
n <- dim(Y)[1]    # number of samples
g <- dim(Y)[2]    # number of genes
m1 <- nlevels(b)    # number of batches
    
################################################################################                
    
X1 <-cbind(rep(1,n),model.matrix(~-1+b,b))             
    
colnames(covariates)<-paste("col",1:ncol(covariates),sep="")
fmla<-as.formula(paste("~",paste(colnames(covariates),collapse="+")))
X2<-model.matrix(fmla,covariates)
X2<-as.matrix(X2[,-1])
m2<-dim(X2)[2]
    
X<-cbind(X1,X2)
    
################################################################################
    
Xinv <- ginv(X)
B <- Xinv %*% Y
    
B1<-B[1:(m1+1),]
B2<-B[(m1+2):(m1+1+m2),]
    
e_1 <- c(1, rep(0, m1 + m2))
e_1_matrix <- matrix(rep(e_1, n), ncol = m1 + m2 + 1, byrow = TRUE)
adjY <- Y - X %*% B 
        
################################################################################
    
res <- Y - X %*% B
res_squared <- res^2
    
################################################################################
    
if(stage2_cov == TRUE){
Bdouble <- Xinv %*% res_squared
Bdouble1 <- Bdouble[1:(m1+1),]
Bdouble2 <- Bdouble[(m1+2):(m1+1+m2),]
    
genes_var <- Bdouble[1,]
genes_var_matrix <- matrix(rep(genes_var, n), byrow = TRUE,
nrow = n)
ones <- matrix(rep(1, n * g), nrow = n)
    
delta <- (X %*% Bdouble)/genes_var_matrix
delta_inv <- ones/delta
    
delta2 <- ones + ((X2 %*% Bdouble2)/genes_var_matrix)
        
adjYdouble <- delta_inv * adjY
adjYdouble <- delta2 * adjYdouble
adjYdouble <- adjYdouble + e_1_matrix %*% B + (X2 %*% B2)
}
    
if(stage2_cov == FALSE){
X1inv <- ginv(X1)
Bdouble <- X1inv %*% res_squared
genes_var <- Bdouble[1,]
genes_var_matrix <- matrix(rep(genes_var, m1), byrow = TRUE, 
nrow = m1)
ones <- matrix(rep(1, m1 * g), nrow = m1)
ScaleFactors <- ones/sqrt(ones + (Bdouble[2:(m1+1),]/genes_var_matrix))
adjYdouble <- ScaleFactors[b,] * adjY
adjYdouble <- adjYdouble + e_1_matrix %*% B + (X2 %*% B2)
}

################################################################################
    
return(adjYdouble)
}
