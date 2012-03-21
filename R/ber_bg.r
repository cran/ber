ber_bg<-function(Y, b, covariates=NULL,stage2=FALSE,partial=TRUE,nSim=150){

################################################################################

if(missing(Y)){stop("Argument 'Y' missing, with no default\n")}

if(missing(b)){stop("Argument 'b' missing, with no default\n")}

if(class(Y)!='matrix'){stop("'Y' must be of class 'matrix'\n")}

if(class(b)!='factor'){stop("'b' must be of class 'factor'\n")}

if(any(is.na(Y))){stop("NA values are not allowed in 'Y'\n")}

if(any(is.na(b))){stop("NA values are not allowed in 'b'\n")}

if(length(b)!=nrow(Y)){stop("length(b) is different from nrow(Y)\n")}

if(any(apply(Y,2,mode)!='numeric')){stop('Array expression columns contain non-numeric values!\n')}

################################################################################

if(!is.null(covariates)){
if(class(covariates)!="data.frame"){stop("'covariates' must be of class 'data.frame'\n")}

col.cov<-ncol(covariates)
for(i in 1:col.cov){if(class(covariates[,i])!="numeric" & class(covariates[,i])!="factor"){
stop("column ", i, " of 'covariates' must be of class 'factor' or 'numeric'\n")}}

if(any(is.na(covariates))){stop("NA values are not allowed in 'covariates'\n")}
}

################################################################################

cnames<-colnames(Y)
rnames<-rownames(Y)

################################################################################
################################################################################
if(is.null(covariates) & partial==TRUE)
{
    library(MASS)
    n <- dim(Y)[1]
    g <- dim(Y)[2]
    m <- nlevels(b)
       
#######################################    
# B bagging estimator     

B <- matrix(rep(0,(m+1)*g), ncol = g, byrow = TRUE)

    for(i in 1:nSim)
    {
    sample_labels <-sample(1:n,n,replace=T)
    Y_bag <- Y[sample_labels,]
    b_bag <- b[sample_labels]
             
    X_bag <- cbind(rep(1, n), model.matrix(~-1 + b_bag, b_bag))
    Xinv_bag <- ginv(X_bag)
    B_bag <- Xinv_bag %*% Y_bag
     
    B<-B+B_bag   
    }
    B<-B/nSim
    
#####################################
# data adjustment

    X <- cbind(rep(1, n), model.matrix(~-1 + b, b))
    Xinv <- ginv(X)
    
    e_1 <- c(1, rep(0, m))                                                      
    e_1_matrix <- matrix(rep(e_1, n), ncol = m + 1, byrow = TRUE)               

    res <- Y - X %*% B
    res_squared <- res^2
    Bdouble <- Xinv %*% res_squared                                             # Bdouble estimator
    genes_var <- Bdouble[1,]
    genes_var_matrix <- matrix(rep(genes_var, m), byrow = TRUE, 
        nrow = m)
    ones <- matrix(rep(1, m * g), nrow = m)
    ScaleFactors <- ones/sqrt(ones + (Bdouble[2:(m + 1),]/genes_var_matrix))
    adjYdouble <- ScaleFactors[b,] * res
    adjYdouble <- adjYdouble + e_1_matrix %*% B
    colnames(adjYdouble)<-cnames
    rownames(adjYdouble)<-rnames
    return(adjYdouble)
}

################################################################################
################################################################################
if(is.null(covariates) & partial==FALSE)
{
    library(MASS)
    n <- dim(Y)[1]
    g <- dim(Y)[2]
    m <- nlevels(b)
       
#######################################    
# bagging estimators     

B <- matrix(rep(0,(m+1)*g), ncol = g, byrow = TRUE)
Bdouble <- matrix(rep(0,(m+1)*g), ncol = g, byrow = TRUE)

    for(i in 1:nSim)
    {
    sample_labels <-sample(1:n,n,replace=T)
    Y_bag <- Y[sample_labels,]
    b_bag <- b[sample_labels]
             
    X_bag <- cbind(rep(1, n), model.matrix(~-1 + b_bag, b_bag))
    Xinv_bag <- ginv(X_bag)
    B_bag <- Xinv_bag %*% Y_bag
    
    adjY_bag <- Y_bag - X_bag %*% B_bag
    res_bag <- Y_bag - X_bag %*% B_bag
    res_squared_bag <- res_bag^2
    Bdouble_bag <- Xinv_bag %*% res_squared_bag
     
    B<-B+B_bag
    Bdouble<-Bdouble+Bdouble_bag   
    }
    B<-B/nSim
    Bdouble<-Bdouble/nSim
    
#####################################
# data adjustment

    X <- cbind(rep(1, n), model.matrix(~-1 + b, b))
    Xinv <- ginv(X)
    
    res <- Y - X %*% B
    
    e_1 <- c(1, rep(0, m))                                                      
    e_1_matrix <- matrix(rep(e_1, n), ncol = m + 1, byrow = TRUE)               

    genes_var <- Bdouble[1,]
    genes_var_matrix <- matrix(rep(genes_var, m), byrow = TRUE, 
        nrow = m)
    ones <- matrix(rep(1, m * g), nrow = m)
    ScaleFactors <- ones/sqrt(ones + (Bdouble[2:(m + 1),]/genes_var_matrix))
    adjYdouble <- ScaleFactors[b,] * res
    adjYdouble <- adjYdouble + e_1_matrix %*% B
    colnames(adjYdouble)<-cnames
    rownames(adjYdouble)<-rnames
    return(adjYdouble)
}


################################################################################
################################################################################
if(!is.null(covariates) & partial==TRUE)
{
    library(MASS)
   
    n <- dim(Y)[1]
    g <- dim(Y)[2]
    m1 <- nlevels(b)
    X1 <- cbind(rep(1, n), model.matrix(~-1 + b, b))
    colnames(covariates) <- paste("col",1:ncol(covariates),sep = "")
    fmla <- as.formula(paste("~", paste(colnames(covariates),collapse = "+")))
    X2 <- model.matrix(fmla, covariates)
    X2 <- as.matrix(X2[,-1])
    X <- cbind(X1, X2)
    m2 <- dim(X2)[2]
    m2double <- m2
    if(stage2==F){m2double <- 0}
    
#####################################
# bagging estimators

    B_Aggregating<-matrix(rep(0,g*(m1+m2+1)),ncol=g,nrow=m1+m2+1)
   
for(i in 1:nSim){
    sample_labels <-sample(1:n,n,replace=T)
    Y_bag <- Y[sample_labels,]
    b_bag <- b[sample_labels]
    covariates_bag <- covariates[sample_labels,]
    covariates_bag <- as.data.frame(covariates_bag)
           
    X1_bag <- cbind(rep(1, n), model.matrix(~-1 + b_bag,b_bag))
    colnames(covariates_bag) <- paste("col",1:ncol(covariates_bag),sep = "")
    fmla <- as.formula(paste("~",paste(colnames(covariates_bag),collapse = "+")))
    X2_bag <- model.matrix(fmla,covariates_bag)
    X2_bag <- as.matrix(X2_bag[,-1])
    m2 <- dim(X2_bag)[2]
    X_bag <- cbind(X1_bag, X2_bag)
    Xinv_bag <- ginv(X_bag)
    B_bag <- Xinv_bag %*% Y_bag
    B1_bag <- B_bag[1:(m1 + 1),]
    B2_bag <- B_bag[(m1 + 2):(m1 + 1 + m2),]
    adjY_bag <- Y_bag - X_bag %*% B_bag
    res_bag <- Y_bag - X_bag %*% B_bag
    res_squared_bag <- res_bag^2

    
   B_Aggregating <- B_Aggregating+B_bag   
   }

   B_Aggregating <- B_Aggregating/nSim
   res_squared <- (Y - X %*% B_Aggregating)^2

   if (stage2 == TRUE) {
        Xinv <-  ginv(X)
        Bdouble_Aggregating <- Xinv %*% res_squared
   }

   if (stage2 == FALSE) {
        X1inv <-  ginv(X1)
        Bdouble_Aggregating <- X1inv %*% res_squared
   }


#####################################
# data adjustment

    X1 <- cbind(rep(1, n), model.matrix(~-1 + b, b))
    colnames(covariates) <- paste("col", 1:ncol(covariates),sep = "")
    fmla <- as.formula(paste("~", paste(colnames(covariates),collapse = "+")))
    X2 <- model.matrix(fmla, covariates)
    X2 <- as.matrix(X2[,-1])
    m2 <- dim(X2)[2]
    X <- cbind(X1, X2)
    Xinv <- ginv(X)
    B <- B_Aggregating
    B1 <- B[1:(m1 + 1),]
    B2 <- B[(m1 + 2):(m1 + 1 + m2),]
    e_1 <- c(1, rep(0, m1 + m2))
    e_1_matrix <- matrix(rep(e_1, n), ncol = m1 + m2 + 1, byrow = TRUE)
    adjY <- Y - X %*% B
    res <- Y - X %*% B
    res_squared <- res^2
    if (stage2 == TRUE) {
        Bdouble <- Bdouble_Aggregating
        Bdouble1 <- Bdouble[1:(m1 + 1),]
        Bdouble2 <- Bdouble[(m1 + 2):(m1 + 1 + m2),]
        genes_var <- Bdouble[1, ]
        genes_var_matrix <- matrix(rep(genes_var, n),byrow = TRUE,nrow = n)
        ones <- matrix(rep(1, n * g), nrow = n)
        delta <- (X %*% Bdouble)/genes_var_matrix
        delta_inv <- ones/delta
        delta2 <- ones + ((X2 %*% Bdouble2)/genes_var_matrix)
        adjYdouble <- delta_inv * adjY
        adjYdouble <- delta2 * adjYdouble
        adjYdouble <- adjYdouble + e_1_matrix %*% B + (X2 %*% B2)
    }
    if (stage2 == FALSE) {
        Bdouble <- Bdouble_Aggregating
        genes_var <- Bdouble[1, ]
        genes_var_matrix <- matrix(rep(genes_var, m1),byrow = TRUE,nrow = m1)
        ones <- matrix(rep(1, m1 * g), nrow = m1)
        ScaleFactors <- ones/sqrt(ones + (Bdouble[2:(m1 + 1),]/genes_var_matrix))
        adjYdouble <- ScaleFactors[b, ] * adjY
        adjYdouble <- adjYdouble + e_1_matrix %*% B + (X2 %*% B2)
    }
    colnames(adjYdouble)<-cnames
    rownames(adjYdouble)<-rnames
    return(adjYdouble)

}

################################################################################
################################################################################
if(!is.null(covariates) & partial==FALSE)
{
    library(MASS)
   
    n <- dim(Y)[1]
    g <- dim(Y)[2]
    m1 <- nlevels(b)
    X1 <- cbind(rep(1, n), model.matrix(~-1 + b, b))
    colnames(covariates) <- paste("col", 1:ncol(covariates), 
        sep = "")
    fmla <- as.formula(paste("~", paste(colnames(covariates), 
        collapse = "+")))
    X2 <- model.matrix(fmla, covariates)
    X2 <- as.matrix(X2[, -1])
    m2 <- dim(X2)[2]
    m2double <- m2
    if(stage2==F){m2double <- 0}
    

    B_Aggregating_bag<-matrix(rep(0,g*(m1+m2+1)),ncol=g,nrow=m1+m2+1)
    Bdouble_Aggregating_bag<-matrix(rep(0,g*(m1+m2double+1)),ncol=g,nrow=m1+m2double+1)

#####################################
# bagging estimators
   
for(i in 1:nSim){
    sample_labels <-sample(1:n,n,replace=T)
    Y_bag <- Y[sample_labels,]
    b_bag <- b[sample_labels]
    covariates_bag <- covariates[sample_labels,]
    covariates_bag <- as.data.frame(covariates_bag)
           
    X1_bag <- cbind(rep(1, n), model.matrix(~-1 + b_bag, b_bag))
    colnames(covariates_bag) <- paste("col",1:ncol(covariates_bag),sep = "")
    fmla <- as.formula(paste("~", paste(colnames(covariates_bag),collapse = "+")))
    X2_bag <- model.matrix(fmla, covariates_bag)
    X2_bag <- as.matrix(X2_bag[,-1])
    m2 <- dim(X2_bag)[2]
    X_bag <- cbind(X1_bag,X2_bag)
    Xinv_bag <- ginv(X_bag)
    B_bag <- Xinv_bag %*% Y_bag
    B1_bag <- B_bag[1:(m1 + 1),]
    B2_bag <- B_bag[(m1 + 2):(m1 + 1 + m2),]
    adjY_bag <- Y_bag - X_bag %*% B_bag
    res_bag <- Y_bag - X_bag %*% B_bag
    res_squared_bag <- res_bag^2

    if (stage2 == TRUE) {
        Bdouble_bag <- Xinv_bag %*% res_squared_bag
    }

    if (stage2 == FALSE) {
        X1inv_bag <- ginv(X1_bag)
        Bdouble_bag <- X1inv_bag %*% res_squared_bag
    }

   B_Aggregating_bag <- B_Aggregating_bag+B_bag 
   Bdouble_Aggregating_bag <- Bdouble_Aggregating_bag+Bdouble_bag  

}

   B_Aggregating <- B_Aggregating_bag/nSim
   Bdouble_Aggregating <- Bdouble_Aggregating_bag/nSim

#####################################
# data adjustment

    X1 <- cbind(rep(1, n), model.matrix(~-1 + b, b))
    colnames(covariates) <- paste("col", 1:ncol(covariates),sep = "")
    fmla <- as.formula(paste("~", paste(colnames(covariates),collapse = "+")))
    X2 <- model.matrix(fmla, covariates)
    X2 <- as.matrix(X2[,-1])
    m2 <- dim(X2)[2]
    X <- cbind(X1, X2)
    Xinv <- ginv(X)
    B <- B_Aggregating
    B1 <- B[1:(m1 + 1),]
    B2 <- B[(m1 + 2):(m1 + 1 + m2),]
    e_1 <- c(1, rep(0, m1 + m2))
    e_1_matrix <- matrix(rep(e_1, n), ncol = m1 + m2 + 1, byrow = TRUE)
    adjY <- Y - X %*% B
    res <- Y - X %*% B
    res_squared <- res^2
    if (stage2 == TRUE) {
        Bdouble <- Bdouble_Aggregating
        Bdouble1 <- Bdouble[1:(m1 + 1),]
        Bdouble2 <- Bdouble[(m1 + 2):(m1 + 1 + m2),]
        genes_var <- Bdouble[1,]
        genes_var_matrix <- matrix(rep(genes_var, n), byrow = TRUE,nrow = n)
        ones <- matrix(rep(1, n * g), nrow = n)
        delta <- (X %*% Bdouble)/genes_var_matrix
        delta_inv <- ones/delta
        delta2 <- ones + ((X2 %*% Bdouble2)/genes_var_matrix)
        adjYdouble <- delta_inv * adjY
        adjYdouble <- delta2 * adjYdouble
        adjYdouble <- adjYdouble + e_1_matrix %*% B + (X2 %*% B2)
    }
    if (stage2 == FALSE) {
        Bdouble <- Bdouble_Aggregating
        genes_var <- Bdouble[1,]
        genes_var_matrix <- matrix(rep(genes_var, m1), byrow = TRUE,nrow = m1)
        ones <- matrix(rep(1, m1 * g), nrow = m1)
        ScaleFactors <- ones/sqrt(ones + (Bdouble[2:(m1 + 1),]/genes_var_matrix))
        adjYdouble <- ScaleFactors[b,] * adjY
        adjYdouble <- adjYdouble + e_1_matrix %*% B + (X2 %*% B2)
    }
    colnames(adjYdouble)<-cnames
    rownames(adjYdouble)<-rnames
    return(adjYdouble)    
}

################################################################################
################################################################################

}
