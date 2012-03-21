ber<-function(Y, b, covariates=NULL, stage2 = FALSE){

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

if(is.null(covariates)){
    library(MASS)
    n <- dim(Y)[1]
    g <- dim(Y)[2]
    m <- nlevels(b)
    X <- cbind(rep(1, n), model.matrix(~-1 + b, b))
    Xinv <- ginv(X)
    B <- Xinv %*% Y
    e_1 <- c(1, rep(0, m))
    e_1_matrix <- matrix(rep(e_1, n), ncol = m + 1, byrow = TRUE)
    adjY <- Y - X %*% B
    res <- Y - X %*% B
    res_squared <- res^2
    Bdouble <- Xinv %*% res_squared
    genes_var <- Bdouble[1,]
    genes_var_matrix <- matrix(rep(genes_var, m), byrow = TRUE,nrow = m)
    ones <- matrix(rep(1, m * g), nrow = m)
    ScaleFactors <- ones/sqrt(ones + (Bdouble[2:(m + 1),]/genes_var_matrix))
    adjYdouble <- ScaleFactors[b,] * adjY
    adjYdouble <- adjYdouble + e_1_matrix %*% B
    colnames(adjYdouble)<-cnames
    rownames(adjYdouble)<-rnames
    return(adjYdouble)
    }else{
    library(MASS)
    n <- dim(Y)[1]
    g <- dim(Y)[2]
    m1 <- nlevels(b)
    X1 <- cbind(rep(1, n), model.matrix(~-1 + b, b))
    colnames(covariates) <- paste("col", 1:ncol(covariates),sep = "")
    fmla <- as.formula(paste("~", paste(colnames(covariates),collapse = "+")))
    X2 <- model.matrix(fmla, covariates)
    X2 <- as.matrix(X2[,-1])
    m2 <- dim(X2)[2]
    X <- cbind(X1, X2)
    Xinv <- ginv(X)
    B <- Xinv %*% Y
    B1 <- B[1:(m1 + 1),]
    B2 <- B[(m1 + 2):(m1 + 1 + m2),]
    e_1 <- c(1, rep(0, m1 + m2))
    e_1_matrix <- matrix(rep(e_1, n), ncol = m1 + m2 + 1, byrow = TRUE)
    adjY <- Y - X %*% B
    res <- Y - X %*% B
    res_squared <- res^2
    if (stage2 == TRUE){
        Bdouble <- Xinv %*% res_squared
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
    if (stage2 == FALSE){
        X1inv <- ginv(X1)
        Bdouble <- X1inv %*% res_squared
        genes_var <- Bdouble[1, ]
        genes_var_matrix <- matrix(rep(genes_var, m1), byrow = TRUE, nrow = m1)
        ones <- matrix(rep(1, m1 * g), nrow = m1)
        ScaleFactors <- ones/sqrt(ones + (Bdouble[2:(m1 + 1),]/genes_var_matrix))
        adjYdouble <- ScaleFactors[b,] * adjY
        adjYdouble <- adjYdouble + e_1_matrix %*% B + (X2 %*% B2)
    }
    colnames(adjYdouble)<-cnames
    rownames(adjYdouble)<-rnames
    return(adjYdouble)
    }
}
