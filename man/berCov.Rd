\name{berCov}
\alias{berCov}
\title{
Batch Effects Removal using Covariates
}
\description{
The function berCov removes batch effects from (microarray) data. This function uses information about the batches, 
the expression levels in the data and other covariates of interest.
}
\usage{
berCov(Y, b, covariates, stage2_cov = FALSE)
}
\arguments{
  \item{Y}{
A matrix with \eqn{n} rows and \eqn{g} columns, where \eqn{n} is the number of objects and \eqn{g} is the number of variables. 
In the case of gene expression data, columns correspond to genes and rows to samples.
}
  \item{b}{
A vector of class factor with the element in position \eqn{i} (\eqn{i=1,\ldots,n}) representing the batch
from which observation \eqn{i} belongs to.
}
  \item{covariates}{
An object of class data.frame where each column corresponds to a quantitative variable (of class numeric) or 
a qualitative variable (of class factor)
}
  \item{stage2_cov}{
logical indicating if the covariates have to be considered only in handling the location effects (stage2_cov = F) or 
also for the scale effects (stage2_cov = T).
}
}
\value{
A matrix of adjusted data with \eqn{n} rows and \eqn{g} columns.
}
\references{
M. Giordan (Submitted). An efficient two-stage procedure for removing batch
effects in high dimensional data experiments.
}
\author{
Marco Giordan, \email{giordan@stat.unipd.it}
}
\seealso{
\code{\link{ber}}
}
\examples{
library(golubEsets)
library(vsn)
data(Golub_Merge)
E<-exprs(vsn2(Golub_Merge))
batch<-Golub_Merge$Source
BMPB<-Golub_Merge$BM.PB
BMPB<-data.frame(BMPB)
Eadj<-berCov(t(E),batch,BMPB)
}
\keyword{model & regression}

