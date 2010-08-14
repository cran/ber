\name{ber}
\alias{ber}
\title{
Batch Effects Removal
}
\description{
The function ber removes the batch effects from a data set.
}
\usage{
ber(Y, b)
}
\arguments{
  \item{Y}{A matrix with \eqn{n} rows and \eqn{g} columns, where \eqn{n} is the number of objects and \eqn{g} is the number of variables. 
In the case of gene expression data, columns correspond to genes and rows to samples.
}
  \item{b}{A vector of integers of length \eqn{n} where the elements correspond to the batches. For \eqn{m} batches 
the elements must be integers between 1 and \eqn{m}; the element in position \eqn{i} (\eqn{i=1,\ldots,n}) represents the batch
from which observation \eqn{i} belongs to.
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
\examples{
library(golubEsets)
library(vsn)
data(Golub_Merge)
b<-Golub_Merge$Source

E<-exprs(vsn2(Golub_Merge))
batch<-rep(1,length(b))
for(i in 1:length(b)){
    if(b[i]=="CALGB"){batch[i]<-2}
    if(b[i]=="CCG"){batch[i]<-3}
    if(b[i]=="St-Jude"){batch[i]<-4}
}

Eadj<-ber(t(E),batch)
}
\keyword{model & regression}
