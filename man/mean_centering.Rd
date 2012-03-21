\name{mean_centering}
\alias{mean_centering}
\title{mean_centering}
\description{
Batch effects are removed using the means of the batches.
}
\usage{
mean_centering(Y, b)
}
\arguments{
\item{Y}{A matrix with \eqn{n} rows and \eqn{g} columns, where \eqn{n} is the number of objects and \eqn{g}
is the number of variables. In the case of gene expression data, columns correspond
to genes (probe sets) and rows to samples.}
\item{b}{A vector of class factor with the element in position \eqn{i} (\eqn{i=1,\ldots,n}) representing
the batch from which observation \eqn{i} belongs to.}
}
\value{A matrix of adjusted data with \eqn{n} rows and \eqn{g} columns.}
\author{Marco Giordan}
\seealso{\code{\link{ber}}, \code{\link{ber_bg}}, \code{\link{combat_np}}, \code{\link{combat_p}},  \code{\link{standardization}}}
\examples{
library(golubEsets)
library(vsn)
data(Golub_Merge)
E<-exprs(vsn2(Golub_Merge))
batch<-Golub_Merge$Source
Eadj<-mean_centering(t(E),batch)
}
