\name{ber_bg}
\alias{ber_bg}
\title{Batch Effects Removal using Bagging}
\description{
Batch effects are removed using a two-stage regression approach and bagging.
}
\usage{
ber_bg(Y, b, covariates = NULL, stage2=FALSE,partial=TRUE,nSim=150)
}
\arguments{
\item{Y}{A matrix with \eqn{n} rows and \eqn{g} columns, where \eqn{n} is the number of objects and \eqn{g}
is the number of variables. In the case of gene expression data, columns correspond
to genes (probe sets) and rows to samples.}
\item{b}{A vector of class factor with the element in position \eqn{i} (\eqn{i=1,\ldots,n}) representing
the batch from which observation \eqn{i} belongs to.}
\item{covariates}{An object of class data.frame where each column corresponds to a quantitative
variable (of class numeric) or a qualitative variable (of class factor).}
\item{stage2}{A logical value indicating if the covariates have to be considered only in handling the
location effects (stage2 = F) or also for the scale effects (stage2 = T).}
\item{partial}{A logical value indicating if partial bagging or full bugging have to be performed. See reference below.}
\item{nSim}{Number of bootstrap samples.}
}
\details{In this implementation NA values are not allowed.}
\value{A matrix of adjusted data with \eqn{n} rows and \eqn{g} columns.}
\references{M. Giordan. A two-stage procedure for the removal of batch effects in microarray studies}
\author{Marco Giordan}
\seealso{\code{\link{ber}}, \code{\link{combat_np}}, \code{\link{combat_p}},  \code{\link{mean_centering}}, \code{\link{standardization}}}
\examples{
library(golubEsets)
library(vsn)
data(Golub_Merge)
E<-exprs(vsn2(Golub_Merge))
batch<-Golub_Merge$Source
BMPB<-Golub_Merge$BM.PB
BMPB<-data.frame(BMPB)
Eadj<-ber_bg(t(E),batch,BMPB)
}
