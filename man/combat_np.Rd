\name{combat_np}
\alias{combat_np}
\title{non-parametric combat}
\description{
Batch effects are removed using a non-parametric empirical Bayes approach, as described in Johnson et al. 2007.
}
\usage{
combat_np(Y, b, covariates = NULL)
}
\arguments{
\item{Y}{A matrix with \eqn{n} rows and \eqn{g} columns, where \eqn{n} is the number of objects and \eqn{g}
is the number of variables. In the case of gene expression data, columns correspond
to genes (probe sets) and rows to samples.}
\item{b}{A vector of class factor with the element in position \eqn{i} (\eqn{i=1,\ldots,n}) representing
the batch from which observation \eqn{i} belongs to.}
\item{covariates}{An object of class data.frame where each column corresponds to a quantitative
variable (of class numeric) or a qualitative variable (of class factor).}
}
\details{In this implementation pre-processing of the data is obtained through Moore-Penrose pseudo inverse.  
Therefore the final matrix of adjusted data can be different from that obtained through the 
original ComBat algorithm (see sva package). NA values are not allowed.}
\value{A matrix of adjusted data with \eqn{n} rows and \eqn{g} columns.}
\references{Johnson et al. 2007 \url{http://www.ncbi.nlm.nih.gov/pubmed/16632515}}
\author{Marco Giordan}
\seealso{\code{\link{ber}}, \code{\link{ber_bg}}, \code{\link{combat_p}}, \code{\link{mean_centering}}, \code{\link{standardization}}}
\examples{
library(golubEsets)
library(vsn)
data(Golub_Merge)
E<-exprs(vsn2(Golub_Merge))
E<-E[1:100,]
batch<-Golub_Merge$Source
BMPB<-Golub_Merge$BM.PB
BMPB<-data.frame(BMPB)
Eadj<-combat_np(t(E),batch,BMPB)
}
