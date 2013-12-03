\name{ber}
\alias{ber}
\title{Batch Effects Removal}
\description{
Batch effects are removed using a two-stage regression approach.
}
\usage{
ber(Y, b, covariates = NULL)
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
\details{In this implementation NA values are not allowed.}
\value{A matrix of adjusted data with \eqn{n} rows and \eqn{g} columns.}
\references{M. Giordan. February 2013. A Two-Stage Procedure for the Removal of Batch Effects in Microarray Studies. \emph{Statistics in Biosciences}.}
\author{Marco Giordan}
\seealso{\code{\link{ber_bg}}, \code{\link{combat_np}}, \code{\link{combat_p}},\code{\link{mean_centering}}, \code{\link{standardization}}}
\examples{
Y<-matrix(rnorm(6000),nrow=12)
class<-gl(2,6,labels=c("Control","Treat"))
class<-data.frame(class)
batch<-rep(gl(2,3,labels=c("Batch1","Batch2")),2)
YEadj<-ber(Y,batch,class)
}
