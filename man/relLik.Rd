\name{relLik}
\alias{relLik}
\title{Relative log-likehood.}
\description{Calculate the relative log-likehood between two models.}
\usage{relLik(x, y, method = c("AIC", "BIC"), ndigits = 6, ...)}
\arguments{
  \item{x,y}{Fitted model objects for which there exists a \code{logLik} method to extract the corresponding log-likelihood, or objects inheriting from class \code{logLik}.}
  \item{method}{Whether to base the comparison on \code{AIC} or \code{BIC}. Defaults to "AIC".}
  \item{ndigits}{An integer denoting the number of decimal digits in the output.}
  \item{\dots}{Further arguments to pass to \code{AIC} or \code{BIC}.}
}
\details{The relative log-likelihood is calculated as \code{exp((abs(AIC(x) - AIC(y)))/2)} or \code{exp((abs(BIC(x) - BIC(y)))/2)}, depending on the method.}
\value{A vector with values:
\item{AIC(x), BIC(x) }{The \code{AIC} or ]code{BIC} value of the first model object.}
\item{AIC(y), BIC(y) }{The \code{AIC} or ]code{BIC} value of the second model object.}
\item{relLik }{The relative likelihood between the two models. Model \code{y} will be that much more likely given the data than model \code{x}.}
}
\author{Antoine Tremblay, Dalhousie University, \email{trea26@gmail.com}}
\seealso{
\code{\link{logLik}; }
\code{\link{AIC}; }
\code{\link{BIC} }
}
\examples{
# see example in LMERConvenienceFunctions help page.
}
