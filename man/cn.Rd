\name{cn}
\alias{cn}
\title{
List the column names of a data frame in matrix format.
}
\description{
The colum names of the specified data frame are listed in matrix format, that is, each one appears in one row preceded by the row number.
}
\usage{
cn(data.frame)
}
\arguments{
  \item{data.frame}{
A data frame.
}
}
\value{
A matrix containing the data frame column names.
}
\author{
Antoine Tremblay, IWK Health Center, <trea26@gmail.com>
}
\seealso{
\code{\link{colnames}}
}

\examples{
data(lexdec)
cn(lexdec)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
