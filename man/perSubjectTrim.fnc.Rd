\name{perSubjectTrim.fnc}
\alias{perSubjectTrim.fnc}
\title{Per-subject Trimming of Response Variable.}
\description{For each subject, removes data points that are, e.g., 2.5 standard deviations above or below the subject mean.}
\usage{perSubjectTrim.fnc(data, response, subject, trim = 2.5)}
\arguments{
  \item{data}{The data frame containing the data to be trimmed.}
  \item{response}{The quoted name of the column containing the to-be-trimmed data.}
  \item{subject}{The quoted name of the column contain subject identifiers.}
  \item{trim}{Threshold at which data points will be removed. Defaults to 2.5 (standard deviations above and below each subject's mean).}
}
\value{
The function returns the following objects:
\item{data}{The data with outliers removed.}
\item{data0}{The original data prior to removing the outliers.}
\item{n.removed}{The number of data points removed.}
\item{percent.removed}{The percentage of removed data points.}
}
\author{Antoine Tremblay, Dalhousie University, \email{trea26@gmail.com}.}

\seealso{
\code{\link{mcp.fnc}}
\code{\link{romr.fnc}}
}

\examples{
\dontrun{
if("LCFdata" \%in\% .packages(all.available=TRUE)){
data(eegWide)
dat<-eegWide
rm(eegWide)
gc(TRUE,TRUE)
# per subject trimming
dat <- perSubjectTrim.fnc(dat, response = "Fz", 
      subject = "Subject", trim = 2.5)$data
# ......
# n.removed = 5130 
# percent.removed = 1.584507 
}
}
}
\keyword{manip}

