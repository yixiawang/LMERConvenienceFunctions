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
\author{Antoine Tremblay, IWK Health Center, \email{trea26@gmail.com}.}

\seealso{
\code{\link{mcp.fnc}}
\code{\link{romr.fnc}}
}

\examples{
if("LCFdata"\%in\%.packages(all.available=TRUE)){
  library(LCFdata)
	data(eeg)
	eeg <- eeg[eeg$Time >= 100 & eeg$Time <= 175, , drop = TRUE]
	
	eeg <- perSubjectTrim.fnc(eeg, "Fz", "Subject", 2.5)
	eeg <- eeg$data
	
	# fit initial model
	m1 = lmer(Fz ~ FreqB * LengthB * WMC + (1 | Subject) + 
	     (1 | Item), data = eeg)
	
	# look at model criticism plots
	mcp.fnc(m1, trim = 2.5)
	
	devAskNewPage(ask=FALSE)
}
}
\keyword{manip}

