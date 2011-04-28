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
\code{\link{romr.fnc}}
}

\examples{
library(lme4)
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct == "correct", ]
m0 = lmer(RT ~ Sex + NativeLanguage + PrevType + Frequency + 
    FamilySize + SynsetCount + Length + DerivEntropy + (1|Subject), 
    data = lexdec3)
mcp.fnc(m0, trim = 2.5)
dat.trimmed = perSubjectTrim.fnc(data = lexdec3, subject = "Subject",
	response = "RT", trim = 2.5)
lexdec3 = dat.trimmed$data
m0 = update(m0)
mcp.fnc(m0, trim = 2.5)
}
\keyword{manip}

