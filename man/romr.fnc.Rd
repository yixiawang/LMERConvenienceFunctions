\name{romr.fnc}
\alias{romr.fnc}
\title{Exclude outliers.}
\description{Exclude outliers with a standardized residual at a distance greater than 2.5 standard deviations from 0.}
\usage{romr.fnc(model, data, trim = 2.5)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}).}
  \item{data}{The data frame on which the \code{mer} object was fitted.}
  \item{trim}{Threshold at which residuals will be removed. Defaults to 2.5 (standard deviations above and below the residuals mean).}
}
\value{
The function returns the following objects:
\item{data}{The data with outliers removed.}
\item{data0}{The original data prior to removing the outliers.}
\item{n.removed}{The number of data points removed.}
\item{percent.removed}{The percentage of removed data points.}
}
\references{
Baayen, R.H. (2008). \emph{Analyzing Linguistic Data. A Practical Introduction to Statistics Using R}. Cambridge, UK: Cambridge University Press. 

Newman, A.J., Tremblay, A., Nichols, E.S., Neville, H.J., and Ullman, M.T. (submitted). The Influence of Language Proficiency on Lexical-Semantic Processing in Native and Late Learners of English: ERP evidence. Submitted to the \emph{Journal of Cognitive Neuroscience}.

Tremblay, A. and Tucker B. V. (submitted). What can the production of four-word sequences tell us about the mental lexicon? Submitted to \emph{The Mental Lexicon}.
}
\author{Antoine Tremblay, IWK Health Center, \email{trea26@gmail.com}.}

\seealso{
\code{\link{perSubjectTrim.fnc}}
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
dat.trimmed = romr.fnc(m0, lexdec3, trim = 2.5)
names(dat.trimmed)
lexdec3 = dat.trimmed$data
m0 = update(m0)
mcp.fnc(m0, trim = 2.5)
}
\keyword{manip}
