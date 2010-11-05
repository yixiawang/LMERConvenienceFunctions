\name{romr.fnc}
\alias{romr.fnc}
\title{
Exclude outliers.
}
\description{
Exclude outliers with a standardized residual at a distance greater than 2.5 standard deviations from 0.
}
\usage{
romr.fnc(model, data)
}
\arguments{
  \item{model}{
A ``mer'' object (fitted by function ``lmer'').
}
  \item{data}{
The data frame on which the ``mer'' object was fitted.
}
}
\value{
The function returns the following objects:
\item{data}{The data with outliers removed.}
\item{data0}{The original data prior to removing the outliers.}
\item{n.removed}{The number of data points removed.}
\item{percent.removed}{The percentage of removed data points.}
}
\references{
Baayen, R.H. (2008). Analyzing Linguistic Data. A Practical Introduction to Statistics Using R. Cambridge, UK: Cambridge University Press. 

Newman, A.J., Tremblay, A., Nichols, E.S., Neville, H.J., and Ullman, M.T. (submitted). The Influence of Language Proficiency on Lexical-Semantic Processing in Native and Late Learners of English: ERP evidence. Submitted to the Journal of Cognitive Neuroscience.

Tremblay, A. and Tucker B. V. (in preparation). What can the production of four-word sequences tell us about the mental lexicon? To be submitted to The Mental Lexicon.
}
\author{
Antoine Tremblay, IWK Health Center, <trea26@gmail.com>
}

\examples{
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct == "correct", ]
m0 = lmer(RT ~ Sex + NativeLanguage + PrevType + Frequency + FamilySize + SynsetCount + Length + DerivEntropy + (1|Subject), data = lexdec3)
mcp.fnc(m0, lexdec3)
dat.trimmed = romr.fnc(m0, lexdec3)
names(dat.trimmed)
lexdec3 = dat.trimmed$data
m0 = update(m0)
mcp.fnc(m0, lexdec3)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
