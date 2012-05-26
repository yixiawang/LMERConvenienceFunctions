\name{romr.fnc}
\alias{romr.fnc}
\title{Exclude outliers.}
\description{Exclude outliers with a standardized residual at a distance greater than 2.5 standard deviations from 0. Note that this function cannot be used with generalized linear mixed-effects models (\code{glmer}s).}
\usage{romr.fnc(model, data, trim = 2.5)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}). Note that this function cannot be used with generalized linear mixed-effects models (\code{glmer}s).}
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
\author{Antoine Tremblay, Dalhousie University, \email{trea26@gmail.com}, with contrbutions from Andy Flies, Michigan State University.}

\seealso{
\code{\link{mcp.fnc}}
\code{\link{perSubjectTrim.fnc}}
}

\examples{
\dontrun{
if("LCFdata" \%in\% .packages(all.available=TRUE)){
  library(LCFdata)
	library(lme4)
	data(eeg)
	eeg <- eeg[eeg$Time >= 100 & eeg$Time <= 175,,drop = TRUE]
	eeg$FreqBd <- as.factor(cut(eeg$FreqB, 3, 
	c("low", "med", "high")))
	eeg$LengthBd <- as.factor(cut(eeg$LengthB, 3, 
	c("short", "med", "long")))
	
	# fit initial model
	m1 = lmer(Fz ~ WMC * (FreqBd + LengthBd) + (1 | Subject) + 
	     (1 | Item), data = eeg)
	
	# look at model criticism plots
	dev.off()
	mcp.fnc(m1, trim = 2.5)
	
	# trim data on initial model residuals
	eeg.trimmed = romr.fnc(m1, eeg, trim = 2.5)
	eeg = eeg.trimmed$data
	
	# update initial model on trimmed data
	mB = update(m1)
	
	# look at model criticism plots
	dev.off()
	mcp.fnc(mB, trim = 2.5)
	
	devAskNewPage(ask=FALSE)
}
}
}
\keyword{manip}
