\name{posthoc.fnc}
\alias{posthoc.fnc}
\title{Posthoc analyses for LMER models.}
\description{The specified factor variable is re-leveled to each one of the factor's levels, the \code{mer} model updated, and, if enabled, MCMC \emph{p}-values returned.}
\usage{posthoc.fnc(model, data, factor, prior.releveling = NULL, 
two.tailed = TRUE, num.comp = NULL, ndigits = 4, MCMC = FALSE, 
nsim = 10000, addPlot = FALSE)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}) or an \code{lm} object (fitted by function \code{lm}).}
  \item{data}{The unquoted name of the data frame on which the model was fitted.}
  \item{factor}{The name of a factor variable (quoted).}
  \item{prior.releveling}{A list of vectors where the first element of a vector is the quoted name of the factor to be releveled prior to performing the posthoc and the second element of the vector is the quoted name of the reference level, e.g., \code{list(c("AnimacyOfRec", "inanimate"), c("DefinOfRec", "indefinite"))}. Defaults to \code{NULL}.}
  \item{two.tailed}{Logical. Whether to perform one- or two-tailed \emph{t}-tests. Defaults to \code{TRUE}, i.e., two-tailed.}
  \item{num.comp}{An integer indicating the number of comparisons made. This value is used to perform Bonferroni correction of the probability values. Defaults to \code{NULL}.}
  \item{ndigits}{Integer indicating the number of decimal places to be used in the \emph{t} tables. Defaults to 4.}
  \item{MCMC}{Whether to run, for each updated model, function \code{pvals.fnc}. defaults to \code{FALSE}.}
  \item{nsim}{Passed to function \code{pvals.fnc}. An integer denoting the required number of Markov Chain Monte Carlo samples.}
  \item{addPlot}{If \code{TRUE} (default), a trellis graph for the posterior distributions of the parameters is saved in \code{.pdf} file \code{posthoc.fnc_MCMC_plots.pdf} in current working directory.}
}
\value{A list containing, for each factor level, a table of coefficients as returned from function \code{summary} as well as lower- and upper-bound \emph{df}, unadjusted lower- and upper-bound \emph{p}-values, and adjusted lower- and upper-bound \emph{p}-values (Bonferroni corrected). If \code{MCMC = TRUE}, an additional table of coefficients as returned by function \code{pvals.fnc}.}
\author{Antoine Tremblay, IWK Health Center, \email{trea26@gmail.com}.}
\seealso{
\code{\link{bfFixefLMER_F.fnc}; }
\code{\link{bfFixefLMER_t.fnc}; }
\code{\link{ffRanefLMER.fnc}; }
\code{\link{fitLMER.fnc}; }
\code{\link{pamer.fnc}; }
\code{\link{mcp.fnc}; }
\code{\link{mcposthoc.fnc}; }
\code{\link{romr.fnc}; }
\code{\link{perSubjectTrim.fnc} }
}
\examples{
if("LCFdata"\%in\%.packages(all.available=TRUE)){
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
	
	m1.posthoc <- posthoc.fnc(mB, eeg, 
		factor = "FreqBd", prior.releveling = 
		list(c("LengthBd", "long")), two.tailed = TRUE, 
	  num.comp = 6)
	m1.posthoc
	
	devAskNewPage(ask=FALSE)
}
}
\keyword{models & regression}
