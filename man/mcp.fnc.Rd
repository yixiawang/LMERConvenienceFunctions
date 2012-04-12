\name{mcp.fnc}
\alias{mcp.fnc}
\title{Model criticism plots.}
\description{A function to graph criticism plots for an LMER model (as in Baayen, 2008, chapter 7). Note that this function cannot be used with generalized linear mixed-effects models (GLMERs).}
\usage{mcp.fnc(model, trim = 2.5, col = "red")}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}). Note that, at the moment, this function cannot be used with generalized linear mixed-effects models (GLMERs).}
  \item{trim}{Used to plot lines in the fitted ~ standardized residuals plot. The lines correspond to the threshold at which residuals would be or were removed. Defaults to 2.5 (standard deviations above and below the residuals mean).}
  \item{col}{Color of the liines added to the quantile-quantile plot and fitted ~ standardized residuals plot. Defaults to red.}
}
\details{The first of the four plots graphs the density of the model residuals. The second plot graphs the quantile-quantile plot (actual standardized residuals versus theoretical quantiles). The third plot illustrates the fitted values versus the standardized residuals. The fourth graph plots the absolute values of the \code{dffits} of the residuals.}
\value{Returns the four plots described above.}
\references{Baayen, R.H. (2008). \emph{Analyzing Linguistic Data. A Practical Introduction to Statistics Using R}. Cambridge, UK: Cambridge University Press.}
\author{Antoine Tremblay, Dalhousie University, \email{trea26@gmail.com}.}

\examples{
\dontrun{
if("LCFdata"\%in\%.packages(all.available=TRUE)){
  library(LCFdata)
	library(lme4)
	data(eeg)
	eeg <- eeg[eeg$Time >= 100 & eeg$Time <= 175, , drop = TRUE]
	
	# fit initial model
	m1 = lmer(Fz ~ FreqB * LengthB * WMC + (1 | Subject) + 
	     (1 | Item), data = eeg)
	
	# look at model criticism plots
	mcp.fnc(m1, trim = 2.5)
	
	devAskNewPage(ask=FALSE)
}
}
}
\keyword{hplot}
