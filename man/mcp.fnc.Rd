\name{mcp.fnc}
\alias{mcp.fnc}
\title{Model criticism plots.}
\description{A function to graph criticism plots for an LMER model (as in Baayen, 2008, chapter 7). Note that this function cannot be used with generalized linear mixed-effects models (GLMERs). Also note that the fourth plot (dffits) is omitted until we can figure out how to calculate dffits for a \code{merMod} object.}
\usage{mcp.fnc(model, trim = 2.5, col = "red")}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}). Note that, at the moment, this function cannot be used with generalized linear mixed-effects models (GLMERs).}
  \item{trim}{Used to plot lines in the fitted ~ standardized residuals plot. The lines correspond to the threshold at which residuals would be or were removed. Defaults to 2.5 (standard deviations above and below the residuals mean).}
  \item{col}{Color of the lines added to the quantile-quantile plot and fitted ~ standardized residuals plot. Defaults to red.}
}
\details{The first of the four plots graphs the density of the model residuals. The second plot graphs the quantile-quantile plot (actual standardized residuals versus theoretical quantiles). The third plot illustrates the fitted values versus the standardized residuals. The fourth graph plots the absolute values of the \code{dffits} of the residuals (not producing this plot as of version 2.2; might come back in future versions).}
\value{Returns the four plots described above.}
\references{Baayen, R.H. (2008). \emph{Analyzing Linguistic Data. A Practical Introduction to Statistics Using R}. Cambridge, UK: Cambridge University Press.}
\author{Antoine Tremblay, Statistics Canada, \email{trea26@gmail.com}.}

\examples{
# see example LMERConvenienceFunctions help page.
}
\keyword{hplot}
