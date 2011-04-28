\name{mcp.fnc}
\alias{mcp.fnc}
\title{Model criticism plots.}
\description{A function to graph criticism plots for an LMER model (as in Baayen, 2008, chapter 7).}
\usage{mcp.fnc(model, trim = 2.5)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}).}
  \item{trim}{Threshold at which residuals are removed. Defaults to 2.5 (standard deviations above and below the residuals mean).}
}
\details{The first of the four plots the density of the model residuals. The second plot graphs the quantile-quantile plot (actual standardized residuals versus theoretical quantiles). The third plot illustrates the fitted values versus the standardized residuals. The fourth graph plots the absolute of the \code{dffits} of the residuals.}
\value{Returns the four plots described above.}
\references{Baayen, R.H. (2008). \emph{Analyzing Linguistic Data. A Practical Introduction to Statistics Using R}. Cambridge, UK: Cambridge University Press.}
\author{Antoine Tremblay, IWK Health Center, \email{trea26@gmail.com}.}

\examples{
# Example 3
# An example with rcs
library(lme4)
library(Design)
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct=="correct", ]

# base model
m0 = lmer(RT ~  rcs(Frequency, 4) * (Sex + NativeLanguage + 
    FamilySize + SynsetCount + Length) + (1|Subject), 
    data = lexdec3)

# look at model criticism plots
mcp.fnc(m0, trim = 2.5)
}
\keyword{hplot}
