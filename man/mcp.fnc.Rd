\name{mcp.fnc}
\alias{mcp.fnc}
\title{
Model criticism plots.
}
\description{
A function to graph criticism plots for an LMER model (as in Baayen, 2008, chapter 7). 
}
\usage{
mcp.fnc(model, data)
}
\arguments{
  \item{model}{
A ``mer'' object (fitted by function ``lmer'').
}
  \item{data}{
The data frame on which the ``mer'' object was fitted.
}
}
\details{
The first of the four plots the density of the model residuals. The second plot graphs the quantile-quantile plot (actual standardized residuals versus theoretical quantiles). The third plot illustrates the fitted values versus the standardized residuals. The fourth graph plots the absolute of the dffits of the residuals.
}
\value{
Returns the four plots described above.
}
\references{
Baayen, R.H. (2008). Analyzing Linguistic Data. A Practical Introduction to Statistics Using R. Cambridge, UK: Cambridge University Press. 
}
\author{
Antoine Tremblay, IWK Health Center, <trea26@gmail.com>
}

\examples{
# Example 3
# An example with rcs
library(Design)
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct=="correct", ]
# base model
m0 = lmer(RT ~  rcs(Frequency, 4) * (Sex + NativeLanguage + 
    FamilySize + SynsetCount + Length) + (1|Subject), 
    data = lexdec3)
# look at model criticism plots
mcp.fnc(m0, lexdec3)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
