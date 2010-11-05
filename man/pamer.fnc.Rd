\name{pamer.fnc}
\alias{pamer.fnc}
\title{
ANOVA with upper-bound p values for LMER model.
}
\description{
Compute upper-bound p values for the analysis of variance (or deviance) of an LMER model.
}
\usage{
pamer.fnc(model, ndigits = 4)
}
\arguments{
  \item{model}{
A ``mer'' object (fitted by function ``lmer'').
}
  \item{ndigits}{
Integer indicating the number of decimal places (`round') or significant digits (`signif') to be used.
}
}
\details{
Probability values are computed by using as denominator df ``nrow(model@frame) - qr(model@X)4rank'', which are anti-conservative. See \url{http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests} for more details. 
}
\value{
This function returns an object of class ``data frame'' with upper-bound (anti-conservative) p values indicating which terms are significant.
}
\references{
\url{http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests}
}
\author{
Antoine Tremblay, IWK Health Center, <trea26@gmail.com>
}

\examples{
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct == "correct", ]
m0 = lmer(RT ~ Sex + NativeLanguage + PrevType + Frequency + FamilySize + SynsetCount + Length + DerivEntropy + (1|Subject), data = lexdec3)
pamer.fnc(m0)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
