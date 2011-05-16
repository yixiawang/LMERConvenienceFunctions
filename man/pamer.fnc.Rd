\name{pamer.fnc}
\alias{pamer.fnc}
\title{ANOVA with upper- and lower-bound \emph{p}-values and R-sqaured values for LMER.}
\description{Compute upper- and lower-bound \emph{p}-values for the analysis of variance (or deviance) as well as R-squared for each fixed-effect of an LMER model. Note that, at the moment, this function cannot be used with generalized linear mixed-effects models (\code{glmer}s).}
\usage{pamer.fnc(model, ndigits = 4)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}). Note that, at the moment, this function cannot be used with generalized linear mixed-effects models (\code{glmer}s).}
  \item{ndigits}{Integer indicating the number of decimal places to be used in the ANOVA table.}
}
\details{Upper-bound \emph{p}-values are computed by using as denominator \emph{df} \code{nrow(model@frame) - qr(model@X)4rank} (i.e., number of data points minus number of fixed effects including the intercept), which are anti-conservative. Lower-bound \emph{p}-values are computed by using as denominator \emph{df} \code{nrow(model@frame) - qr(model@X)4rank - number of random effects} (e.g., if by-subject intercepts and slopes, and there are 10 subjects, \code{10 * 2 = 20}). See \url{http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests} for more details. R-squared values are calculated as follows: sum of squares of a fixed-effect divided by the sum of squares total, which is calculated here as the variance of the independent variable divided by its length.}
\value{This function returns an object of class \code{data frame} with upper- and lower-bound (anti-conservative and conservative, respectively) \emph{df}s, \emph{p}-values, and R-squared (deviance explained) for each model term.}
\references{\url{http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests}}
\author{Antoine Tremblay, IWK Health Center, \email{trea26@gmail.com}}
\examples{
library(lme4)
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct == "correct", ]
m0 = lmer(RT ~ Sex + NativeLanguage + PrevType + Frequency + 
    FamilySize + SynsetCount + Length + DerivEntropy + (1|Subject), 
    data = lexdec3)
pamer.fnc(m0)
}
\keyword{models & regression}
