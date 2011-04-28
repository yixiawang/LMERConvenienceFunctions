\name{pamer.fnc}
\alias{pamer.fnc}
\title{ANOVA with upper- and lower-bound \emph{p}-values for LMER models.}
\description{Compute upper- and lower-bound \emph{p}-values for the analysis of variance (or deviance) of an LMER model.}
\usage{pamer.fnc(model, ndigits = 4)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}).}
  \item{ndigits}{Integer indicating the number of decimal places to be used in the ANOVA table.}
}
\details{Upper-bound \emph{p}-values are computed by using as denominator \emph{df} \code{nrow(model@frame) - qr(model@X)4rank} (i.e., number of data points minus number of fixed effects including the intercept), which are anti-conservative. Lower-bound \emph{p}-values are computed by using as denominator \emph{df} \code{nrow(model@frame) - qr(model@X)4rank - number of random effects} (e.g., if by-subject intercepts and slopes, and there are 10 subjects, \code{10 * 2 = 20}). See \url{http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests} for more details.}
\value{This function returns an object of class \code{data frame} with upper- and lower-bound (anti-conservative and conservative, respectively) \emph{df}s and \emph{p}-values for each model term.}
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
