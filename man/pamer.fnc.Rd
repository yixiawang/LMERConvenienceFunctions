\name{pamer.fnc}
\alias{pamer.fnc}
\title{ANOVA with upper- and lower-bound \emph{p}-values and R-sqaured values for LMER.}
\description{Compute upper- and lower-bound \emph{p}-values for the analysis of variance (or deviance) as well as the amount of deviance explained (\%) for each fixed-effect of an LMER model. Note that for \code{glmer} models, there is no deviance explained column.}
\usage{pamer.fnc(model, ndigits = 4)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}). Note that, at the moment, this function cannot be used with generalized linear mixed-effects models (\code{glmer}s).}
  \item{ndigits}{Integer indicating the number of decimal places to be used in the ANOVA table.}
}
\details{Upper-bound \emph{p}-values are computed by using as denominator \emph{df} \code{nrow(model@frame) - qr(model@X)4rank} (i.e., number of data points minus number of fixed effects including the intercept), which are anti-conservative. Lower-bound \emph{p}-values are computed by using as denominator \emph{df} \code{nrow(model@frame) - qr(model@X)4rank - number of random effects} (e.g., if by-subject intercepts and slopes, and there are 10 subjects, \code{10 * 2 = 20}). The amount of deviance explained by each model term (i.e., eta squared) is calculated as \code{[Sum of Squares for the effect] / [Sum of Squares total]}. More specifically: \code{as.data.frame(anova(model))[,2] / sum((model@frame[, dv]-mean(model@frame[, dv]))^2)} where \code{dv} is a vector of the names of the independent variables in the model.} 
\value{This function returns an object of class \code{data frame} with upper- and lower-bound (anti-conservative and conservative, respectively) \emph{df}s, \emph{p}-values, and deviance explained (\%) for each model term. Note that for \code{glmer} models, there is no deviance explained column.}
\author{Antoine Tremblay, Statistics Canada, \email{trea26@gmail.com}}
\references{\code{[R] lmer, p-values and all that} available at https://stat.ethz.ch/pipermail/r-help/2006-May/094765.html.}
\examples{
# see example LMERConvenienceFunctions help page.
}
\keyword{models & regression}
