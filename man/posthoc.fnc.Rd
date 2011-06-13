\name{posthoc.fnc}
\alias{posthoc.fnc}
\title{Posthoc analyses for LMER models.}
\description{The specified factor variable is re-leveled to each one of the factor's levels, the \code{mer} model updated, and, if enabled, MCMC \emph{p}-values returned.}
\usage{posthoc.fnc(model, factor, two.tailed = TRUE, num.comp = NULL,
ndigits = 4, MCMC = FALSE, nsim = 10000, addPlot = FALSE)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}) or an \code{lm} object (fitted by function \code{lm}).}
  \item{factor}{The name of a factor variable (quoted).}
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
\code{\link{romr.fnc}; }
\code{\link{perSubjectTrim.fnc} }
}
\examples{
library(lme4)
library(languageR)
data(dative)
m0 <- glmer(RealizationOfRecipient ~ SemanticClass + 
	LengthOfRecipient + AnimacyOfRec + DefinOfRec + 
	(1 | Speaker) + (1 | Verb), data = dative)
m0.posthoc <- posthoc.fnc(m0, factor = "SemanticClass",
       two.tailed = TRUE, num.comp = 6)
m0.posthoc
}
\keyword{models & regression}
