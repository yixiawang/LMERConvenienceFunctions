\name{bfFixefLMER_F.fnc}
\alias{bfFixefLMER_F.fnc}
\title{Back-fits an LMER model on p-values from ANOVA and, optionally, on LLRT.}
\description{This function back-fits an initial LMER model on upper- or lower-bound \emph{p}-values obtained from function \code{pamer.fnc} and, if enabled, log-likelihood ratio testing.}
\usage{bfFixefLMER_F.fnc(model = as.character(), data = as.character(), 
item = FALSE, alpha = 0.05, llrt = FALSE, p.value = "Upper", 
set.REML.FALSE = TRUE, reset.REML.TRUE = TRUE, log.file = 
paste("fixef_backfit_log_", gsub(":", "_", gsub(" ", "_", date())), 
".txt", sep=""))}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}). Note that logistic models are not supported at the moment.}
  \item{data}{The data frame on which the \code{mer} object was fitted.}
  \item{item}{Whether or not to evaluate the addition of by-item random intercepts to the model, evaluated by way of log-likelihood ratio test. Either \code{FALSE} (the default, does not evaluate this addition) or the column name (quoted) of the item identifier (e.g., \code{"Item"}, or \code{"Word"}).}
  \item{alpha}{Level of significance (for ANOVA and log-likelihood ratio test). Defaults to 0.05.}
  \item{llrt}{Logical. Whether log-likelihood ratio testing between the complex and simple models should be performed. Defaults to \code{FALSE}.}
  \item{p.value}{Whether to use upper-bound (``Upper''; the default) or lower-bound (``Lower'') \emph{p}-values during backfitting.}
  \item{set.REML.FALSE}{Logical. Whether or not to set \code{REML} to \code{FALSE}. Defaults to \code{FALSE}.}
  \item{reset.REML.TRUE}{Logical. Whether or not to re-set the back-fitted model to \code{REML = TRUE}.}
  \item{log.file}{Whether a back-fitting log should be saved. By default, it is saved in a temporary folder (platform dependent) as \code{paste("fixef_backfit_log_", gsub(":", "_", gsub(" ", "_", date())), ".txt", sep = "")}. The path and file name of the log can be changed to whatever the use wishes. Set to \code{FALSE} to disable.}
}
\details{The back-fitting process works as follows:
\enumerate{
\item If argument \code{llrt} is set to \code{TRUE}, \code{REML} is set to \code{FALSE};
\item First consider only highest-order interaction model terms:
  \enumerate{
       \item The model term with the highest ANOVA \emph{p}-value is identified. If this \emph{p}-value is smaller than the alpha value and argument \code{llrt} is set to \code{FALSE}, the model term is removed and a new model is fitted. This is repeated for each model term for which the \emph{p}-value is greater than the alpha value. The algorithm then moves on to step (b). If this \emph{p}-value is greater than the alpha value and \code{llrt} is set to \code{TRUE}, then the following is evaluated:
	\enumerate{
           \item If this \emph{p}-value is greater than the alpha level (set by deault at 0.05), a new model without this model term is fitted; 
	   \item The more complex and simpler models are compared by way of a log-likelihood ratio test (LLRT). If the results of the LLRT determine that the term under consideration does not increase model fit, it is removed; otherwise it is kept.
           \item Move on to the next model term with the largest \emph{p}-value greater than the alpha value and repeat steps (i)--(iii).
	}
        \item Once all highest-order interaction terms have been evaluated, go down to the the second highest order interactions: Repeat steps (ai)--(aiii) with the following addition: If a term would be removed from the model, but it is part of a high-order interaction, keep it. Once all terms of the interaction level have been evaluated, move down to the next lower-order level until main effects have been evaluated, after which the process stops.
  }
\item If argument \code{llrt} is set to \code{TRUE}, set \code{REML} to \code{TRUE} (default) unless otherwise specified.
}

In brief, if \code{llrt} is set to \code{FALSE}, a term remains in the model if its \emph{p}-value from the ANOVA is equal to or smaller than alpha; if \code{llrt} is set to \code{TRUE}, a term remains in the model if
\enumerate{
\item its \emph{p}-value from the ANOVA is equal to or smaller than alpha or, in the case that it is not, if 
\item it significantly increases model fit as determined by a log-likelihood ratio test.
}

If factorial terms are included in the initial model, back-fitting on F is recommended. This backfitting method was used in Newman, Tremblay, Nichols, Neville, and Ullman (submitted).
}
\value{A \code{mer} model with back-fitted fixed effects (on F values) is returned and a log of the back-fitting process is printed on screen and (by default) in a log file in a temporary file (platform dependent).}
\references{
Newman, A.J., Tremblay, A., Nichols, E.S., Neville, H.J., and Ullman, M.T. (submitted). The Influence of Language Proficiency on Lexical-Semantic Processing in Native and Late Learners of English: ERP evidence. Submitted to the \emph{Journal of Cognitive Neuroscience}.
}
\author{Antoine Tremblay, IWK Health Center, \email{trea26@gmail.com}.}

\section{Warnings}{Upper-bound \emph{p}-values can be anti-conservative, while lower-bound \emph{p}-values can be conservative. See \url{http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests} and function \code{pamer.fnc}.}

\seealso{
\code{\link{bfFixefLMER_t.fnc}}
\code{\link{ffRanefLMER.fnc}}
\code{\link{fitLMER.fnc}}
\code{\link{posthoc.fnc}}
\code{\link{pamer.fnc}}
\code{\link{mcp.fnc}}
\code{\link{romr.fnc}}
\code{\link{perSubjectTrim.fnc}}
}
\examples{
# bfFixefLMER -- backfit on F
library(lme4)
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct == "correct", ]

# fit initial model
m1 = lmer(RT ~ Frequency * FamilySize * Length * 
    NativeLanguage + (1|Subject), data = lexdec3)

# look at model criticism plots
mcp.fnc(m1, trim = 2.5)

# trim data on initial model residuals
dat.trimmed = romr.fnc(m1, lexdec3, trim = 2.5)
lexdec3 = dat.trimmed$data

# update initial model on trimmed data
mB = update(m1)

# backfit fixed effects on F
m2 = bfFixefLMER_F.fnc(mB, lexdec3, log.file = FALSE,
     llrt = FALSE)
pamer.fnc(m2)

### An example with rcs
library(Design)
library(lme4)
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct=="correct", ]

# base model
m0 = lmer(RT ~  rcs(Frequency, 4) * (Sex + NativeLanguage + 
    FamilySize + Length) + (1|Subject), data = lexdec3)

# look at model criticism plots
mcp.fnc(m0, trim = 2.5)

# trim data on model residuals
lexdec3 = romr.fnc(m0, lexdec3, trim = 2.5)$data

# update model
m0 = update(m0)

# backfit fixed effects on F
m1 = bfFixefLMER_F.fnc(m0, lexdec3, log.file = FALSE,
	llrt = FALSE)
pamer.fnc(m1)

### An example with poly
library(Design)
library(lme4)
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct=="correct", ]

# base model
m0 = lmer(RT ~  poly(Frequency, 4, raw = TRUE) * (Sex + 
    NativeLanguage + FamilySize + Length) + (1|Subject), 
    data = lexdec3)

# look at model criticism plots
mcp.fnc(m0, trim = 2.5)

# trim data on model residuals
lexdec3 = romr.fnc(m0, lexdec3, trim = 2.5)$data

# update model
m0 = update(m0)

# backfit fixed effects on F
m1 = bfFixefLMER_F.fnc(m0, lexdec3, log.file = FALSE,
	llrt = TRUE)
pamer.fnc(m1)
}
\keyword{models & regression}
