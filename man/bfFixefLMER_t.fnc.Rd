\name{bfFixefLMER_t.fnc}
\alias{bfFixefLMER_t.fnc}
\title{
Back-fits an LMER model on absolute t value and on LLRT.
}
\description{
This function back-fits an inital LMER model on (absolute) t values and log-likelihood ratio test.

It works as follows:
  (1)  Set REML to FALSE;
  (2)  First consider only highest-order interaction model terms:
       (a) The model term with the lowest absolute t value is identified;
           (i) If this t value smaller than the threshold t value (set by 
               deault at 2), a new model without this model term is fitted;
          (ii) The more complex and simpler models are compared by way of a 
               log-likelihood ratio test (LLRT). If the results of the LLRT 
               determine that the term under consideration does not increase
               model fit, it is removed; otherwise it is kept.
         (iii) Move on to the next model term with the next smallest t value 
               below the threshold value and repeat steps (i)--(iii).
        (b) Once all highest-order interaction terms have been evaluated, 
            go down to the the second highest order of interaction: Repeat 
            steps (ai)--(aiii) with the following addition: If a term would be 
            removed from the model, but it is part of a high-order interaction,
            keep it. Once all terms of the interaction level have been 
            evaluated, move down to the next lower-order level until main 
            effects have been evaluated, after which the process stops.
  (3)  Set REML to TRUE (default) unless otherwise specified.
}
\usage{
bfFixefLMER_t.fnc(model = as.character(), data = as.character(), item = FALSE, alpha = 0.05, t.threshold = 2, set.REML.FALSE = TRUE, reset.REML.TRUE = TRUE, log.file = paste("./fixef_backfit_log_", gsub(" ", "_", date()), ".txt", sep = ""))
}
\arguments{
  \item{model}{
A ``mer'' object (fitted by function ``lmer'').
}
  \item{data}{
The data frame on which the ``mer'' object was fitted.
}
  \item{item}{
Whether or not to evaluate the addition of by-item random intercepts to the model, evaluated by way of log-likelihood ratio test. Either FALSE (the default, does not evaluate this addition) or the columnname (quoted) of the item identifier (e.g., "Item", "Word").
}
  \item{alpha}{
Level of significance (for log-likelihood ratio test). It is 0.05 by default.
}
  \item{t.threshold}{
Threshold for t value below which the log-likelihood ratio test is performed between a model with the term under consideration and a simpler model without it. Defaults to 2.
}
  \item{set.REML.FALSE}{
Logical. Whether or not to set REML to FALSE. Defaults to FALSE.
}
  \item{reset.REML.TRUE}{
Logical. Whether or not to re-set the back-fitted model to REML = TRUE.
}
  \item{log.file}{
Should the back-fitting log be saved? defaults to ``paste("./fixef_backfit_log_", gsub(" ", "_", date()), ".txt", sep = "")'', that is, the path and file name of the log, which can be changed to whatever the use wishes. Set to FALSE to disable.
}
}
\details{
Basically, a term remains in the model if (1) its t value is equal to or larger than the threshold t value or, in the case that it is not, if (2) it significantly increases model fit as determined by a log-likelihood ratio test (as in Pinheiro and Bates, 2000).

This backfitting method was used in Tremblay and Tucker (in preparation).
}
\value{
A ``mer'' model with back-fitted fixed effects (on t values) is returned and a log of the back-fitting process is printed on screen and (by default) in a log file.
}
\references{
Pinheiro, J.C. and Bates, D.M. (2000). Mixed Effects Models in S and S-Plus. New York: Springer.

Baayen, R.H., Davidson, D.J. and Bates, D.M. (2008). Mixed-effects modeling with crossed random effects for subjects and items. Journal of Memory and Language, 59, 390--412.

Tremblay, A. and Tucker B. V. (in preparation). What can the production of four-word sequences tell us about the mental lexicon? To be submitted to The Mental Lexicon.
}
\author{
Antoine Tremblay, IWK Health Center, <trea26@gmail.com>
}

\seealso{
\code{\link{bfFixefLMER_F.fnc}}
\code{\link{ffRanefLMER.fnc}}
\code{\link{fitLMER.fnc}}
\code{\link{pamer.fnc}}
\code{\link{mcp.fnc}}
\code{\link{romr.fnc}}
}
\examples{
# bfFixefLMER -- backfit on F
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct == "correct", ]
# fit initial model
m1 = lmer(RT ~ Frequency * FamilySize * Length * DerivEntropy + (1|Subject), data = lexdec3)
# look at model criticism plots
mcp.fnc(m1, lexdec3)
# trim data on initial model residuals
dat.trimmed = romr.fnc(m1, lexdec3)
lexdec3 = dat.trimmed$data
# update initial model on trimmed data
mB = update(m1)
# backfit fixed effects on F
m2 = bfFixefLMER_t.fnc(mB, lexdec3)
pamer.fnc(m2)

# An example with rcs
library(Design)
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct=="correct", ]
# base model
m0 = lmer(RT ~  rcs(Frequency, 4) * (Sex + NativeLanguage + Length) + (1|Subject), data = lexdec3)
# look at model criticism plots
mcp.fnc(m0, lexdec3)
# trim data on model residuals
lexdec3 = romr.fnc(m0, lexdec3)$data
# update model
m0 = update(m0)
# backfit fixed effects on F
m1 = bfFixefLMER_t.fnc(m0, lexdec3)
pamer.fnc(m1)

# An example with poly
library(Design)
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct=="correct", ]
# base model
m0 = lmer(RT ~  poly(Frequency, 4, raw = TRUE) * (Sex + NativeLanguage + FamilySize + Length) + (1|Subject), data = lexdec3)
# look at model criticism plots
mcp.fnc(m0, lexdec3)
# trim data on model residuals
lexdec3 = romr.fnc(m0, lexdec3)$data
# update model
m0 = update(m0)
# backfit fixed effects on F
m1 = bfFixefLMER_t.fnc(m0, lexdec3)
pamer.fnc(m1)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
