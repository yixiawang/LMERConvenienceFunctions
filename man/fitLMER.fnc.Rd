\name{fitLMER.fnc}
\alias{fitLMER.fnc}
\title{Back-fit fixed effects and forward-fit random effects of an LMER model.}
\description{The function follows these steps: (1) If \code{llrt} is set to \code{TRUE}, set \code{REML} to \code{FALSE} (unless specified otherwise); (2) back-fit initial model either on \emph{F}- (by default) or on \emph{t}-values; (3) forward-fit random effects; (4) re-back-fit fixed effects; (5) if \code{llrt} is set to \code{TRUE}, set \code{REML} to \code{TRUE} (unless specified otherwise).}
\usage{
fitLMER.fnc(model = as.character(), data = as.character(), 
backfit.on = "F", item = FALSE, ran.effects = list(ran.intercepts= 
as.character(), slopes = as.character(), by.vars = as.character()), 
alpha = 0.05, llrt = FALSE, p.value = "Upper", t.threshold = 2, 
set.REML.FALSE = TRUE, reset.REML.TRUE = TRUE, file.name = 
paste("full_backfitLMER_log_", gsub(":", "_", gsub(" ", "_", date())), 
".txt", sep=""))
}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}).}
  \item{data}{The data frame on which the \code{mer} object was fitted.}
  \item{backfit.on}{Either back-fit on \code{F}- (default) or on \code{t}-values.}
  \item{item}{Whether or not to evaluate the addition of by-item random intercepts to the model, evaluated by way of log-likelihood ratio test. Either \code{FALSE} (the default, does not evaluate this addition) or the columnname (quoted) of the item identifier (e.g., \code{"Item"}, \code{"Word"}).}
  \item{ran.effects}{Can be either a vector or a list. In the former case, the random effects to be evaluated are provided. For example \code{c("(1 + Frequency | Subject)", "(0 + Length | Subject)", "(1 + NSynSet | Subject)")}. In the latter case, the list can be composed of (i) a vector of random intercepts to be evaluated, (ii) a vector of random slopes to be evaluated, and (iii) a vector of ``by'' variables for the random slopes. For example \code{list(ran.intercepts = "WordType", slopes = c("Frequency", "Length", "NSynSet"), by.vars = "Subject")}. In this case, correlations between random intercepts and slopes are automatically included.}
  \item{alpha}{Level of significance (for ANOVA and log-likelihood ratio test). It is 0.05 by default.}
  \item{llrt}{Logical. Whether log-likelihood ratio testing between the complex and simple models should be performed. Defaults to \code{FALSE}.}
  \item{p.value}{Whether to use upper-bound (``Upper''; the default) or lower-bound (``Lower'') \emph{p}-values during backfitting.}
  \item{t.threshold}{Threshold for \emph{t}-value below which the log-likelihood ratio test is performed between a model with the term under consideration and a simpler model without it. Defaults to 2.}
  \item{set.REML.FALSE}{Logical. Whether or not to set \code{REML} to \code{FALSE}. Defaults to \code{FALSE}.}
  \item{reset.REML.TRUE}{Logical. Whether or not to re-set the back-fitted model to \code{REML = TRUE}.}
  \item{file.name}{Should the back-fitting log be saved? By default, it is saved in a temporary folder (platform dependent) as \code{paste("full_backfitLMER_log_", gsub(":", "_", gsub(" ", "_", date())), ".txt", sep = "")}. The path and file name of the log can be changed to whatever the use wishes. Set to \code{FALSE} to disable.}
}
\details{The process has three stages. In the first stage, either \code{bfFixefLMER_F.fnc} or \code{bfFixefLMER_t.fnc} is called (depending on the user's choice) and the fixed effects are back-fitted accordingly. In the second stage, \code{ffRanefLMER.fnc} is called and random effects are forward-fitted. In the third stage, the fixed effects are back-fitted again. This is done because the inclusion of certain random effects sometimes renders certain fixed effects non-significant. Note that if this results in a random effect being present in the model with no corresponding fixed effect (an undesirable situation), the user will have to manually remove this random effect. This process was used in Tremblay and Tucker (submitted) and in Newman, Tremblay, Nichols, Neville, and Ullman (submitted).}
\value{A \code{mer} object with fixed effects (on \emph{F}- or \code{t}-values) and back-fitted forward-fitted random effect structure as well as a log of the process is printed on screen and, optionally, printed in a log file saved in a temporary file (platform specific).}
\references{
Baayen, R.H., Davidson, D.J. and Bates, D.M. (2008). Mixed-effects modeling with crossed random effects for subjects and items. Journal of Memory and Language, 59, 390--412.

Newman, A.J., Tremblay, A., Nichols, E.S., Neville, H.J., and Ullman, M.T. (submitted). The Influence of Language Proficiency on Lexical-Semantic Processing in Native and Late Learners of English: ERP evidence. Submitted to the \emph{Journal of Cognitive Neuroscience}.

Pinheiro, J.C. and Bates, D.M. (2000). \emph{Mixed Effects Models in S and S-Plus}. New York: Springer.

Tremblay, A. and Tucker B. V. (submitted). What can the production of four-word sequences tell us about the mental lexicon? Submitted to \emph{The Mental Lexicon}.
}
\author{Antoine Tremblay, IWK Health Center, \email{trea26@gmail.com}}

\section{Warnings}{Upper-bound \emph{p}-values can be anti-conservative, while lower-bound \emph{p}-values can be conservative. See \url{http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests} and function \code{pamer.fnc}.}

\seealso{
\code{\link{bfFixefLMER_F.fnc}}
\code{\link{bfFixefLMER_t.fnc}}
\code{\link{ffRanefLMER.fnc}}
\code{\link{posthoc.fnc}}
\code{\link{pamer.fnc}}
\code{\link{mcp.fnc}}
\code{\link{romr.fnc}}
\code{\link{perSubjectTrim.fnc}}
}
\examples{
# Example 1
library(lme4)
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct == "correct", ]

# fit initial model
m0 = lmer(RT ~ Sex + NativeLanguage + PrevType + Frequency +
    FamilySize + Length + DerivEntropy + (1|Subject), 
    data = lexdec3)

# look at model criticism plots
mcp.fnc(m0, trim = 2.5)

# trim data on initial model residuals
dat.trimmed = romr.fnc(m0, lexdec3, trim = 2.5)
lexdec3 = dat.trimmed$data

# update initial model on trimmed data
m0 = update(m0)

# fitLMER #1 -- backfit on t
mA = fitLMER.fnc(m0, lexdec3,item = "Word", ran.effects = 
    list(ran.intercepts = "Class", slopes = c("Frequency", 
    "Length", "FamilySize"), by.vars = "Subject"), 
    file.name = FALSE, backfit.on = "t")
pamer.fnc(mA)

# fitLMER #2 -- backfit on F
mB = fitLMER.fnc(m0, lexdec3, item = "Word", ran.effects = 
    c("(1+Length|Subject)", "(1+Frequency|Subject)", 
    "(1+FamilySize|Subject)"), file.name = FALSE)
pamer.fnc(mB)

# Example 2
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct=="correct", ]

# fit initial model
m0 = lmer(RT ~ Frequency * FamilySize * Length * 
    DerivEntropy + (1|Subject), data = lexdec3)

# look at criticism plots
mcp.fnc(m0, trim = 2.5)

# trim data on initial model residuals
dat.trimmed = romr.fnc(m0, lexdec3, trim = 2.5)
lexdec3 = dat.trimmed$data

# update initial model on trimmed data
m0 = update(m0)

# fitLMER #3 -- backfit on F
mA = fitLMER.fnc(m0, lexdec3, item = "Word", ran.effects = 
    list(slopes = c("Frequency", "FamilySize", "Length", 
    "DerivEntropy"), by.vars = "Subject"), file.name = FALSE)

# look at model criticism plots
mcp.fnc(mA, trim = 2.5)
pamer.fnc(mA)
}
\keyword{models & regression}
