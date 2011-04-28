\name{LMERConvenienceFunctions-package}
\alias{LMERConvenienceFunctions-package}
\alias{LMERConvenienceFunctions}
\docType{package}
\title{An suite of functions to facilitate modeling with LMER.}
\description{The main functions are fixed effect back-fitting functions (\code{bfFixefLMER_F.fnc} to back-fit on \emph{p}-values from ANOVA or \code{bfFixefLMER.fnc_t.fnc} \emph{t}-values from \code{summary(model)} and, optionally, log-likelihood ratio testing) and a random effect forward-fitting function (\code{ffRanefLMER.fnc}; log-likelihood ratio testing). There is also a function to first back-fit fixed effects from an initial model, then forward-fit random effects, and finally re-backfit fixed effects (\code{fitLMER.fnc}). Other functions include a function to compute ANOVAs with upper- or lower-bound \emph{p}-values; \code{pamer.fnc}), a function to graph model criticism plots (\code{mcp.fnc}), a function to trim data on model residuals (\code{romr.fnc}), one to perform per-subject trimming on the response variable (\code{perSubjectTrim.fnc}, and a function to perform posthoc analyses (\code{posthoc.fnc}). Additionally, there are functions to list files in current directory in matrix format (\code{f}, with easily readable numbers for each file/directory), change directory using the \code{f} function (\code{cd}), and to automatically list files in new directory (\code{cdf}), to go up one directory and automatically list files (\code{cdup}), and finally a function to list in matrix format (\code{cn}, with easily readable numbers for each column name).
}
\details{
\tabular{ll}{
Package: \tab LMERConvenienceFunctions\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 28-04-2011\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
bfFixefLMER_F.fnc
bfFixefLMER_t.fnc
ffRanefLMER.fnc
fitLMER.fnc 
posthoc.fnc
pamer.fnc
mcp.fnc
romr.fnc
perSubjecttrim.fnc
cn
f
cd
cdf
cdup
}
\author{Antoine Tremblay, IWK Health Center 

Maintainer: Antoine Tremblay <trea26@gmail.com>
}

\references{
Baayen, R.H. (2008). \emph{Analyzing Linguistic Data. A Practical Introduction to Statistics Using R}. Cambridge, UK: Cambridge University Press. 

Baayen, R.H., Davidson, D.J. and Bates, D.M. (2008). Mixed-effects modeling with crossed random effects for subjects and items. \emph{Journal of Memory and Language}, \emph{59}, 390--412.

Newman, A.J., Tremblay, A., Nichols, E.S., Neville, H.J., and Ullman, M.T. (submitted). The Influence of Language Proficiency on Lexical-Semantic Processing in Native and Late Learners of English: ERP evidence. Submitted to the \emph{Journal of Cognitive Neuroscience}.

Pinheiro, J.C. and Bates, D.M. (2000). \emph{Mixed Effects Models in S and S-Plus}. New York: Springer.

Tremblay, A. and Tucker B. V. (submitted). What can the production of four-word sequences tell us about the mental lexicon? Submitted to \emph{The Mental Lexicon}.

\url{http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests}
}
\keyword{ package }
\seealso{
\code{\link{bfFixefLMER_F.fnc}}
\code{\link{bfFixefLMER_t.fnc}}
\code{\link{ffRanefLMER.fnc}}
\code{\link{fitLMER.fnc}}
\code{\link{posthoc.fnc}}
\code{\link{pamer.fnc}}
\code{\link{mcp.fnc}}
\code{\link{romr.fnc}}
\code{\link{perSubjectTrim.fnc}}
\code{\link{cn}}
\code{\link{f}}
\code{\link{cd}}
\code{\link{cdf}}
\code{\link{cdup}}
}
\examples{
# Example 1
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

# fitLMER #1 -- backfit on F
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
m0 = lmer(RT ~ Frequency * FamilySize * Length * DerivEntropy + 
    (1|Subject), data = lexdec3)

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

# bfFixefLMER -- backfit on F
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct == "correct", ]

# fit initial model
m0 = lmer(RT ~ Frequency * FamilySize * DerivEntropy +
    (1|Subject), data = lexdec3)

# look at model criticism plots
mcp.fnc(m0, trim = 2.5)

# trim data on initial model residuals
dat.trimmed = romr.fnc(m0, lexdec3, trim = 2.5)
lexdec3 = dat.trimmed$data

# update initial model on trimmed data
m0 = update(m0)

# backfit fixed effects on F
mB = bfFixefLMER_F.fnc(m0, lexdec3, log.file = FALSE)
pamer.fnc(mB)

# bfFixefLMER -- backfit on t
m1 = bfFixefLMER_F.fnc(m0, lexdec3, log.file = FALSE)
pamer.fnc(m1)

# Example 3
# An example with rcs
library(Design)
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
m1 = bfFixefLMER_F.fnc(m0, lexdec3, log.file = FALSE)
pamer.fnc(m1)

# Example 4
# An example with poly
library(Design)
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
m1 = bfFixefLMER_F.fnc(m0, lexdec3, log.file = FALSE)
pamer.fnc(m1)

data(lexdec)
cn(lexdec)

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
mA = ffRanefLMER.fnc(m0, lexdec3, ran.effects = list(ran.intercepts = 
    "Class", slopes = c("Frequency", "Length", "FamilySize"), 
    by.vars = "Subject"), log.file = FALSE)
print(mA, corr = FALSE)
mB = ffRanefLMER.fnc(m0, lexdec3, ran.effects = c("(1+Length|Subject)", 
    "(1+Frequency|Subject)", "(1+FamilySize|Subject)"), log.file = FALSE)
print(mB, corr = FALSE)

library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct == "correct", ]
m0 = lmer(RT ~ Sex + NativeLanguage + PrevType + Frequency + FamilySize + 
    Length + DerivEntropy + (1|Subject), data = lexdec3)
pamer.fnc(m0)

library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct == "correct", ]
m0 = lmer(RT ~ Sex + NativeLanguage + PrevType + Frequency + FamilySize + 
    Length + DerivEntropy + (1|Subject), data = lexdec3)
out=romr.fnc(m0,lexdec3, trim = 2.5)
names(out)

library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct == "correct", ]
m0 = lmer(RT ~ Sex + NativeLanguage + PrevType + Frequency + FamilySize + 
    Length + DerivEntropy + (1|Subject), data = lexdec3)
mcp.fnc(m0, trim = 2.5)
dat.trimmed = romr.fnc(m0, lexdec3, trim = 2.5)
names(dat.trimmed)
lexdec3 = dat.trimmed$data
m0 = update(m0)
mcp.fnc(m0, trim = 2.5)
}
