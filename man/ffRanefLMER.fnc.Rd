\name{ffRanefLMER.fnc}
\alias{ffRanefLMER.fnc}
\title{Forward-fit the random effect structure of an LMER model.}
\description{Forward-fit an LMER model's random effect structure by comparing a model without one of the specified random effects and a model with it by way of log-likelihood ratio testing. If the more complex model is a significantly better fit, the random effect is kept, otherwise it is dropped. Note that, at the moment, this function cannot be used with generalized linear mixed-effects models (\code{glmer}s).}
\usage{ffRanefLMER.fnc(model, data, ran.effects = 
list(ran.intercepts = as.character(), slopes = 
as.character(), by.vars = as.character()), alpha = 0.05, 
log.file = file.path(tempdir(), paste("ffRanefLMER_log_", 
gsub(":", "-", gsub(" ", "_", date())), ".txt", sep="")))}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}). Note that, at the moment, this function cannot be used with generalized linear mixed-effects models (\code{glmer}s).}
  \item{data}{The data frame on which the \code{mer} object was fitted.}
  \item{ran.effects}{Can be either a vector or a list. In the former case, the random effects to be evaluated are provided. For example \code{c("(1 + Frequency | Subject)", "(0 + Length | Subject)", "(1 + NSynSet | Subject)")}. In the latter case, the list can be composed of (i) a vector of random intercepts to be evaluated, (ii) a vector of random slopes to be evaluated, and (iii) a vector of ``by'' variables for the random slopes. For example \code{list(ran.intercepts = "WordType", slopes = c("Frequency", "Length", "NSynSet"), by.vars = "Subject")}. In this case, correlations between random intercepts and slopes are automatically included.
}
  \item{alpha}{Level of significance for log-likelihood ratio test. Defaults to 0.05.}
  \item{log.file}{Should the back-fitting log be saved? By default, it is saved in a temporary folder as ``paste("ffRanefLMER_log_", gsub(":", "-", gsub(" ", "_", date())), ".txt", sep = "")''. The path and file name of the log can be changed to whatever the use wishes. Set to \code{FALSE} to disable.}
}
\value{A \code{mer} object with forward-fitted random effect structure as well as a log of the process is printed on screen and, optionally, printed in a log file.}
\references{Pinheiro, J.C. and Bates, D.M. (2000). \emph{Mixed Effects Models in S and S-Plus}. New York: Springer.}
\author{Antoine Tremblay, IWK Health Center, \code{trea26@gmail.com}.}

\seealso{
\code{\link{bfFixefLMER_F.fnc}; }
\code{\link{bfFixefLMER_t.fnc}; }
\code{\link{fitLMER.fnc}; }
\code{\link{posthoc.fnc}; }
\code{\link{pamer.fnc}; }
\code{\link{mcp.fnc}; }
\code{\link{romr.fnc}; }
\code{\link{perSubjectTrim.fnc}; }
}
\examples{
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
mA = ffRanefLMER.fnc(m0, lexdec3, ran.effects = list(ran.intercepts = 
    "Class", slopes = c("Frequency", "Length", "FamilySize"), 
    by.vars = "Subject"), log.file = FALSE)
print(mA, corr = FALSE)
mB = ffRanefLMER.fnc(m0, lexdec3, ran.effects = c("(1+Length|Subject)", 
    "(1+Frequency|Subject)", "(1+FamilySize|Subject)"), log.file = FALSE)
print(mB, corr = FALSE)
}
\keyword{models & regression}
