\name{fitLMER.fnc}
\alias{fitLMER.fnc}
\title{Back-fit fixed effects and forward-fit random effects of an LMER model.}
\description{The function follows these steps: (1) If \code{llrt} is set to \code{TRUE}, set \code{REML} to \code{FALSE} (unless specified otherwise); (2) back-fit initial model either on \emph{F}- (by default) or on \emph{t}-values; (3) forward-fit random effects; (4) re-back-fit fixed effects; (5) if \code{llrt} is set to \code{TRUE}, set \code{REML} to \code{TRUE} (unless specified otherwise). Note that, at the moment, this function cannot be used with generalized linear mixed-effects models (\code{glmer}s).}
\usage{fitLMER.fnc(model, data, backfit.on = "F", item = FALSE, 
ran.effects = list(ran.intercepts= as.character(), 
slopes = as.character(), corr = as.character(), 
by.vars = as.character()), alpha = 0.05, if.warn.not.add = TRUE, 
llrt = FALSE, prune.ranefs = TRUE, p.value = "upper", 
t.threshold = 2, set.REML.FALSE = TRUE, reset.REML.TRUE = TRUE, 
log.file.name = file.path(tempdir(), paste("fitLMER_log_", 
gsub(":", "-", gsub(" ", "_", date())), ".txt", sep="")))}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}). Note that, at the moment, this function cannot be used with generalized linear mixed-effects models (\code{glmer}s).}
  \item{data}{The data frame on which the \code{mer} object was fitted.}
  \item{backfit.on}{Either back-fit on \code{F}- (default) or on \code{t}-values.}
  \item{item}{Whether or not to evaluate the addition of by-item random intercepts to the model, evaluated by way of log-likelihood ratio test. Either \code{FALSE} (the default, does not evaluate this addition) or the columnname (quoted) of the item identifier (e.g., \code{"Item"}, \code{"Word"}).}
  \item{ran.effects}{Can be either a vector or a list. In the former case, the random effects to be evaluated are provided. For example \code{c("(1 + Frequency | Subject)", "(0 + Length | Subject)", "(1 + NSynSet | Subject)")}. In the latter case, the list can be composed of (i) a vector of random intercepts to be evaluated (\code{ran.intercepts}), (ii) a vector of random slopes to be evaluated (\code{slopes}), (iii) a vector specifying, for each element of \code{slopes}, whether the correlation between the slope and by-variables specified in \code{by.vars} should be added (\code{corr}), and (iv) a vector of ``by'' variables for the random slopes (\code{by.vars}). Values that can be supplied to the \code{corr} argument are \code{1} (add correlation), \code{0} (do not add correlation), and \code{NA} (for when the "slope" is a factor variable). Note that if a term in \code{slopes} is a factor variable, the \code{corr} value tied to it will be automatically set to \code{NA}. Also note that if no values are supplied to \code{corr}, a vector of 0 as long as the \code{slopes} vector will be automatically supplied. For example \code{list(ran.intercepts = "Word", slopes = c("Frequency", "Length", "NSynSet","Class"), corr = c(0, 0, 1, NA), by.vars = "Subject")}. Another example is \code{list(slopes = c("Trial", "Class"), by.vars = "Subject")}, where the \code{corr} argument will be equal to \code{c(0, NA)}.}
  \item{alpha}{Level of significance (for ANOVA and log-likelihood ratio test). It is 0.05 by default.}
  \item{if.warn.not.add}{Logical. If a warning is issued after fitting a model with a new random effect (e.g., \code{false convergence} or the like), should the random effect nevertheless be evaluated? Defaults to \code{TRUE}, meaning that if such a warning is issued, the random effect will not be added to the random effects structure of the model. If set to \code{FALSE}, the random effect will be evaluated for inclusion as any other random effects would be via log likelihood ratio testing even if a warning is issued.}
  \item{llrt}{Logical. Whether log-likelihood ratio testing between the complex and simple models should be performed. Defaults to \code{FALSE}.}
  \item{prune.ranefs}{Logical. Whether to remove any random effect for which its variable is not also present in the fixed effects structure (with the exception of the grouping variables such as \code{"Subjects"} and \code{"Items"}). Defaults to \code{TRUE}. For example, if the random effects structure contains the terms \code{Condition + ROI + Group}, and the random effects structure contains the terms \code{(1 | Subject) + (0 + TrialNum | Subject)}, the ranedom effect \code{(0 + TrialNum | Subject)} will be pruned from the model given that it is not in the model's fixed effects structure.}
  \item{p.value}{Whether to use upper-bound (``upper''; the default) or lower-bound (``lower'') \emph{p}-values during back-fitting.}
  \item{t.threshold}{Threshold for \emph{t}-value below which either a model terms is dropped or a log-likelihood ratio test is performed between a model with the term under consideration and a simpler model without it. Defaults to 2. If data frame has fewer than ~ 100 observations, this values should be adjusted.}
  \item{set.REML.FALSE}{Logical. Whether or not to set \code{REML} to \code{FALSE}. Defaults to \code{FALSE}.}
  \item{reset.REML.TRUE}{Logical. Whether or not to re-set the back-fitted model to \code{REML = TRUE}.}
  \item{log.file.name}{Should the back-fitting log be saved? By default, it is saved in a temporary folder (platform dependent) as \code{file.path(tempdir(), paste("fitLMER_log_", gsub(":", "-", gsub(" ", "_", date())), ".txt", sep = ""))}. The path and file name of the log can be changed to whatever the use wishes. Set to \code{FALSE} to disable.}
}
\details{The process has three stages. In the first stage, either \code{bfFixefLMER_F.fnc} or \code{bfFixefLMER_t.fnc} is called (depending on the user's choice) and the fixed effects are back-fitted accordingly. In the second stage, \code{ffRanefLMER.fnc} is called and random effects are forward-fitted. In the third stage, the fixed effects are back-fitted again. This is done because the inclusion of certain random effects sometimes renders certain fixed effects non-significant. This process was used in Tremblay and Tucker (Accepted) and in Newman, Tremblay, Nichols, Neville, and Ullman (Resubmitted).}
\value{A \code{mer} object with fixed effects (on \emph{F}- or \code{t}-values) and back-fitted forward-fitted random effect structure as well as a log of the process is printed on screen and, optionally, printed in a log file.}
\references{
Baayen, R.H., Davidson, D.J. and Bates, D.M. (2008). Mixed-effects modeling with crossed random effects for subjects and items. Journal of Memory and Language, 59, 390--412.

Newman, A.J., Tremblay, A., Nichols, E.S., Neville, H.J., and Ullman, M.T. (Resubmitted). The Influence of Language Proficiency on Lexical-Semantic Processing in Native and Late Learners of English: ERP evidence. Submitted to the \emph{Journal of Cognitive Neuroscience}.

Pinheiro, J.C. and Bates, D.M. (2000). \emph{Mixed Effects Models in S and S-Plus}. New York: Springer.

Tremblay, A. and Tucker B. V. (Accepted). The Effects of N-gram Probabilistic Measures on the Processing and Production of Four-word Sequences. \emph{The Mental Lexicon}.
}
\author{Antoine Tremblay, IWK Health Center, \email{trea26@gmail.com}}

\section{Warnings}{Upper-bound \emph{p}-values can be anti-conservative, while lower-bound \emph{p}-values can be conservative. See \url{http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests} and function \code{pamer.fnc}.}

\seealso{
\code{\link{bfFixefLMER_F.fnc}; }
\code{\link{bfFixefLMER_t.fnc}; }
\code{\link{ffRanefLMER.fnc}; }
\code{\link{posthoc.fnc}; }
\code{\link{pamer.fnc}; }
\code{\link{mcp.fnc}; }
\code{\link{romr.fnc}; }
\code{\link{perSubjectTrim.fnc}; }
}
\examples{
if("LCFdata"\%in\%.packages(all.available=TRUE)){
  library(LCFdata)
	library(lme4)
	# bfFixefLMER -- backfit on t
	data(eeg)
	eeg <- eeg[eeg$Time >= 100 & eeg$Time <= 175, , drop = TRUE]
	
	# fit initial model
	m1 = lmer(Fz ~ FreqB * LengthB * WMC + (1 | Subject) + 
	     (1 | Item), data = eeg)
	
	# look at model criticism plots
	mcp.fnc(m1, trim = 2.5)
	
	# trim data on initial model residuals
	eeg.trimmed = romr.fnc(m1, eeg, trim = 2.5)
	eeg = eeg.trimmed$data
	
	# update initial model on trimmed data
	mB = update(m1)
	dev.new()
	mcp.fnc(mB, trim = 2.5)
	
	# forward and backward fit model on t
	m2<-fitLMER.fnc(mB, eeg, backfit.on = "t", item = FALSE, 
	ran.effects = c("(FreqB | Subject)",
	"(LengthB | Subject)", "(WMC | Item)"))
	
	devAskNewPage(ask=FALSE)
}
}
\keyword{models & regression}
