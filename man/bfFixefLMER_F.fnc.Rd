\name{bfFixefLMER_F.fnc}
\alias{bfFixefLMER_F.fnc}
\title{Back-fits an LMER model on p-values from ANOVA, llrt, AIC, BIC, relLik.AIC or relLik.BIC.}
\description{This function back-fits an initial LMER model either on upper- or lower-bound \emph{p}-values obtained from function \code{pamer.fnc}, log-likelihood ratio testing (LLRT), AIC, BIC, relLik.AIC, or relLik.BIC. Note that this function CANNOT be used with generalized linear mixed-effects models (\code{glmer}s).}
\usage{bfFixefLMER_F.fnc(model, item = FALSE, 
method = c("F", "llrt", "AIC", "BIC", "relLik.AIC", 
"relLik.BIC"), threshold = NULL, alpha = NULL, 
alphaitem = NULL, prune.ranefs = TRUE, 
p.value = "upper", set.REML.FALSE = TRUE, 
keep.single.factors=FALSE, reset.REML.TRUE = TRUE, 
log.file = NULL)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}). Note that
	  this function cannot be used with generalized linear mixed-effects models
	  (\code{glmer}s).}
  \item{item}{Whether or not to evaluate the addition of by-item random
	  intercepts to the model, evaluated by way of log-likelihood ratio test.
	  Either \code{FALSE} (the default) or the column name (quoted) of the item
	  identifier (e.g., \code{"Item"}, or \code{"Word"}).}
  \item{method}{Backfitting method. One of "F" (\emph{p}-value), "llrt", "AIC",
	  "BIC", "relLik.AIC", or "relLik.BIC" (relative likelihood, see function
	  \code{relLik}). Defaults to \code{F}. You can find information regarding
	  differences between AIC and BIC from
	  \code{http://methodology.psu.edu/eresources/ask/sp07}.}
  \item{threshold}{Method-specific threshold for parameter selection. It refers
	  to alpha in the case of "F" and "llrt", to the minimum reduction in
	  likelihood in the case of "AIC" and "BIC", or to the minimum difference
	  in probability in the case of "relLik.AIC" and "relLik.BIC". Defaults
	  \code{NULL}, which means \code{0.05} for "F" and "llrt", 5 for "AIC" and
	  "BIC", and 4 for "relLik.AIC" and "relLik.BIC".} \item{alpha}{If the
	  method is \code{F}, it is the \emph{p}-value (from \code{pamer.fnc})
	  above which a model term is dropped. In this case, it defaults to the
	  value passed to argument \code{threshold}, i.e., 0.05. Otherwise it is
	  the \emph{p}-value threshold above which a test (see \code{method}) is
	  performed between a model with the term under consideration and a simpler
	  model without it (in this case, defaults to \code{0}, i.e. all terms will
	  be tested).} \item{alphaitem}{Alpha value for the evaluation of by-item
	  random intercepts. Defaults to \code{0.05} or to the specified
	  threshold.} \item{prune.ranefs}{Logical. Whether to remove any random
	  effect for which its variable is not also present in the fixed effects
	  structure (with the exception of the grouping variables such as
	  \code{"Subjects"} and \code{"Items"}). Defaults to \code{TRUE}. For
	  example, if the random effects structure contains the terms
	  \code{Condition + ROI + Group}, and the random effects structure contains
	  the terms \code{(1 | Subject) + (0 + TrialNum | Subject)}, the random
	  effect \code{(0 + TrialNum | Subject)} will be pruned from the model
	  given that it is not in the model's fixed effects structure.}
  \item{p.value}{If \code{method = "F"}, whether to use upper-bound (``upper'';
	  the default) or lower-bound (``lower'') \emph{p}-values during
	  backfitting.} \item{set.REML.FALSE}{Logical. Whether or not to set
	  \code{REML} to \code{FALSE}. Defaults to \code{TRUE}.}
  \item{keep.single.factors}{Logical. Whether or not main effects are kept (not
	  subjected to testing and reduction). Defaults to \code{FALSE}.}
  \item{reset.REML.TRUE}{Logical. Whether or not to re-set the back-fitted
	  model to \code{REML = TRUE}.} \item{log.file}{Whether a back-fitting log
	  should be saved. Defaults to \code{NULL}, which means that a log is saved
	  in a temporary folder with the file name \code{file.path(tempdir(),
		  paste("bfFixefLMER_F_log_", gsub(":", "-", gsub(" ", "_", date())),
		  ".txt", sep = ""))}. The path and file name of the log can be changed
	  to whatever the use wishes. Set to \code{FALSE} to disable.}
}
\details{The back-fitting process works as follows: \enumerate{ \item If
		argument \code{method} is not set to \code{F}, \code{REML} is set to
		\code{FALSE}; \item First consider only highest-order interaction model
		terms: \enumerate{ \item If \code{method} is \code{F}, the model term
			with the highest ANOVA \emph{p}-value is identified. If this
			\emph{p}-value is higher than \code{alpha},the model term is
			removed and a new model is fitted. This is repeated for each model
			term that has a \emph{p}-value higher than the \code{alpha} value.
			The algorithm then moves on to step (b). If \code{method} is not
			\code{F}, the model term with the lowest \emph{p}-value is
			identified and the following is evaluated: \enumerate{ \item A new
				model without this model term is fitted; \item The more complex
				and simpler models are compared by way of a log-likelihood
				ratio test in case \code{method} is "llrt", by way of AIC or
				BIC values in case \code{method} is "AIC" or "BIC", or by
				calculating the \code{relLik} based on AIC or BIC in case
				\code{method} is "relLik.AIC" or "relLik.BIC". If the result
				determines that the term under consideration does not increase
				model fit, it is removed; otherwise it is kept.  \item Move on
				to the next model term with the smallest \emph{p}-value smaller
				than \code{alpha} and repeat steps (i)--(iii).  } \item Once
			all highest-order interaction terms have been evaluated, go down to
			the second highest order interactions: Repeat steps (ai)--(aiii)
			with the following addition: If a term would be removed from the
			model, but it is part of a high-order interaction, keep it. Once
			all terms of the interaction level have been evaluated, move down
			to the next lower-order level until main effects have been
			evaluated, after which the process stops. If \code{keep.single
				factors = TRUE}, the process stops after the evaluation of all
			interaction terms.  } \item If argument \code{method} is set to
		something else other than "F", set \code{reset.REML.TRUE} to
		\code{TRUE} (default) unless otherwise specified.  }

In brief, if \code{method} is set to "F", a term remains in the model if its
\emph{p}-value is equal to or greater than \code{alpha}; if \code{method} is
set to something else, a term remains in the model if \enumerate{ \item its
	\emph{p}-value from the ANOVA is equal to or smaller than \code{alpha};
	\item it significantly increases model fit as determined by the specified
	method; \item it is part of a significant higher-order interaction term.  }

This backfitting method was used in Newman, Tremblay, Nichols, Neville, and
Ullman (2012). If factorial terms are included in the initial model,
back-fitting on F is recommended.  } \note{If you get this error:
\preformatted{

Error in model.frame.default(data = ..2, formula = log_Segment_Duration ~  :
The ... list does not contain 2 elements }

It is probably because you updated the model using function \code{update} and
the data now appears as \code{data = ..2} or something similar to this. You can
check this by typing \code{model@call}. If this is the case, re-fit your model
as \code{lmer(DV ~ IV + IV + (RANEF), data = dat)}.} \value{A \code{mer} model
with back-fitted fixed effects is returned and a log of the back-fitting
process is printed on screen and (by default) in a log file in a temporary
file.} \references{ Newman, A.J., Tremblay, A., Nichols, E.S., Neville, H.J.,
and Ullman, M.T. (2012). The Influence of Language Proficiency on Lexical
Semantic Processing in Native and Late Learners of English. \emph{Journal of
	Cognitive Neuroscience}, \emph{25}, 1205--1223.  } \author{Antoine
Tremblay, Statistics Canada, \email{trea26@gmail.com} and Johannes Ransijn
\email{johannesransijn@gmail.com}.}

\section{Warnings}{Upper-bound \emph{p}-values can be anti-conservative, while
	lower-bound \emph{p}-values can be conservative. See function
	\code{pamer.fnc}.}

\seealso{
\code{\link{bfFixefLMER_t.fnc}; }
\code{\link{ffRanefLMER.fnc}; }
\code{\link{fitLMER.fnc}; }
\code{\link{mcposthoc.fnc}; }
\code{\link{pamer.fnc}; }
\code{\link{mcp.fnc}; }
\code{\link{relLik}; }
\code{\link{romr.fnc}; }
\code{\link{perSubjectTrim.fnc}.}
}
\examples{
# see example in LMERConvenienceFunctions help page.
}
\keyword{models & regression}
