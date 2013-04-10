\name{bfFixefLMER_t.fnc}
\alias{bfFixefLMER_t.fnc}
\title{Back-fits an LMER model on absolute t-value and, optionally, on LLRT.}
\description{This function back-fits an initial LMER model on \emph{t}-values, and, if enabled, log-likelihood ratio testing. Note that, this function CAN be used with generalized linear mixed-effects models (\code{glmer}s).}
\usage{bfFixefLMER_t.fnc(model, item = FALSE,
method = c("t", "llrt", "AIC", "BIC", "relLik.AIC", 
"relLik.BIC"), threshold = NULL, t.threshold = NULL, 
alphaitem = NULL, prune.ranefs = TRUE, set.REML.FALSE = TRUE, 
keep.single.factors=FALSE, reset.REML.TRUE = TRUE, 
log.file = NULL)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}). Note that this function can be used with generalized linear mixed-effects models (\code{glmer}s).}
  \item{item}{Whether or not to evaluate the addition of by-item random intercepts to the model, evaluated by way of log-likelihood ratio test. Either \code{FALSE} (the default) or the column name (quoted) of the item identifier (e.g., \code{"Item"}, or \code{"Word"}).}
  \item{method}{Backfitting method. One of "t" (\emph{p}-value), "llrt", "AIC", "BIC", "relLik.AIC", or "relLik.BIC" (relative likelihood, see function \code{relLik}). Defaults to "t".}  
  \item{threshold}{Method-specific threshold for parameter selection. It refers to the minimum \emph{t}-value in the case of "t", to the alpha value in the case of "llrt", to the minimum reduction in likelihood in the case of "AIC" and "BIC", or to the minimum difference in probability in the case of "relLik.AIC" and "relLik.BIC". Defaults \code{NULL}, which means \code{2} for "t", \code{0.05} for "llrt", 5 for "AIC" and "BIC", and 4 for "relLik.AIC" and "relLik.BIC".}
  \item{t.threshold}{Defaults to \code{NULL}. If the \code{method = "t"}, it is the \emph{t}-value below which a model term is dropped (if \code{t.threshold = NULL}, it will be set to 2). Otherwise it is the threshold for \emph{t}-value below which a test (see \code{method}) is performed between a model with the term under consideration and a simpler model without it (if \code{t.threshold = NULL}, it is set to \code{Inf}, which means that all terms are tested.} 
  \item{alphaitem}{Alpha value for the evaluation of by-item random intercepts. Defaults to \code{0.05} or to the specified threshold in case \code{method} is \code{llrt}.}
  \item{prune.ranefs}{Logical. Whether to remove any random effect for which its variable is not also present in the fixed effects structure (with the exception of the grouping variables such as \code{"Subjects"} and \code{"Items"}). Defaults to \code{TRUE}. For example, if the random effects structure contains the terms \code{Condition + ROI + Group}, and the random effects structure contains the terms \code{(1 | Subject) + (0 + TrialNum | Subject)}, the random effect \code{(0 + TrialNum | Subject)} will be pruned from the model given that it is not in the model's fixed effects structure.}
  \item{set.REML.FALSE}{Logical. Whether or not to set REML to \code{FALSE}. Defaults to \code{TRUE}.}
  \item{reset.REML.TRUE}{Logical. Whether or not to re-set the back-fitted model to \code{REML = TRUE}.}
  \item{keep.single.factors}{Logical. Whether or not main effects are kept (not subjected to testing and reduction). Defaults to \code{FALSE}.} 
  \item{log.file}{Whether a back-fitting log should be saved. Defaults to \code{NULL}, which means that a log is saved in a temporary folder with the file name \code{file.path(tempdir(), paste("bfFixefLMER_F_log_", gsub(":", "-", gsub(" ", "_", date())), ".txt", sep = ""))}. The path and file name of the log can be changed to whatever the use wishes. Set to \code{FALSE} to disable.}
}
\details{The back-fitting process works as follows:
\enumerate{
\item If argument \code{method} is not set to "t", \code{REML} is set to \code{FALSE};
\item First consider only highest-order interaction model terms:
	\enumerate{
		\item If \code{method} is "t", the model term with the lowest \emph{t}-value is identified. If this \emph{t}-value is smaller than \code{threshold}, the model term is removed and a new model is fitted. This is repeated for each model term for term that has a \emph{t}-value smaller than the threshold value. The algorithm then moves on to step (b). If \code{method} is not "t", the model term with the lowest \emph{t}-value is identified and the following is evaluated:
	\enumerate{
    	\item A new model without this model term is fitted; 
		\item The more complex and simpler models are compared by way of a log-likelihood ratio test in case \code{method} is "llrt", by way of AIC or BIC comparison if \code{method} is "AIC" "BIC", or by calculating the \code{relLik} based on AIC or BIC in case \code{method} is "relLik.AIC" or "relLik.BIC". If the result determines that the term under consideration does not increase model fit, it is removed; otherwise it is kept.
		\item Move on to the next model term with the smallest \emph{t}-value smaller than \code{threshold} and repeat steps (i)--(iii).
	}
    \item Once all highest-order interaction terms have been evaluated, go down to the second highest order interactions: Repeat steps (ai)--(aiii) with the following addition: If a term would be removed from the model, but it is part of a high-order interaction, keep it. Once all terms of the interaction level have been evaluated, move down to the next lower-order level until main effects have been evaluated, after which the process stops. If \code{keep.single factors = TRUE}, the process stops after the evaluation of all interaction terms.
  }
\item If argument \code{method} is set to something other than \code{t}, set \code{reset.REML.TRUE} to \code{TRUE} (default) unless otherwise specified.
}

In brief, if \code{method} is set to "t", a term remains in the model if its \emph{t}-value is equal to or greater than \code{threshold}; if \code{method} is set to something else, a term remains in the model if
\enumerate{
	\item its \emph{t}-value is equal to or greater than \code{threshold};
	\item it significantly increases model fit as determined by the specified method;
	\item it is part of a significant interaction term.
}

This backfitting method was used in Tremblay & Tucker (2011). If factorial terms with more than two levels are included in the initial model, back-fitting on F is recommended. 
}
\note{If you get this error:
\preformatted{

Error in model.frame.default(data = ..2, formula = log_Segment_Duration ~  : 
  The ... list does not contain 2 elements
}

It is probably because you updated the model using function \code{update} and the data now appears as \code{data = ..2} or something similar to this. You can check this by typing \code{model@call}. If this is the case, re-fit your model as \code{lmer(DV ~ IV + IV + (RANEF), data = dat)}.}
\value{A \code{mer} model with back-fitted fixed effects (on \code{t}-values) is returned and a log of the back-fitting process is printed on screen and (by default) in a log file.}
\references{
Tremblay, A. and Tucker B. V. (2011). The Effects of N-gram Probabilistic Measures on the Processing and Production of Four-word Sequences. \emph{The Mental Lexicon}, \emph{6(2)}, 302--324.
}
\author{Antoine Tremblay, Dalhousie University, \email{trea26@gmail.com} and Johannes Ransijn \email{johannesransijn@gmail.com}.}
\seealso{
\code{\link{bfFixefLMER_F.fnc}; }
\code{\link{ffRanefLMER.fnc}; }
\code{\link{fitLMER.fnc}; }
\code{\link{mcposthoc.fnc}; }
\code{\link{pamer.fnc}; }
\code{\link{mcp.fnc}; }
\code{\link{relLik}; }
\code{\link{romr.fnc}; }
\code{\link{perSubjectTrim.fnc}. }
}

\examples{
# see example in LMERConvenienceFunctions help page.
}
\keyword{models & regression}
