\name{bfFixefLMER_F.fnc}
\alias{bfFixefLMER_F.fnc}
\title{Back-fits an LMER model on p-values from ANOVA and, optionally, on LLRT.}
\description{This function back-fits an initial LMER model on upper- or lower-bound \emph{p}-values obtained from function \code{pamer.fnc} and, if enabled, log-likelihood ratio testing. Note that, at the moment, this function cannot be used with generalized linear mixed-effects models (\code{glmer}s).}
\usage{bfFixefLMER_F.fnc(model, item = FALSE, 
alpha = 0.05, llrt = FALSE, prune.ranefs = TRUE, 
p.value = "upper", set.REML.FALSE = TRUE, 
reset.REML.TRUE = TRUE, log.file = file.path(tempdir(), 
paste("bfFixefLMER_F_log_", gsub(":", "-", 
gsub(" ", "_", date())), ".txt", sep="")))}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}). Note that, at the moment, this function cannot be used with generalized linear mixed-effects models (\code{glmer}s).}
  \item{item}{Whether or not to evaluate the addition of by-item random intercepts to the model, evaluated by way of log-likelihood ratio test. Either \code{FALSE} (the default) or the column name (quoted) of the item identifier (e.g., \code{"Item"}, or \code{"Word"}).}
  \item{alpha}{Level of significance (for ANOVA and log-likelihood ratio test). Defaults to 0.05.}
  \item{llrt}{Logical. Whether log-likelihood ratio testing between the complex and simple models should be performed. Defaults to \code{FALSE}.}
  \item{prune.ranefs}{Logical. Whether to remove any random effect for which its variable is not also present in the fixed effects structure (with the exception of the grouping variables such as \code{"Subjects"} and \code{"Items"}). Defaults to \code{TRUE}. For example, if the random effects structure contains the terms \code{Condition + ROI + Group}, and the random effects structure contains the terms \code{(1 | Subject) + (0 + TrialNum | Subject)}, the ranedom effect \code{(0 + TrialNum | Subject)} will be pruned from the model given that it is not in the model's fixed effects structure.}
  \item{p.value}{Whether to use upper-bound (``upper''; the default) or lower-bound (``lower'') \emph{p}-values during backfitting.}
  \item{set.REML.FALSE}{Logical. Whether or not to set \code{REML} to \code{FALSE}. Defaults to \code{FALSE}.}
  \item{reset.REML.TRUE}{Logical. Whether or not to re-set the back-fitted model to \code{REML = TRUE}.}
  \item{log.file}{Whether a back-fitting log should be saved. By default, it is saved in a temporary folder as \code{file.path(tempdir(), paste("bfFixefLMER_F_log_", gsub(":", "-", gsub(" ", "_", date())), ".txt", sep = ""))}. The path and file name of the log can be changed to whatever the use wishes. Set to \code{FALSE} to disable.}
}
\details{The back-fitting process works as follows:
\enumerate{
\item If argument \code{llrt} is set to \code{TRUE}, \code{REML} is set to \code{FALSE};
\item First consider only highest-order interaction model terms:
  \enumerate{
       \item The model term with the highest ANOVA \emph{p}-value is identified. If this \emph{p}-value is smaller than the alpha value the model term is removed and a new model is fitted. This is repeated for each model term for which the \emph{p}-value is greater than the alpha value. The algorithm then moves on to step (b). If this \emph{p}-value is greater than the alpha value and \code{llrt} is set to \code{TRUE}, then the following is evaluated:
	\enumerate{
           \item A new model without this model term is fitted; 
	   \item The more complex and simpler models are compared by way of a log-likelihood ratio test (LLRT). If the result of the LLRT determines that the term under consideration does not increase model fit, it is removed; otherwise it is kept.
           \item Move on to the next model term with the largest \emph{p}-value greater than the alpha value and repeat steps (i)--(iii).
	}
        \item Once all highest-order interaction terms have been evaluated, go down to the the second highest order interactions: Repeat steps (ai)--(aiii) with the following addition: If a term would be removed from the model, but it is part of a high-order interaction, keep it. Once all terms of the interaction level have been evaluated, move down to the next lower-order level until main effects have been evaluated, after which the process stops.
  }
\item If argument \code{llrt} is set to \code{TRUE}, set \code{REML} to \code{TRUE} (default) unless otherwise specified.
}

In brief, if \code{llrt} is set to \code{FALSE}, a term remains in the model if its \emph{p}-value from the ANOVA is equal to or smaller than \code{alpha}; if \code{llrt} is set to \code{TRUE}, a term remains in the model if
\enumerate{
\item its \emph{p}-value from the ANOVA is equal to or smaller than \code{alpha} or, in the case that it is not, if 
\item it significantly increases model fit as determined by a log-likelihood ratio test.
}

This backfitting method was used in Newman, Tremblay, Nichols, Neville, and Ullman (accepted). If factorial terms are included in the initial model, back-fitting on F is recommended. 
}
\note{At the moment, it is absolutely necessary that names of variables not occur as parts of other variable names, e.g. "c"  and "cv". The script will recognise the presence of "c" where ever "cv" is present.}
\value{A \code{mer} model with back-fitted fixed effects (on F values) is returned and a log of the back-fitting process is printed on screen and (by default) in a log file in a temporary file.}
\references{
Newman, A.J., Tremblay, A., Nichols, E.S., Neville, H.J., and Ullman, M.T. (Resubmitted). The Influence of Language Proficiency on Lexical-Semantic Processing in Native and Late Learners of English: ERP evidence. Submitted to the \emph{Journal of Cognitive Neuroscience}.
}
\author{Antoine Tremblay, Dalhousie University, \email{trea26@gmail.com}.}

\section{Warnings}{Upper-bound \emph{p}-values can be anti-conservative, while lower-bound \emph{p}-values can be conservative. See \url{http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests} and function \code{pamer.fnc}.}

\seealso{
\code{\link{bfFixefLMER_t.fnc}; }
\code{\link{ffRanefLMER.fnc}; }
\code{\link{fitLMER.fnc}; }
\code{\link{mcposthoc.fnc}; }
\code{\link{pamer.fnc}; }
\code{\link{mcp.fnc}; }
\code{\link{romr.fnc}; }
\code{\link{perSubjectTrim.fnc} }
}
\examples{
\dontrun{
if("LCFdata"\%in\%.packages(all.available=TRUE)){
  library(LCFdata)
  # bfFixefLMER -- backfit on F
  library(lme4)
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
  
  # backfit fixed effects on F
  m2 = bfFixefLMER_F.fnc(mB, log.file = FALSE,
      llrt = FALSE)
  pamer.fnc(m2)

  plotLMER3d.fnc(m2, "LengthB", "WMC")
  devAskNewPage(ask=FALSE)
}
}
}
\keyword{models & regression}
