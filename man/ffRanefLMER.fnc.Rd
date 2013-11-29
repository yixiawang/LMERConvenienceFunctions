\name{ffRanefLMER.fnc}
\alias{ffRanefLMER.fnc}
\title{Forward-fit the random effect structure of an LMER model.}
\description{Forward-fit an LMER model's random effect structure by comparing a model without one of the specified random effects and a model with it by way of log-likelihood ratio testing. If the more complex model is a significantly better fit, the random effect is kept, otherwise it is dropped. This function can now be used with generalized linear mixed-effects models (\code{glmer}s).}
\usage{ffRanefLMER.fnc(model, ran.effects = list(ran.intercepts = 
as.character(), slopes = as.character(), corr = as.character(), 
by.vars = as.character()), alpha = 0.05, if.warn.not.add = TRUE, 
log.file = NULL)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}). This function can now be used with generalized linear mixed-effects models (\code{glmer}s).}
  \item{ran.effects}{Can be either a vector or a list. In the former case, the random effects to be evaluated are provided. For example \code{c("(1 + Frequency | Subject)", "(0 + Length | Subject)", "(1 + NSynSet | Subject)")}. In the latter case, the list can be composed of (i) a vector of random intercepts to be evaluated (\code{ran.intercepts}), (ii) a vector of random slopes to be evaluated (\code{slopes}), (iii) a vector specifying, for each element of \code{slopes}, whether the correlation between the slope and by-variables specified in \code{by.vars} should be added (\code{corr}), and (iv) a vector of ``by'' variables for the random slopes (\code{by.vars}). Values that can be supplied to the \code{corr} argument are \code{1} (add correlation), \code{0} (do not add correlation), and \code{NA} (for when the "slope" is a factor variable). Note that if a term in \code{slopes} is a factor variable, the \code{corr} value tied to it will be automatically set to \code{NA}. Also note that if no values are supplied to \code{corr}, a vector of 0 as long as the \code{slopes} vector will be automatically supplied. For example \code{list(ran.intercepts = "Word", slopes = c("Frequency", "Length", "NSynSet","Class"), corr = c(0, 0, 1, NA), by.vars = "Subject")}. Another example is \code{list(slopes = c("Trial", "Class"), by.vars = "Subject")}, where the \code{corr} argument will be equal to \code{c(0, NA)}.}
  \item{alpha}{Level of significance for log-likelihood ratio test. Defaults to 0.05.}
  \item{if.warn.not.add}{Logical. If a warning is issued after fitting a model with a new random effect (e.g., \code{false convergence} or the like), should the random effect nevertheless be evaluated? Defaults to \code{TRUE}, meaning that if such a warning is issued, the random effect will not be added to the random effects structure of the model. If set to \code{FALSE}, the random effect will be evaluated for inclusion as any other random effects would be via log likelihood ratio testing even if a warning is issued.}
  \item{log.file}{Should the back-fitting log be saved? Defaults to \code{NULL}, which means that a log file is saved in a temporary folder as \code{paste("ffRanefLMER_log_", gsub(":", "-", gsub(" ", "_", date())), ".txt", sep = "")}. The path and file name of the log can be changed to whatever the use wishes. Set to \code{FALSE} to disable.}
}
\value{A \code{mer} object with forward-fitted random effect structure as well as a log of the process is printed on screen and, optionally, printed in a log file.}
\references{Pinheiro, J.C. and Bates, D.M. (2000). \emph{Mixed Effects Models in S and S-Plus}. New York: Springer.}
\author{Antoine Tremblay, Dalhousie University, \code{trea26@gmail.com}.}

\seealso{
\code{\link{bfFixefLMER_F.fnc}; }
\code{\link{bfFixefLMER_t.fnc}; }
\code{\link{fitLMER.fnc}; }
\code{\link{mcposthoc.fnc}; }
\code{\link{pamer.fnc}; }
\code{\link{mcp.fnc}; }
\code{\link{romr.fnc}; }
\code{\link{perSubjectTrim.fnc}. }
}
\examples{
# see example in LMERConvenienceFunctions help page.
}
\keyword{models & regression}
