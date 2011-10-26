\name{LMERConvenienceFunctions-package}
\alias{LMERConvenienceFunctions-package}
\alias{LMERConvenienceFunctions}
\docType{package}
\title{An suite of functions to facilitate modeling with LMER (other miscellanea).}
\description{The main functions are fixed effect back-fitting functions (\code{bfFixefLMER_F.fnc} to back-fit on \emph{p}-values from ANOVA or \code{bfFixefLMER.fnc_t.fnc} to back-fit on \emph{t}-values from \code{summary(model)} and, optionally, log-likelihood ratio testing) and a random effect forward-fitting function (\code{ffRanefLMER.fnc}; log-likelihood ratio testing). There is also a function to first back-fit fixed effects from an initial model, then forward-fit random effects, and finally re-backfit fixed effects (\code{fitLMER.fnc}). Other functions include a function to compute ANOVAs with upper- or lower-bound \emph{p}-values and R-squared values for each model term; \code{pamer.fnc}), a function to graph model criticism plots (\code{mcp.fnc}), a function to trim data on model residuals (\code{romr.fnc}), one to perform per-subject trimming on the response variable (\code{perSubjectTrim.fnc}, functions to perform posthoc analyses (\code{posthoc.fnc} and \code{mcposthoc.fnc}), a function to generate (dynamic) 3d plots of \code{mer} objects (\code{plotLMER3d.fnc}), a function to generate (dynamic) 3d plots of the raw data as a function of an interaction between two numeric variables (\code{plotRaw3d.fnc}), and a function to plot (dynamic) 3d kernel estimates of two numeric variables (\code{plotDensity3d.fnc}). Additionally, there are functions to list files in the current directory in matrix format (\code{f}, with easily readable numbers for each file/directory), change directory (\code{cd}), and to change directory and automatically list files in new directory (\code{cdf}), to go up one directory and automatically list files (\code{cdup}), and a function to list in matrix format the column names of a data frame (\code{cn}).
}
\details{
\tabular{ll}{
Package: \tab LMERConvenienceFunctions\cr
Type: \tab Package\cr
Version: \tab 1.6.5\cr
Date: \tab 26-10-2011\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{Antoine Tremblay, IWK Health Center 

Maintainer: Antoine Tremblay <trea26@gmail.com>
}

\references{
Baayen, R.H. (2008). \emph{Analyzing Linguistic Data. A Practical Introduction to Statistics Using R}. Cambridge, UK: Cambridge University Press. 

Baayen, R.H., Davidson, D.J. and Bates, D.M. (2008). Mixed-effects modeling with crossed random effects for subjects and items. \emph{Journal of Memory and Language}, \emph{59}, 390--412.

Newman, A.J., Tremblay, A., Nichols, E.S., Neville, H.J., and Ullman, M.T. (Accepted). The Influence of Language Proficiency on Lexical-Semantic Processing in Native and Late Learners of English: ERP evidence. Submitted to the \emph{Journal of Cognitive Neuroscience}.

Pinheiro, J.C. and Bates, D.M. (2000). \emph{Mixed Effects Models in S and S-Plus}. New York: Springer.

Tremblay, Antoine. (2009). \emph{Processing Advantages of Lexical Bundles: Evidence from Self-paced Reading, Word and Sentence Recall, and Free Recall with Event-related Brain Potential Recordings}. Ph.D. Dissertation. University of Alberta, Edmonton, Canada. Available for download at \url{http://www.ualberta.ca/~antoinet/ThesisDraft_10_B.pdf}.

Tremblay, A. and Tucker B. V. (2011). The Effects of N-gram Probabilistic Measures on the Processing and Production of Four-word Sequences. \emph{The Mental Lexicon}, \emph{6(2)}, 302--324.

Tremblay, Antoine, and Newman, Aaron J. (In Preparation). \emph{The Analysis of Event-related Potentials using Linear Mixed-effects Models with Complex Random-effect Structures}.

\url{http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests}
}
\keyword{ package }
\seealso{
\code{\link{bfFixefLMER_F.fnc}};
\code{\link{bfFixefLMER_t.fnc}};
\code{\link{ffRanefLMER.fnc}};
\code{\link{fitLMER.fnc}};
\code{\link{posthoc.fnc}};
\code{\link{mcposthoc.fnc}};
\code{\link{pamer.fnc}};
\code{\link{mcp.fnc}};
\code{\link{romr.fnc}};
\code{\link{plotLMER3d.fnc}};
\code{\link{plotDensity3d.fnc}};
\code{\link{plotRaw3d.fnc}};
\code{\link{perSubjectTrim.fnc}};
\code{\link{cn}};
\code{\link{f}};
\code{\link{cd}};
\code{\link{cdf}};
\code{\link{cdup}};
\code{\link{erpFz}};
\code{\link{z}}
}
