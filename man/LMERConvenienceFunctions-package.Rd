\name{LMERConvenienceFunctions-package}
\alias{LMERConvenienceFunctions-package}
\alias{LMERConvenienceFunctions}
\docType{package}
\title{An suite of functions to facilitate modeling with LMER (other miscellanea).}
\description{The main functions are fixed effect back-fitting functions (\code{bfFixefLMER_F.fnc} to back-fit on \emph{p}-values from ANOVA or \code{bfFixefLMER.fnc_t.fnc} to back-fit on \emph{t}-values from \code{summary(model)} and, optionally, log-likelihood ratio testing) and a random effect forward-fitting function (\code{ffRanefLMER.fnc}; log-likelihood ratio testing). There is also a function to first back-fit fixed effects from an initial model, then forward-fit random effects, and finally re-backfit fixed effects (\code{fitLMER.fnc}). Other functions include a function to compute ANOVAs with upper- or lower-bound \emph{p}-values and R-squared values for each model term; \code{pamer.fnc}), a function to graph model criticism plots (\code{mcp.fnc}), a function to trim data on model residuals (\code{romr.fnc}), one to perform per-subject trimming on the response variable (\code{perSubjectTrim.fnc}, functions to perform posthoc analyses (\code{mcposthoc.fnc}), a function to generate summaries of \code{mcposthoc} objects (\code{summary.mcposthoc}), a function to generate (dynamic) 3d plots of \code{mer} objects (\code{plotLMER3d.fnc}), a function to generate (dynamic) 3d plots of the raw data as a function of an interaction between two numeric variables (\code{plotRaw3d.fnc}), a function to plot (dynamic) 3d kernel estimates of two numeric variables (\code{plotDensity3d.fnc}), and a function to calculate the relative log-likehood between two models (\code{relLik}).. Additionally, there are functions to list files in the current directory in matrix format (\code{f}, with easily readable numbers for each file/directory), change directory (\code{cd}), and to change directory and automatically list files in new directory (\code{cdf}), to go up one directory and automatically list files (\code{cdup}), and a function to list in matrix format the column names of a data frame (\code{cn}). The data to run examples is contained in package \code{LCF_data}.
}
\details{
\tabular{ll}{
Package: \tab LMERConvenienceFunctions\cr
Type: \tab Package\cr
Version: \tab 1.7\cr
Date: \tab 2012-12-31\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{Antoine Tremblay, Dalhousie University, and Johannes Ransijn, University of Copenhagen

Maintainer: Antoine Tremblay <trea26@gmail.com>
}

\references{
Baayen, R.H. (2008). \emph{Analyzing Linguistic Data. A Practical Introduction to Statistics Using R}. Cambridge, UK: Cambridge University Press. 


Baayen, R.H., Davidson, D.J. and Bates, D.M. (2008). Mixed-effects modeling with crossed random effects for subjects and items. \emph{Journal of Memory and Language}, \emph{59}, 390--412.


Newman, A.J., Tremblay, A., Nichols, E.S., Neville, H.J., and Ullman, M.T. (2012). The Influence of Language Proficiency on Lexical Semantic Processing in Native and Late Learners of English. \emph{Journal of Cognitive Neuroscience}, \emph{25}, 1205--1223.


Newman, A.J., Tremblay, A., Neville, H.J., and Ullman, M.T. (In preparation).  The relationship between proficiency and ERP components evoked by grammatical violations in native and late learners of English.


Pinheiro, J.C. and Bates, D.M. (2000). \emph{Mixed Effects Models in S and S-Plus}. New York: Springer.


Tremblay, Antoine. (2009). \emph{Processing Advantages of Lexical Bundles: Evidence from Self-paced Reading, Word and Sentence Recall, and Free Recall with Event-related Brain Potential Recordings}. Ph.D. Dissertation. University of Alberta, Edmonton, Canada.


Tremblay, A. and Tucker B. V. (2011). The Effects of N-gram Probabilistic Measures on the Processing and Production of Four-word Sequences. \emph{The Mental Lexicon}, \emph{6(2)}, 302--324.


Tremblay, Antoine, & Newman, Aaron J. (2013). Modelling event-related potentials using linear mixed-effects regression. \emph{NeuroImage}. 55 pages.


\url{http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests}
}
\keyword{ package }
\seealso{
\code{\link{bfFixefLMER_F.fnc}};
\code{\link{bfFixefLMER_t.fnc}};
\code{\link{ffRanefLMER.fnc}};
\code{\link{fitLMER.fnc}};
\code{\link{mcposthoc.fnc}};
\code{\link{summary.mcposthoc}};
\code{\link{pamer.fnc}};
\code{\link{mcp.fnc}};
\code{\link{relLik}};
\code{\link{romr.fnc}};
\code{\link{plotLMER3d.fnc}};
\code{\link{plotDensity3d.fnc}};
\code{\link{plotRaw3d.fnc}};
\code{\link{perSubjectTrim.fnc}};
\code{\link{cn}};
\code{\link{f}};
\code{\link{cd}};
\code{\link{cdf}};
\code{\link{cdup}}
}
\examples{
\dontrun{
if(try(require(LCFdata,quietly=TRUE))){
  ### Load and format data
  require(LCFdata)
  data(eeg)

  # restrict to electrode Fz and 80--180 ms window
  eeg <- eeg[eeg$Time >= 80 & eeg$Time <= 180, ]
  eeg <- eeg[, c("Subject", "Item", "Time", "Fz",
    "FreqB", "LengthB", "WMC")]

  # mean center FreqB
  eeg$FreqBc <- eeg$FreqB - mean(eeg$FreqB)
  # split FreqBc into 3 categories. Doesn't make sense, 
  # but it's merely for example
  eeg$FreqBdc <- "high"
  eeg$FreqBdc[eeg$FreqBc<=quantile(eeg$FreqBc)[3]] <- "mid"
  eeg$FreqBdc[eeg$FreqBc<=quantile(eeg$FreqBc)[2]] <- "low"
  eeg$FreqBdc <- as.factor(eeg$FreqBdc)
  eeg$FreqBdc <- relevel(eeg$FreqBdc, "low")

  # mean center LengthB
  eeg$LengthBc <- eeg$LengthB - mean(eeg$LengthB)

  # mean center WMC
  eeg$WMCc <- eeg$WMC - mean(eeg$WMC)

  ### demonstrate plotDensity3d,fnc
  plotDensity3d.fnc(x = sort(unique(eeg$WMCc)), 
    y = sort(unique(eeg$LengthBc)))

  ### demonstrate plotDensity3d,fnc
  plotRaw3d.fnc(data = eeg, response = "Fz", pred = "WMCc",
    intr = "LengthBc", plot.type = "persp", theta = 150)

  ### Analyze data
  # fit initial model
  m0 <- lmer(Fz ~ (FreqBdc + LengthBc + WMCc)^2 + (1 | Subject), 
    data = eeg)
  m1 <- lmer(Fz ~ (FreqBdc + LengthBc + WMCc)^2 + (1 | Subject) +
    (1 | Item), data = eeg)

  # which model to choose?
  relLik(m0, m1)

  # choose m1
  # check model assumptions 
  mcp.fnc(m1)

  # remove outliers
  eeg <- romr.fnc(m1, eeg, trim = 2.5)
  eeg$n.removed
  eeg$percent.removed
  eeg<-eeg$data

  # update model
  m1 <- lmer(Fz ~ (FreqBdc + LengthBc + WMCc)^2 + (1 | Subject) +
    (1 | Item), data = eeg)
  
  # re-check model assumptions 
  mcp.fnc(m1)

  # forward-fit random effect structure (simple for the purposes
  # of the example).
  m2 <- ffRanefLMER.fnc(model = m1, data = eeg, ran.effects = 
    c("(0 + LengthBc | Subject)", "(0 + WMCc | Item)"), 
    log.file = FALSE)

  # backfit model m2. In this case, could use bfFixefLMER_t.fnc instead.
  m3 <- bfFixefLMER_F.fnc(m2, log.file = FALSE)

  # The calls to ffRanefLMER.fnc and bfFixefLMER_F.fnc could 
  # be replaced by a call to fitLMER.fnc. In this latter case, however, 
  # bfFixefLMER_F.fnc would be called first, then the random effect 
  # structure would be forward fitted, and finally teh fixed effects
  # would be backfitted again.
  m3b <- fitLMER.fnc(model = m1, ran.effects = c("(0 + LengthBc | Subject)",
    "(0 + WMCc | Item)"), backfit.on = "F", log.file = FALSE)
  pamer.fnc(m3b)
  # The results are the same. This may not necessarily be the case
  # elsewhere. First forward fitting the random effect structure and
  # then backfitting the fixed effects, potentially pruning irrelevant 
  # random effects, is probably the best approach. Nonetheless, there is 
  # no hard evidence to this effect.

  # check model assumptions 
  mcp.fnc(m3)

  # check significance of model terms
  pamer.fnc(m3)

  # Only the intercept is significant. For purposes of the 
  # example, let's perform a posthoc analysis on FreqBdc on
  # model m2.
  m2.ph <- mcposthoc.fnc(model = m2, var = list(ph1 = "FreqBdc"))


  if(try(require(languageR,quietly=TRUE))){
  # Generate plot for Length X Freq with function plotLMER.fnc
  require(languageR)
  plotLMER.fnc(m2, pred = "LengthBc", intr = list("FreqBdc", 
    levels(eeg$FreqBdc), "beg", list(1 : 3, 1 : 3)))
  }
  # effect of Length at each FreqB bin
  summary(object = m2.ph, term = "LengthBc")
  # Length effect different Freq bins? Start with low 
  # versus mid and high
  smry <- summary(object = m2.ph, term = "FreqBdclow:LengthBc")
  # mid versus low and high
  smry <- summary(object = m2.ph, term = "FreqBdcmid:LengthBc")
  }
}
}
