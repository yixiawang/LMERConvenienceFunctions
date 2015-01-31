\name{LMERConvenienceFunctions-package}
\alias{LMERConvenienceFunctions-package}
\alias{LMERConvenienceFunctions}
\docType{package}
\title{Model Selection and Post-hoc Analysis for (G)LMER Models}
\description{The main function of the package is to perform backward selection of fixed effects, forward fitting of the random effects, and post-hoc analyses using parallel capabilities. Other functionality includes the computation of ANOVAs with upper- or lower-bound p-values and R-squared values for each model term, model criticism plots, data trimming on model residuals, and data visualization. The data to run examples is contained in package \code{LCF_data}.}
\details{
\tabular{ll}{
Package: \tab LMERConvenienceFunctions\cr
Type: \tab Package\cr
Version: \tab 2.10\cr
Date: \tab 2015-1-31\cr
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


Quene, H., & van den Bergh, H. (2008). Examples of mixed-effects modeling with crossed random effects and with binomial data. \emph{Journal of Memory and Language}, 59, 413--425. doi: 10.1016/j.jml.2008.02.002.


Symonds, M.R.E and Moussalli, A. (2011). \emph{A brief guide to model selection, multimodel inference and model averaging in behavioural ecology using Akaike's information criterion.} \emph{Behavioral Ecology and Sociobiology}, 65, 13--21. doi: 10.1007/s00265-010-1037-6


Tremblay, Antoine. (2009). \emph{Processing Advantages of Lexical Bundles: Evidence from Self-paced Reading, Word and Sentence Recall, and Free Recall with Event-related Brain Potential Recordings}. Ph.D. Dissertation. University of Alberta, Edmonton, Canada.


Tremblay, A. and Tucker B. V. (2011). The Effects of N-gram Probabilistic Measures on the Processing and Production of Four-word Sequences. \emph{The Mental Lexicon}, \emph{6(2)}, 302--324.
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
\code{\link{plotLMER.fnc}};
\code{\link{plotLMER3d.fnc}};
\code{\link{plotDensity3d.fnc}};
\code{\link{plotRaw3d.fnc}};
\code{\link{perSubjectTrim.fnc}};
\code{\link{cn}};
\code{\link{f}};
\code{\link{cd}};
\code{\link{cdf}};
\code{\link{cdup}}.
}
\examples{
\dontrun{
  ############################################
  #            Load and format data.         #
  ############################################
  library(LCFdata)
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

  ############################################
  #      Demonstrate plotDensity3d.fnc.      #
  ############################################
  plotDensity3d.fnc(x = sort(unique(eeg$WMCc)), 
    y = sort(unique(eeg$LengthBc)))

  ############################################
  #        Demonstrate plotRaw3d.fnc.        #
  ############################################
  plotRaw3d.fnc(data = eeg, response = "Fz", pred = "WMCc",
    intr = "LengthBc", plot.type = "persp", theta = 150)

  ############################################
  #       Analyze data. Demonstrate model    #
  #       selection, and diagnostic plots.   #
  #       Also demonstrate forward fitting   #
  #       of random effects and back fitting #
  #       of fixed effects. Finally,         #
  #       demonstrate pamer.fnc.             # 
  ############################################
  library(lme4)
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
  m2 <- ffRanefLMER.fnc(model = m1, ran.effects = 
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

  ############################################
  #       Demonstrate mcposthoc.fnc and      #
  #       summary.mcposthoc.                 #
  ############################################
  # Only the intercept is significant. For purposes of the 
  # example, let's perform a posthoc analysis on FreqBdc on
  # model m2.
  m2.ph <- mcposthoc.fnc(model = m2, var = list(ph1 = "FreqBdc"))

  # Now check if and how the different levels differ between
  # each other. First check high vs mid and high vs low:
  summary(m2.ph, term = "FreqBdchigh") 
  # Then low vs mid (the low vs high row is redundant from the 
  # above summary):
  summary(m2.ph, term = "FreqBdclow")
  # Note that none of the levels differ from each other. Indeed, 
  # the backfitting process indicated that the model only has an 
  # intercept (i.e., the FreqBc factor variable was not significant).

  # Just to show how one would look at posthocs for interactions. Let's 
  # look at the effect of Length at each FreqB bin:
  summary(object = m2.ph, term = "LengthBc")
  # Does Length effect different Freq bins? Start with low 
  # versus mid and high
  smry <- summary(object = m2.ph, term = "FreqBdclow:LengthBc")
  # then mid versus low and high
  smry <- summary(object = m2.ph, term = "FreqBdcmid:LengthBc")

  ############################################
  #       Demonstrate `revived' version of   #
  #       plotLMER.fnc and plotLMER3d.fnc.   #
  ############################################
  # Generate plot for Length X Freq with function plotLMER.fnc.
  plotLMER.fnc(m2, pred = "LengthBc", intr = list("FreqBdc", 
    levels(eeg$FreqBdc), "beg", list(1 : 3, 1 : 3)))

  # Plotting the Length:WMC interaction with plotLMER3d.fnc. It'll
  # take a little bit of time.
  plotLMER3d.fnc(m2,"LengthBc","WMCc")
  # Plot it a second time to demonstrate caching. You can notice the 
  # speed-up.
  plotLMER3d.fnc(m2,"LengthBc","WMCc")


  ############################################
  #       Demonstrate modeling and           #
  #       backfitting of glmer.              #
  ############################################
  # Split FreqBc into 2 categories.
  eeg$FreqBdc <- "high"
  eeg$FreqBdc[eeg$FreqBc<=median(eeg$FreqBc)] <- "low"
  eeg$FreqBdc <- as.factor(eeg$FreqBdc)
  eeg$FreqBdc <- relevel(eeg$FreqBdc, "low")

  # Fit glmer model.
  m4 <- glmer(FreqBdc ~ (Fz + LengthBc + WMCc)^2 + (1 | Subject),
	family = "binomial", data = eeg)
  summary(m4)

  # Back fit fixed effects, forward fit random effects, and then
  # re-back fit fixed effects. Need to set argument backfit.on to "t".
  m5 <- fitLMER.fnc(model = m4, ran.effects = "(0 + LengthBc | Subject)",
	backfit.on = "t", log.file = FALSE)
  summary(m5)

  # Plot the 2-way interaction.
  plotLMER.fnc(m5, pred = "Fz", intr = list("LengthBc", 
	quantile(eeg$LengthBc), "med",list(1:5,1:5)))

  # Look at the same plot, but in 3d.
  plotLMER3d.fnc(m5, pred = "Fz", intr = "LengthBc")

  ############################################
  #       Test backfitting on AIC,           #
  #       BIC, llrt, relLik.AIC, and         #
  #       relLik.BIC.                        #
  ############################################
  # AIC
  m.test <- bfFixefLMER_F.fnc(m2, method = "AIC",
	log.file = FALSE)
  m.test <- bfFixefLMER_t.fnc(m2, method = "AIC",
	log.file = FALSE)
  m.test <- bfFixefLMER_t.fnc(m4, method = "AIC",
	log.file = FALSE)
  m.test <- bfFixefLMER_F.fnc(m4, method = "AIC",
	log.file = FALSE)

  # BIC
  m.test <- bfFixefLMER_F.fnc(m2, method = "BIC",
	log.file = FALSE)
  m.test <- bfFixefLMER_t.fnc(m2, method = "BIC",
	log.file = FALSE)
  m.test <- bfFixefLMER_t.fnc(m4, method = "BIC",
	log.file = FALSE)

  # llrt
  m.test <- bfFixefLMER_F.fnc(m2, method = "llrt",
	log.file = FALSE)
  m.test <- bfFixefLMER_t.fnc(m2, method = "llrt",
	log.file = FALSE)
  m.test <- bfFixefLMER_t.fnc(m4, method = "llrt",
	log.file = FALSE)

  # relLik.AIC
  m.test <- bfFixefLMER_F.fnc(m2, method = "relLik.AIC",
	log.file = FALSE)
  m.test <- bfFixefLMER_t.fnc(m2, method = "relLik.AIC",
	log.file = FALSE)
  m.test <- bfFixefLMER_t.fnc(m4, method = "relLik.AIC",
	log.file = FALSE)

  # relLik.BIC
  m.test <- bfFixefLMER_F.fnc(m2, method = "relLik.BIC",
	log.file = FALSE)
  m.test <- bfFixefLMER_t.fnc(m2, method = "relLik.BIC",
	log.file = FALSE)
  m.test <- bfFixefLMER_t.fnc(m4, method = "relLik.BIC",
	log.file = FALSE)
}
}
