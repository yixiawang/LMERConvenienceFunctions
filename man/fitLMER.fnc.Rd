\name{fitLMER.fnc}
\alias{fitLMER.fnc}
\title{Back-fit fixed effects and forward-fit random effects of an LMER model.}
\description{The function follows these steps: (1) If \code{llrt} is set to \code{TRUE}, set \code{REML} to \code{FALSE} (unless specified otherwise); (2) back-fit initial model either on \emph{F}- (by default) or on \emph{t}-values; (3) forward-fit random effects; (4) re-back-fit fixed effects; (5) if \code{llrt} is set to \code{TRUE}, set \code{REML} to \code{TRUE} (unless specified otherwise). Note that, this function CAN be used with generalized linear mixed-effects models (\code{glmer}s).}
\usage{fitLMER.fnc(model, item = FALSE, backfit.on = c("F",
"t"), method = c("F", "t", "llrt", "AIC", "BIC", "relLik.AIC", 
"relLik.BIC"), threshold = NULL, t.threshold = NULL, 
ran.effects = list(ran.intercepts = as.character(), 
slopes = as.character(), corr = as.character(), 
by.vars = as.character()), alpha = NULL, alphaitem = NULL, 
if.warn.not.add = TRUE, prune.ranefs = TRUE, p.value = "upper", 
set.REML.FALSE = TRUE, keep.single.factors = FALSE, 
reset.REML.TRUE = TRUE, log.file.name = NULL)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}). This function can be used with generalized linear mixed-effects models (\code{glmer}s) if argument \code{backfit.on} is set to "t", but not if it is set to "F".}
  \item{item}{Whether or not to evaluate the addition of by-item random intercepts to the model, evaluated by way of log-likelihood ratio test. Either \code{FALSE} (the default, does not evaluate this addition) or the column name (quoted) of the item identifier (e.g., \code{"Item"}, \code{"Word"}).}
  \item{backfit.on}{Either "F" (default) or "t". Refers to the statistic which will be used to determine which term to test and potentially remove from the model. If you are backfitting a generalized linear mixed-effects model (\code{glmer}), make sure to set \code{backfit.on} to "t".}
  \item{method}{Backfitting method. One of "F" (\emph{p}-value), "t" (\emph{t} statistic), "llrt" (log-likelihood ratio test), "AIC", "BIC", "relLik.AIC", or "relLik.BIC" (relative likelihood, see function \code{relLik}). Defaults to \code{F}. You can find information regarding differences between AIC and BIC from \code{http://methodology.psu.edu/eresources/ask/sp07}.}
  \item{threshold}{Method-specific threshold for parameter selection. It refers to alpha in the case of "F" and "llrt", to the \emph{t}-value in case of "t", to the minimum reduction in likelihood in the case of "AIC" and "BIC", or to the minimum difference in probability in the case of "relLik.AIC" and "relLik.BIC". Defaults \code{NULL}, which means \code{0.05} for "F" and "llrt", 2 for "t", 5 for "AIC" and "BIC", and 4 for "relLik.AIC" and "relLik.BIC".}
    \item{t.threshold}{Defaults to \code{NULL}. If the \code{method = "t"}, it is the \emph{t}-value below which a model term is dropped (if \code{t.threshold = NULL}, it will be set to 2). Otherwise it is the threshold for \emph{t}-value below which a test (see \code{method}) is performed between a model with the term under consideration and a simpler model without it (if \code{t.threshold = NULL}, it is set to \code{Inf}, which means that all terms are tested.} 
  \item{ran.effects}{Can be either a vector or a list. In the former case, the random effects to be evaluated are provided. For example \code{c("(1 + Frequency | Subject)", "(0 + Length | Subject)", "(1 + NSynSet | Subject)")}. In the latter case, the list can be composed of (i) a vector of random intercepts to be evaluated (\code{ran.intercepts}), (ii) a vector of random slopes to be evaluated (\code{slopes}), (iii) a vector specifying, for each element of \code{slopes}, whether the correlation between the slope and by-variables specified in \code{by.vars} should be added (\code{corr}), and (iv) a vector of ``by'' variables for the random slopes (\code{by.vars}). Values that can be supplied to the \code{corr} argument are \code{1} (add correlation), \code{0} (do not add correlation), and \code{NA} (for when the "slope" is a factor variable). Note that if a term in \code{slopes} is a factor variable, the \code{corr} value tied to it will be automatically set to \code{NA}. Also note that if no values are supplied to \code{corr}, a vector of 0 as long as the \code{slopes} vector will be automatically supplied. For example \code{list(ran.intercepts = "Word", slopes = c("Frequency", "Length", "NSynSet","Class"), corr = c(0, 0, 1, NA), by.vars = "Subject")}. Another example is \code{list(slopes = c("Trial", "Class"), by.vars = "Subject")}, where the \code{corr} argument will be equal to \code{c(0, NA)}.}
  \item{alpha}{If the method is \code{F}, it is the \emph{p}-value (from \code{pamer.fnc}) above which a model term is dropped. In this case, it defaults to the value passed to argument \code{threshold}, i.e., 0.05. Otherwise it is the \emph{p}-value threshold above which a test (see \code{method}) is performed between a model with the term under consideration and a simpler model without it (in this case, defaults to \code{0}, i.e. all terms will be tested).}
  \item{alphaitem}{Alpha value for the evaluation of by-item random intercepts. Defaults to \code{0.05} or to the specified threshold.}
  \item{if.warn.not.add}{Logical. If a warning is issued after fitting a model with a new random effect (e.g., \code{false convergence} or the like), should the random effect nevertheless be evaluated? Defaults to \code{TRUE}, meaning that if such a warning is issued, the random effect will not be added to the random effects structure of the model. If set to \code{FALSE}, the random effect will be evaluated for inclusion as any other random effects would be via log likelihood ratio testing even if a warning is issued.}
  \item{prune.ranefs}{Logical. Whether to remove any random effect for which its variable is not also present in the fixed effects structure (with the exception of the grouping variables such as \code{"Subjects"} and \code{"Items"}). Defaults to \code{TRUE}. For example, if the random effects structure contains the terms \code{Condition + ROI + Group}, and the random effects structure contains the terms \code{(1 | Subject) + (0 + TrialNum | Subject)}, the ranedom effect \code{(0 + TrialNum | Subject)} will be pruned from the model given that it is not in the model's fixed effects structure.}
  \item{p.value}{Whether to use upper-bound (``upper''; the default) or lower-bound (``lower'') \emph{p}-values when back-fitting with method "F".}
  \item{set.REML.FALSE}{Logical. Whether or not to set \code{REML} to \code{FALSE}. Defaults to \code{FALSE}.}
  \item{reset.REML.TRUE}{Logical. Whether or not to re-set the back-fitted model to \code{REML = TRUE}.}
  \item{keep.single.factors}{Logical. Whether or not main effects are kept (not subjected to testing and reduction). Defaults to \code{FALSE}.}     
  \item{log.file.name}{Should the back-fitting log be saved? Defaults to \code{NULL}, which means that a log file is saved in a temporary folder (platform dependent) as \code{file.path(tempdir(), paste("fitLMER_log_", gsub(":", "-", gsub(" ", "_", date())), ".txt", sep = ""))}. The path and file name of the log can be changed to whatever the use wishes. Set to \code{FALSE} to disable.}
}
\details{The process has three stages. In the first stage, either \code{bfFixefLMER_F.fnc} or \code{bfFixefLMER_t.fnc} is called (depending on the user's choice) and the fixed effects are back-fitted accordingly. In the second stage, \code{ffRanefLMER.fnc} is called and random effects are forward-fitted. In the third stage, the fixed effects are back-fitted again. This is done because the inclusion of certain random effects sometimes renders certain fixed effects non-significant. This process was used in Tremblay and Tucker (2011) and in Newman, Tremblay, Nichols, Neville, and Ullman (2012).

If, for example, you have many analyses to run and a cluster is available, write a bash script that will create (1) \code{.R} files that will relevel the conditions and update the model, and (2) an associated \code{.sh} job submission script to submit the \code{.R} files. For example, let's consider two ERP analyses all in a time window ranging from 100 to 250 ms. Two three-way interactions were considered: Position (factor; 1 to 6) X  Length of the second word of a four-word sequence (e.g., \emph{in the middle of}) X Working Memory Capacity score (continuous, from 0 to 100) and Trial (continuous; 1 to 432) X Length X Working Memory Capacity. Analyses were performed at electrodes Fp1 Fp2 AF3 AF4 F7 F3 Fz F4 F8 FC5 FC1 FC2 FC6 T7 C3 Cz C4 T8 CP5 CP1 CP2 CP6. See Tremblay and Newman (In preparation) for more details. The analysis script named \code{Fp1-CP6_100250.sh} we used on the ACEnet cluster is as follows: 
\preformatted{
      electrodes=(Fp1 Fp2 AF3 AF4 F7 F3 Fz F4 F8 FC5 FC1 FC2 FC6 T7 C3 Cz C4 T8 CP5 CP1 CP2 CP6)
      for e in ${electrodes[*]}; do 
            export E=$e;
            # create .R script to load data, perform necessary manipulations
            # and perform the analysis using fitLMER.fnc
            echo 'e<-Sys.getenv("E")' > $e".R"
            echo 'load("../data/eeg600_trim_v2.rda")' >> $e".R"
            echo 'dat0<-dat' >> $e".R"
            echo 'rm(dat);gc(T,T)' >> $e".R"
            echo 'dat <- dat0[dat0$Time >= 100 & dat0$Time <= 250, , drop = TRUE]' >> $e".R"
            echo 'dat <- dat[dat$Electrode == e, , drop = TRUE]' >> $e".R"
            echo 'subj<-sort(unique(dat$Subject))' >> $e".R"
            echo 'for(i in subj){' >> $e".R"
            echo 'tmp<-dat[dat$Subject==i,,drop=TRUE]' >> $e".R"
            echo 'tmp$newfact<-paste(tmp$Block,tmp$Position,sep="_")' >> $e".R"
            echo 'newvec<-vector("numeric")' >> $e".R"
            echo 'for(j in 1:length(unique(tmp$newfact))){' >> $e".R"
            echo 'newvec<-c(newvec,rep(j,nrow(tmp[tmp$newfact==unique(tmp$newfact)[j],])))' >> $e".R"
            echo '}' >> $e".R"
            echo 'tmp$Trial<-newvec' >> $e".R"
            echo 'if(grep(i,subj)[1]==1){' >> $e".R"
            echo 'newdat<-tmp' >> $e".R"
            echo '}else{' >> $e".R"
            echo 'newdat<-rbind(newdat,tmp)' >> $e".R"
            echo '}' >> $e".R"
            echo '}' >> $e".R"
            echo 'dat<-newdat' >> $e".R"
            echo 'dat$Position<-as.factor(dat$Position)' >> $e".R"
            echo 'm7 <- lmer(Amplitude ~ (Position + Trial)*(LengthBc * WMCc) + ' >> $e".R"
            echo '(1 | Subject), data = dat)' >> $e".R"
            echo 'm7b<-fitLMER.fnc(m7,item="Item",ran.effects=c("(0+Trial|Subject)",' >> $e".R"
            echo '"(0+LengthBc|Subject)","(0+Trial|Item)","(0+WMCc|Item)",' >> $e".R"
            echo '"(Position|Subject)"))' >> $e".R"
            echo 'smry<-pamer.fnc(m7b)' >> $e".R"
            echo 'save(m7b,file=file.path("..","models",paste("m7b_",e,"_100250.rda",sep="")))' >> $e".R"
            echo 'save(smry,file=file.path("..","summaries",paste("smry_m7b_",e,"_100250.rda",sep="")))' >> $e".R"
  
            ### create the job submission script for the .R file created above
            echo '#$ -S /bin/bash' > "job."$e".sh"
            echo '#$ -cwd' >> "job."$e".sh"
            echo '#$ -j y' >> "job."$e".sh"
            echo '#$ -l h_rt=48:00:00' >> "job."$e".sh"
            echo '#$ -l h_vmem=8G' >> "job."$e".sh"
            echo '#$ -R y' >> "job."$e".sh"
            echo '#$ -N '$e >> "job."$e".sh"
            echo 'R -q -f '$e'.R' >> "job."$e".sh"
            
            ### submit the job
            qsub  "job."$e".sh"
      done;

}
and then type in the console
\preformatted{
      . Fp1-CP6_100250.sh

}
On the ACEnet cluster, this results in 22 independent analyses, simultaneously using a total of 22 cores and 176 GB of RAM. This analysis completes in about 30 minutes to 1 hour.
}
\value{A \code{mer} object with back-fitted fixed effects and forward-fitted random effects, as well as a log of the process, which is printed on screen and, optionally, printed in a log file.}
\references{
Baayen, R.H., Davidson, D.J. and Bates, D.M. (2008). Mixed-effects modeling with crossed random effects for subjects and items. Journal of Memory and Language, 59, 390--412.

Newman, A.J., Tremblay, A., Nichols, E.S., Neville, H.J., and Ullman, M.T. (2012). The Influence of Language Proficiency on Lexical Semantic Processing in Native and Late Learners of English. \emph{Journal of Cognitive Neuroscience}, \emph{25}, 1205--1223.

Pinheiro, J.C. and Bates, D.M. (2000). \emph{Mixed Effects Models in S and S-Plus}. New York: Springer.

Tremblay, A. and Tucker B. V. (2011). The Effects of N-gram Probabilistic Measures on the Processing and Production of Four-word Sequences. \emph{The Mental Lexicon}, \emph{6(2)}, 302--324.
}
\author{Antoine Tremblay, Dalhousie University, \email{trea26@gmail.com}}

\section{Warnings}{Upper-bound \emph{p}-values can be anti-conservative, while lower-bound \emph{p}-values can be conservative. See \url{http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests} and function \code{pamer.fnc}.}

\seealso{
\code{\link{bfFixefLMER_F.fnc}; }
\code{\link{bfFixefLMER_t.fnc}; }
\code{\link{ffRanefLMER.fnc}; }
\code{\link{mcposthoc.fnc}; }
\code{\link{pamer.fnc}; }
\code{\link{mcp.fnc}; }
\code{\link{relLik}; }
\code{\link{romr.fnc}; }
\code{\link{perSubjectTrim.fnc}. }
}
\examples{
# see example LMERConvenienceFunctions help page.
}
\keyword{models & regression}
