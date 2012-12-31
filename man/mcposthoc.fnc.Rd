\name{mcposthoc.fnc}
\alias{mcposthoc.fnc}
\title{Posthoc analyses for LMER models using parallel capabilities.}
\description{This function uses the \code{parallel} package. For each factor level, a slave process is sent to one of the computer's cores unsing function \code{mclapply} where the specified factor variables are re-leveled to each one of their levels, the \code{mer} model updated, and summaries returned. \emph{MCMC p-value calculation is now implemented}. R will wait until all slave processes have finished running. See package \code{parallel} for more information about parallel computing. Note that tradional sequential computing can be achieved by specifying \code{mc.cores = 1}. Posthoc results can be viewed with function \code{summary.mcposthoc}.}
\usage{mcposthoc.fnc(model, var, two.tailed = TRUE, 
mcmc = FALSE, nsim = 10000, ndigits = 4, mc.cores = 1, 
verbosity = 1, ...)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}) or an \code{lm} object (fitted by function \code{lm}).}
  \item{var}{A named list of variable on which to perform the posthoc analysis. For example \code{list(ph1 = c("PronomOfTheme", "AnimacyOfRec", "DefinOfRec"), ph2 = c("SemanticClass"))}.}
  \item{two.tailed}{Logical. Whether to perform one- or two-tailed \emph{t}-tests. Defaults to \code{TRUE}, i.e., two-tailed.}
  \item{mcmc}{Logical. Whether to calculate \emph{p}-values using function \code{pamer.fnc} (the default) or using function \code{pvals.fnc} from package \code{languageR}.}
  \item{nsim}{An integer denoting the required number of Markov chain Monte Carlo samples. Defaults to 10000.}
  \item{ndigits}{Integer indicating the number of decimal places to be used in the \emph{t} tables. Defaults to 4.}
  \item{mc.cores}{The number of cores to use, i.e. how many processes will be spawned (at most).}
  \item{verbosity}{Numeric. The amount of information printed to screen during the modeling process. The higher the number, the more information is printed. \code{0} turns this option off. Defaults to \code{1}.}
  \item{...}{Further arguments to pass to "mclapply".}
}
\details{If \code{var = list(ph1 = c("PronomOfTheme", "AnimacyOfRec", "DefinOfRec"))}, for example, the function will re-level and update the model on each combination of the variable levels as follows:
\preformatted{
(1) data$PronomOfTheme <- relevel(data$PronomOfTheme = "nonpronominal")
    data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "animate")
    data$DefinOfTheme <- relevel(data$DefinOfTheme = "definite")

(2) data$PronomOfTheme <- relevel(data$PronomOfTheme = "nonpronominal")
    data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "inanimate")
    data$DefinOfTheme <- relevel(data$DefinOfTheme = "definite")

(3) data$PronomOfTheme <- relevel(data$PronomOfTheme = "nonpronominal")
    data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "animate")
    data$DefinOfTheme <- relevel(data$DefinOfTheme = "indefinite")

(4) data$PronomOfTheme <- relevel(data$PronomOfTheme = "pronominal")
    data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "animate")
    data$DefinOfTheme <- relevel(data$DefinOfTheme = "definite")

(5) data$PronomOfTheme <- relevel(data$PronomOfTheme = "nonpronominal")
    data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "inanimate")
    data$DefinOfTheme <- relevel(data$DefinOfTheme = "indefinite")

(6) data$PronomOfTheme <- relevel(data$PronomOfTheme = "pronominal")
    data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "animate")
    data$DefinOfTheme <- relevel(data$DefinOfTheme = "indefinite")

(7) data$PronomOfTheme <- relevel(data$PronomOfTheme = "pronominal")
    data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "inanimate")
    data$DefinOfTheme <- relevel(data$DefinOfTheme = "indefinite")

(8) data$PronomOfTheme <- relevel(data$PronomOfTheme = "pronominal")
    data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "inanimate")
    data$DefinOfTheme <- relevel(data$DefinOfTheme = "definite")
}

On a cluster, instead of using \code{mcposthoc.fnc} it is better (faster and less complicated) to write a bash script that will create (1) \code{.R} files that will relevel the conditions and update the model, and (2) an associated \code{.sh} job submission script to submit the \code{.R} files. For example, let's consider two ERP analyses (regular past tense inflection and phrase structure) with three time windows each (300--400 ms, 550--700 ms, 750--850 in the regular past tense analysis, and 300--400 ms, 400--600 ms, and 750--850 ms in the phrase structure analysis). We investigated the effects of proficiency on ERP amplitudes. The initial models included a four-way interaction between Region of Interest (ROI) -- with levels left anterior, left central, left posterior, midline anterior, midline central, midline posterior, right anterior, right central, and right posterior) -- Group (with levels L1 and L2), Condition (wth levels control and violation), and Proficiency. After back-fitting the fixed effects, forward-fitting randomg effects, and reback-fitting the fixed effects as per \code{fitLMER.fnc}, the four-way interaction remained in every model. See Newman et al. (In preparation) for more details. The posthoc analysis script named \code{posthocs.sh} we used on the ACEnet cluster is as follows: 
\preformatted{
      time=(Reg300400 Reg550700 Reg750850 PS300400 PS400600 PS750850)
      condition=(Good Bad)
      group=(L1 L2)
      roi=(Lant Lcent Lpost Mant Mcent Mpost Rant Rcent Rpost)
      
      for t in ${time[*]}; do for i in ${condition[*]}; do for j in ${group[*]}; do for k in ${roi[*]}; do 
          ### create .R file where the modell is updated on the data where 
          ### re-leveld on each possible combination of variable levels
          export CONDITION=$i;
          export GROUP=$j;
          export ROI=$k; 
          echo 'condition<-Sys.getenv("CONDITION")' > "ph"$t$CONDITION$GROUP$ROI".R"
          echo 'group<-Sys.getenv("GROUP")' >> "ph"$t$CONDITION$GROUP$ROI".R"
          echo 'roi<-Sys.getenv("ROI")' >> "ph"$t$CONDITION$GROUP$ROI".R"
          echo 'load("models/m1'$t'.rda")' >> "ph"$t$CONDITION$GROUP$ROI".R"
          echo 'dat<-m1@frame' >> "ph"$t$CONDITION$GROUP$ROI".R"
          echo 'dat$Condition<-relevel(dat$Condition,'condition')' >> "ph"$t$CONDITION$GROUP$ROI".R"
          echo 'dat$Group<-relevel(dat$Group,'group')' >> "ph"$t$CONDITION$GROUP$ROI".R"
          echo 'dat$ROI<-relevel(dat$ROI,'roi')' >> "ph"$t$CONDITION$GROUP$ROI".R"
          echo 'm1<-update(m1,.~.,data=dat)' >> "ph"$t$CONDITION$GROUP$ROI".R"
          echo 'save(m1,file="ph'$t$CONDITION$GROUP$ROI'.rda")' >> "ph"$t$CONDITION$GROUP$ROI".R"
      
          ### create the job submission script for the .R file created above
          echo '#$ -S /bin/bash' > "job.ph"$t$CONDITION$GROUP$ROI".sh"
          echo '#$ -cwd' >> "job.ph"$t$CONDITION$GROUP$ROI".sh"
          echo '#$ -j y' >> "job.ph"$t$CONDITION$GROUP$ROI".sh"
          echo '#$ -l h_rt=48:00:00' >> "job.ph"$t$CONDITION$GROUP$ROI".sh"
          echo '#$ -l h_vmem=8G' >> "job.ph"$t$CONDITION$GROUP$ROI".sh"
          echo '#$ -R y' >> "job.ph"$t$CONDITION$GROUP$ROI".sh"
          echo '#$ -N "ph'$t$CONDITION$GROUP$ROI'"' >> "job.ph"$t$CONDITION$GROUP$ROI".sh"
          echo 'R -q -f ph'$t$CONDITION$GROUP$ROI'.R' >> "job.ph"$t$CONDITION$GROUP$ROI".sh"
      
          ### submit the job
          qsub  "job.ph"$t$CONDITION$GROUP$ROI".sh"
      done; done; done; done

}
and then type in the console
\preformatted{
      . posthocs.sh

}
On the ACEnet cluster, this results in 2 * 3 * 9 * 2 * 2 = 216 independent analyses, simultaneously using a total of 216 cores and 1728 GB of RAM. This posthoc analysis completes in about 3-6 hours.
}
\section{warning}{Parallel computing capabilities will not be available on Windows because \code{mclapply} relies on forking. Sequential computing, however, will work on Windows if \code{mc.cores = 1} (the default).}
\value{
An object of class "mcposthoc" with the following slots:
\item{n }{The number of data points in data frame \code{data}.}        
\item{var }{A named list containing the names of the variables used in the posthoc.}
\item{summaries }{A named list containing the posthoc summaries for each factor re-leveling. If \code{mcmc = FALSE}, data frames with upper- and lower-bound (anti-conservative and conservative, respectively) \emph{df}s, \emph{p}-values, and deviance explained (\%) for each model term. If \code{mcmc = TRUE}, data frames with the estimated coefficients, their MCMC mean, the HPD 95 and the probability based on the \emph{t} distribution with the number of observations minus the number of fixed-effects coefficients as degrees of freedom.  This last \emph{p}-value is anti-conservative, especially for small data sets.}
}
\author{Antoine Tremblay, Dalhousie University, \email{trea26@gmail.com}.}
\references{
Newman, A.J., Tremblay, A., Neville, H.J., and Ullman, M.T. (In preparation).  The relationship between proficiency and ERP components evoked by grammatical violations in native and late Learners of English.
}
\seealso{
\code{\link{summary.mcposthoc}}
}
\examples{
# see example in LMERConvenienceFunctions help page.
}
\keyword{models & regression}
