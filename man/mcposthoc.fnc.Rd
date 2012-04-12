\name{mcposthoc.fnc}
\alias{mcposthoc.fnc}
\title{Posthoc analyses for LMER models using multicore capabilities.}
\description{This function uses the \code{multicore} package. For each factor level, a slave process is sent to one of the computer's cores unsing function \code{mclapply} where the specified factor variables are re-leveled to each one of their levels, the \code{mer} model updated, and summaries returned. \emph{MCMC p-value calculation is now implemented}. R will wait until all slave processes have finished running. See package \code{multicore} for more information about multicore computing. Note that tradional sequential computing can be achieved by specifying \code{mc.cores = 1}. Posthoc results can be viewed with function \code{summary.mcposthoc}.}
\usage{mcposthoc.fnc(model, var, two.tailed = TRUE, 
mcmc = FALSE, nsim = 10000, ndigits = 4, mc.cores = 1, ...)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}) or an \code{lm} object (fitted by function \code{lm}).}
  \item{var}{A named list of variable on which to perform the posthoc analysis. For example \code{list(ph1 = c("PronomOfTheme", "AnimacyOfRec", "DefinOfRec"), ph2 = c("SemanticClass"))}.}
  \item{two.tailed}{Logical. Whether to perform one- or two-tailed \emph{t}-tests. Defaults to \code{TRUE}, i.e., two-tailed.}
  \item{mcmc}{Logical. Whether to calculate \emph{p}-values using function \code{pamer.fnc} (the default) or using function \code{pvals.fnc} from package \code{languageR}.}
  \item{nsim}{An integer denoting the required number of Markov chain Monte Carlo samples. Defaults to 10000.}
  \item{ndigits}{Integer indicating the number of decimal places to be used in the \emph{t} tables. Defaults to 4.}
  \item{mc.cores}{The number of cores to use, i.e. how many processes will be spawned (at most).}
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
}
\value{
An object of class "mcposthoc" with the following slots:
\item{n }{The number of data points in data frame \code{data}.}        
\item{var }{A named list containing the names of the variables used in the posthoc.}
\item{summaries }{A named list containing the posthoc summaries for each factor re-leveling. If \code{mcmc = FALSE}, data frames with upper- and lower-bound (anti-conservative and conservative, respectively) \emph{df}s, \emph{p}-values, and deviance explained (\%) for each model term. If \code{mcmc = TRUE}, data frames with the estimated coefficients, their MCMC mean, the HPD 95 and the probability based on the \emph{t} distribution with the number of observations minus the number of fixed-effects coefficients as degrees of freedom.  This last \emph{p}-value is anti-conservative, especially for small data sets.}
}
\author{Antoine Tremblay, Dalhousie University, \email{trea26@gmail.com}.}
\seealso{
\code{\link{bfFixefLMER_F.fnc}; }
\code{\link{bfFixefLMER_t.fnc}; }
\code{\link{ffRanefLMER.fnc}; }
\code{\link{fitLMER.fnc}; }
\code{\link{pamer.fnc}; }
\code{\link{mcp.fnc}; }
\code{\link{romr.fnc}; }
\code{\link{perSubjectTrim.fnc}. }
}
\examples{
\dontrun{
if("multicore" \%in\% .packages(all.available = TRUE)){
  library(multicore)
  if("languageR" \%in\% .packages(all.available = TRUE)){
    if("LCFdata"\%in\% .packages(all.available=TRUE)){
      library(LCFdata)
      library(lme4)
		  data(eeg)
		  eeg <- eeg[eeg$Time >= 100 & eeg$Time <= 175,,drop = TRUE]
		  eeg$FreqBd <- as.factor(cut(eeg$FreqB, 3, 
		  c("low", "med", "high")))
		  eeg$LengthBd <- as.factor(cut(eeg$LengthB, 3, 
		  c("short", "med", "long")))
		
		  # fit initial model
		  m1 <- lmer(Fz ~ WMC * (FreqBd + LengthBd) + (1 | Subject) + 
		      (1 | Item), data = eeg)

      pamer.fnc(m1)
		
      # setting nsim to 100 simply to make the example run time
      # short. In other cases, nsim should at least be 10000, the
      # default in pvals.fnc().
		  m1.posthoc <- mcposthoc.fnc(m1, var = list(ph1 = "FreqBd", 
            ph2 = "LengthBd"), mcmc = TRUE, nsim = 100, 
            mc.cores = 2)
		
		  names(m1.posthoc)
		  as.matrix(names(m1.posthoc$var))
		  as.matrix(names(m1.posthoc$summaries))
    }
  }
}
}
}
\keyword{models & regression}
