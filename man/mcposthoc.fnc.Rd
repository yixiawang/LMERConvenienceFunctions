\name{mcposthoc.fnc}
\alias{mcposthoc.fnc}
\title{Posthoc analyses for LMER models using multicore capabilities.}
\description{This function is a essentially the same as \code{posthoc.fnc} with hte exception that it uses the \code{multicore} package. For each factor level, a slave process is sent to one of the computer's cores where the specified factor variable is re-leveled to each one of the factor's levels, and the \code{mer} model updated. MCMC \emph{p}-value calculation is not yet implemented.  R will wait until all slave processes finished. At that time, the results from each slave process are collected and returned. See package \code{multicore} for more information about multicore computing.}
\usage{mcposthoc.fnc(model, data, var, two.tailed = TRUE, 
num.comp = NULL, ndigits = 4, verbose = FALSE)}
\arguments{
  \item{model}{A \code{mer} object (fitted by function \code{lmer}) or an \code{lm} object (fitted by function \code{lm}).}
  \item{data}{The unquoted name of the data frame on which the model was fitted.}
  \item{var}{A named list of variable on which to perform the posthoc analysis. For example \code{list(ph1 = c("PronomOfTheme", "AnimacyOfRec", "DefinOfRec"), ph2 = c("SemanticClass"))}.}
  \item{two.tailed}{Logical. Whether to perform one- or two-tailed \emph{t}-tests. Defaults to \code{TRUE}, i.e., two-tailed.}
  \item{num.comp}{An integer indicating the number of comparisons made. This value is used to perform Bonferroni correction of the probability values. Defaults to \code{NULL}.}
  \item{ndigits}{Integer indicating the number of decimal places to be used in the \emph{t} tables. Defaults to 4.}
  \item{verbose}{Logical. If \code{TRUE} verbose output is generated.}
}
\details{For example, if \code{var = list(ph1 = c("PronomOfTheme", "AnimacyOfRec", "DefinOfRec"))} the function will re-level and update the model on each combination of the variable levels as follows:

(1) \code{data$PronomOfTheme <- relevel(data$PronomOfTheme = "nonpronominal")}
    \code{data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "animate")}
    \code{data$DefinOfTheme <- relevel(data$DefinOfTheme = "definite")}

(2) \code{data$PronomOfTheme <- relevel(data$PronomOfTheme = "nonpronominal")}
    \code{data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "inanimate")}
    \code{data$DefinOfTheme <- relevel(data$DefinOfTheme = "definite")}

(3) \code{data$PronomOfTheme <- relevel(data$PronomOfTheme = "nonpronominal")}
    \code{data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "animate")}
    \code{data$DefinOfTheme <- relevel(data$DefinOfTheme = "indefinite")}

(4) \code{data$PronomOfTheme <- relevel(data$PronomOfTheme = "pronominal")}
    \code{data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "animate")}
    \code{data$DefinOfTheme <- relevel(data$DefinOfTheme = "definite")}

(5) \code{data$PronomOfTheme <- relevel(data$PronomOfTheme = "nonpronominal")}
    \code{data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "inanimate")}
    \code{data$DefinOfTheme <- relevel(data$DefinOfTheme = "indefinite")}

(6) \code{data$PronomOfTheme <- relevel(data$PronomOfTheme = "pronominal")}
    \code{data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "animate")}
    \code{data$DefinOfTheme <- relevel(data$DefinOfTheme = "indefinite")}

(7) \code{data$PronomOfTheme <- relevel(data$PronomOfTheme = "pronominal")}
    \code{data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "inanimate")}
    \code{data$DefinOfTheme <- relevel(data$DefinOfTheme = "indefinite")}

(8) \code{data$PronomOfTheme <- relevel(data$PronomOfTheme = "pronominal")}
    \code{data$AnimacyOfTheme <- relevel(data$AnimacyOfTheme = "inanimate")}
    \code{data$DefinOfTheme <- relevel(data$DefinOfTheme = "definite")}
}
\value{
\item{n }{The number of data points in data frame \code{data}.}        
\item{var }{A named list containing the names of the variables used in the posthoc.}
\item{summaries }{A named list containing the posthoc summaries for each re-leveling. The tables of coefficients contain the information returned from function \code{summary} as well as lower- and upper-bound \emph{df}, unadjusted lower- and upper-bound \emph{p}-values, and adjusted lower- and upper-bound \emph{p}-values if (Bonferroni) correction was enabled via argument \code{num.comp}.}
}
\author{Antoine Tremblay, IWK Health Center, \email{trea26@gmail.com}.}
\seealso{
\code{\link{posthoc.fnc} }
}
\examples{
if("multicore" \%in\% .packages(all.available = TRUE)){
library(lme4)
library(languageR)
library(multicore)
data(dative)

m0 <- glmer(RealizationOfRecipient ~ SemanticClass * 
	LengthOfRecipient + AnimacyOfRec * DefinOfRec 
	* PronomOfTheme + (1 | Speaker) + (1 | Verb), 
	data = dative)

m0.posthoc <- mcposthoc.fnc(model = m0, data = dative, 
	var = list(ph1 = c("PronomOfTheme", "AnimacyOfRec", 
	"DefinOfRec"), ph2 = c("SemanticClass")), 
	verbose = TRUE)

names(m0.posthoc)
#      [,1]       
# [1,] "n"        
# [2,] "var"      
# [3,] "summaries"

as.matrix(names(m0.posthoc$var))
#      [,1] 
# [1,] "ph1"
# [2,] "ph2"

as.matrix(names(m0.posthoc$summaries))
#     [,1]                                    
# [1,] "ph1_nonpronominal_animate_definite"    
# [2,] "ph1_nonpronominal_inanimate_definite"  
# [3,] "ph1_nonpronominal_animate_indefinite"  
# [4,] "ph1_pronominal_animate_definite"       
# [5,] "ph1_nonpronominal_inanimate_indefinite"
# [6,] "ph1_pronominal_animate_indefinite"     
# [7,] "ph1_pronominal_inanimate_indefinite"   
# [8,] "ph1_pronominal_inanimate_definite"     
# [9,] "ph2_t"                                 
# [10,] "ph2_a"                                 
# [11,] "ph2_c"                                 
# [12,] "ph2_p"                                 
# [13,] "ph2_f"     
}
}
\keyword{models & regression}
