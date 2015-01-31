\name{plotLMER.fnc}
\alias{plotLMER.fnc}
\title{plot a mer object}
\description{
  Plot partial effects of a (generalized) linear mixed-effects model fit with
  \code{lmer} (compatible with package \code{lme4} version > 1.0).
}
\usage{
plotLMER.fnc(model, xlabel = NA, xlabs = NA, ylabel = NA, 
ylimit = NA, ilabel = NA, fun = NA, pred = NA, control = NA, 
ranefs = NA, n = 100, intr = NA,lockYlim = TRUE, addlines = FALSE, 
withList = FALSE, cexsize = 0.5, linecolor = 1, 
addToExistingPlot = FALSE, verbose = TRUE, ...)
}
\arguments{
  \item{model}{ a \code{mer} model object }
  \item{xlabel}{ label for X-axis (if other than the variable name in the
    original model formula)}
  \item{xlabs}{ character vector with labels for X-axes in multipanel plot (if 
  other than the variable names in the original model formula); if used, 
  \code{xlabel} should not be specified}
  \item{ylabel}{ label for Y-axis (if other than the variable name of 
    the dependent variable in the original model formula) }
  \item{ylimit}{ range for vertical axis; if not specified, this range will 
    be chosen such that all data points across all subplots, including HPD intervals, will be accommodated }
  \item{ilabel}{ label for the interaction shown in the lower right-hand margin of the plot, overriding the original variable name in the model formula}
  \item{fun}{ a function to be applied for transforming the dependent variable, 
    if \code{NA}, no transformation is applied; for models with \code{family = "binomial"}, 
      fun is set to \code{plogis} by default; this can be disabled by setting
      \code{fun=function(x)return(x)}.}
  \item{pred}{ character string with name of predictor; 
    if specified, a single plot will produced for the partial effect of this
    specific predictor}
  \item{control}{ a two-element list \code{list(predictor, val)} specifying a \code{predictor} the value of which has to be set to \code{val} in the partial effect plot(s); the predictor name should be exactly as specified in \code{names(model@fixef)}.  It is up to the user to make sure that name and value make sense, the code here hands full 'control' to the user.} 
  \item{ranefs}{a four-element list \code{Group, Level, Predictor}, specifying a random-effect Group (e.g. \code{Subject}), a level (e.g., \code{S10}) and a value (e.g., \code{LogFrequency}) for which partial effects have to be calibrated.}
  \item{n}{ integer denoting number of points for the plot, chosen at equally
    spaced intervals across the empirical range of the predictor variable }
  \item{intr}{ a list specifying an interaction to be graphed; obligatory
  arguments are (1) the name of the interaction variable, followed by (2) 
  a vector of values for that variable, followed by (3) the position for 
  interaction labels ('"beg"', '"mid"', or '"end"', or 'NA' if no labels are 
  desired), optionally followed by (4) a list with as first element
  a vector of colors and as second element a vector of line types. The number
  of elements in both vectors should match the number of values specified 
  under (2) for the interaction predictor.}
  \item{lockYlim}{ logical specifying whether all subplots should have the same
  range of values for the vertical axis; if \code{TRUE}, this range will be
  chosen to accomodate all fitted values including HDP intervals for all
  predictors across all plots}
  \item{addlines}{ if TRUE, adds line(s) between levels of same factor(s)}
  \item{withList}{ logical, if \code{TRUE}, a list will be output with all data 
  frames for the subplots}
  \item{cexsize}{ character expansion size (cex) for additional information in
  the plot for interactions}
  \item{linecolor}{ color of lines in the plot, by default set to 1 (black) }
  \item{addToExistingPlot}{ default FALSE, if set to TRUE, plot will be added to previous plot, but only if pred is specified}
  \item{verbose}{ if TRUE (default), effect sizes and default transformations are reported}
  \item{...}{ further graphical parameters to be passed down; warning: \code{col}, 
  \code{pch}, \code{lty} and \code{cex} will often generate an error as they are 
  internally already fully specified for specialized subplots}
}
\details{

  When no predictor is specified, a series of plots is produced for the partial
  effects of each predictor.  The graphs are shown for the reference level for
  factors and are adjusted for the median value for the other numerical
  predicors in the model.  Interactions are not shown.  The user should set up the
  appropriate number of subplots on the graphics device before running plotLMER.fnc().

  Instead of showing all predictors jointly, plotLMER.fnc() can also be used to
  plot the partial effect of a specific predictor.  When a specific predictor
  is specified (with \code{pred = ...}), a single plot is produced for that
  predictor.  In this case, the \code{intr} argument can be used to specify a
  single second predictor that enters into an interaction with the selected
  main predictor.  
  
  Polynomials have to be fitted with \code{poly(..., degree, raw=TRUE)} and
  restricted cubic splines with \code{rcs()} from the \code{rms} package.

  Note that any MCMC capabilities available in the \code{languageR} version of this function are not available in this version. 
}
\value{
  A plot is produced on the graphical device.
}
\author{R. H. Baayen, tweaked by Antoine Tremblay}
\note{ This code needs much more work, including (i) extension to \code{poly} with \code{raw=FALSE}, and (ii) general clean-up of the code.}
\seealso{
\code{\link{plotLMER3d.fnc}}.
}
\examples{
# see example in LMERConvenienceFunctions help page.
}
\keyword{ regression }
