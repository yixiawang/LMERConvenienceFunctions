\name{summary.mcposthoc}
\alias{summary.mcposthoc}
\title{Summarize a "mcposthoc" object.}
\description{This function extracts the desired portions of an "mcposthoc" object.}
\usage{\method{summary}{mcposthoc}(object, ph.list = NULL, 
term = NULL, print = TRUE, ...)}
\arguments{
  \item{object}{An "mcposthoc" object as returned by function \code{mcposthoc.fnc}.}
  \item{ph.list}{The name of the posthoc analysis for which results are desired. For example, if, in function \code{mcposthoc.fnc}, argument \code{var} was set to \code{list(ph1 = c("PronomOfTheme", "AnimacyOfRec", "DefinOfRec"), ph2 = c("SemanticClass"))}, \code{ph.list} should be one of "ph1" or "ph2". Defaults to \code{NULL}. If \code{ph.list = NULL} and more than one posthoc analysis was performed, the user will be prompted to select one of the analyses.}
  \item{term}{The model term for which posthoc results are desired. Defaults to \code{NULL}, in which case the user will be prompted to select a term.}
  \item{print}{Whether to print to screen the posthoc summary. Defaults to \code{TRUE}.}
  \item{\dots}{Not used.}
}
\details{The function creates a summary data frame from statistics obtained from an "mcposthoc" object for the specified term. It goes through each element of the \code{ph.list} -- each list element is the summary of the model re-leveled on one factor level (or combination of factor levels) -- extracts the row corresponding to the \code{term}, and binds it to the other extracted rows.}
\value{
\item{ph.list }{The posthoc list in the "mcposthoc" object from which the summary originates.}
\item{term }{The term from the posthoc list for which a summary is desired.}
\item{summary }{The posthoc summary.}
}
\author{Antoine Tremblay, Dalhousie University, \email{trea26@gmail.com}}
\seealso{
\code{\link{mcposthoc.fnc}; }
\code{\link{pamer.fnc}. }
}
\examples{
### See examples from mcposthoc.fnc() help page.
}
\keyword{models & regression}
