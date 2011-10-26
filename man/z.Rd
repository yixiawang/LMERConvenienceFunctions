\name{z}
\alias{z}
\docType{data}
\title{Plotting data generated from a linear mixed-effects model from Tremblay & Newman (In Preparation).}
\description{
The model call is:
	lmer(formula = Amplitude ~ FreqBc * LengthBc * WMCc + 
	(1 | Subject) + (1 | Item) + (0 + WMCc | Item), data = dat)

(see \emph{details} for more).

The plotting data was generated from this model with function \code{plotLMER3d.fnc}.
}
\usage{data(z)}
\format{
  The format is:
 num [1:30, 1:10] -1.83 -1.95 -2.07 -2.19 -2.32 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:30] "-2.88262910798122" "-2.50331876315363" "-2.12400841832605" "-1.74469807349846" ...
  ..$ : chr [1:10] "-0.1857142855" "-0.1410714285" "-0.0517857144999999" "-0.0428571424999999" ...
}
\details{The model includes a three-way interaction between \emph{WMCc} (mean-centered working memory capacity), \emph{FreqBc} (the frequency of use of the second word of a four-word sequence), and \emph{LengthBc} (the length in number of letters of the second work of a four-word sequence) in addition to by-subject and by-item random intercepts and by-item random slopes for \emph{WMCc}.}
\source{
The data is from:

Tremblay, Antoine. (2009). \emph{Processing Advantages of Lexical Bundles: Evidence from Self-paced Reading, Word and Sentence Recall, and Free Recall with Event-related Brain Potential Recordings}. Ph.D. Dissertation. University of Alberta, Edmonton, Canada. Available for download at \url{http://www.ualberta.ca/~antoinet/ThesisDraft_10_B.pdf}.


The model is from:

Tremblay, Antoine, and Newman, Aaron J. (In Preparation). \emph{The Analysis of Event-related Potentials using Linear Mixed-effects Models with Complex Random-effect Structures}.
}
\examples{
data(z)
str(z)
}
\keyword{datasets}
