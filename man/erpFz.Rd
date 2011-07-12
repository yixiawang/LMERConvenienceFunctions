\name{erpFz}
\alias{erpFz}
\docType{data}
\title{ERP amplitudes at electrode Fz restricted to the 100 to 175 millisecond time window.}
\description{Event-related brain potential data recorded from 10 participants at electrode Fz (100--175 ms) in the context of Antoine Tremblay's PhD thesis (Tremblay, 2009). The actual raw data was averaged over variables \emph{LengthBc} and \emph{WMCc}.}
\usage{data(erpFz)}
\format{
  A data frame with 120 observations on the following 3 variables.
  \describe{
    \item{\code{LengthBc}}{A numeric vector. The mean-centered length (in number of letters) of the second word of a four-word sequence.}
    \item{\code{WMCc}}{A numeric vector. The mean-centered working memory capacity of the participants.}
    \item{\code{Amplitude}}{A numeric vector. The avereaged amplitude of the recorded ERPs.}
  }
}
\source{Tremblay, Antoine. (2009). \emph{Processing Advantages of Lexical Bundles: Evidence from Self-paced Reading, Word and Sentence Recall, and Free Recall with Event-related Brain Potential Recordings. Ph.D. Dissertation. University of Alberta, Edmonton, Canada. Available for download at \url{http://www.ualberta.ca/~antoinet/ThesisDraft_10_B.pdf}.}
}
\examples{
data(erpFz)
summary(erpFz)
}
\keyword{datasets}
