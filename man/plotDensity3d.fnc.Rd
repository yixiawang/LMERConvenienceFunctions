\name{plotDensity3d.fnc}
\alias{plotDensity3d.fnc}
\title{Kernel density estimation for two continuous variables.}
\description{The densities of two continuous variables is first computed using the \code{density} function from package \code{stats}. The outer product of the two densities is computed, which can be plotted as a contour map or a perspective plot.}
\usage{plotDensity3d.fnc(x, y, plot.type = "contour", color = "terrain", 
xlab = NULL, ylab = NULL, zlab = NULL, main = NULL, cex = 1, 
lit = TRUE, theta = 0, phi = 0, bw = "nrd0", adjust = 1, kernel = c("gaussian", 
"epanechnikov", "rectangular", "triangular", "biweight", "cosine", 
"optcosine"), weights = NULL, window = kernel, width, give.Rkern = FALSE, 
n = 50, from, to, cut = 3, na.rm = FALSE, ...)}
\arguments{
  \item{x, y}{Numeric vectors.}
  \item{plot.type}{The type of plot to make. Can be any of \code{"contour"} (default) or \code{"persp"}.}
  \item{color}{The colour scheme to use for plots. One of ``\code{topo}'', ``\code{heat}'', ``\code{cm}'', ``\code{terrain}'', ``\code{gray}'' or ``\code{bw}''.  Schemes ``\code{gray}'' and ``\code{bw}'' also modify the colors used.}
  \item{xlab, ylab, zlab}{Titles for the axes. N.B. These must be character strings; expressions are not accepted. Numbers will be coerced to character strings.}
  \item{main}{The main title on top of the plot.}
  \item{cex}{The size of label and main text.}
  \item{lit}{Logical, specifying if lighting calculation should take place on geometry.}
  \item{theta}{Angle defining the viewing direction. \code{theta} gives the azimuthal direction. Used only if \code{plot.type = "persp"}.}
  \item{phi}{Angle defining the viewing direction. \code{phi} gives the colatitude. Used only if \code{plot.type = "persp"}.}
  \item{bw, adjust, kernel, weights, window, width, give.Rkern, n, from, to, cut, na.rm}{See help page to function \code{density}.}
  \item{...}{Further arguments passed to functions \code{image}, \code{contour} or \code{persp}.}
}
\details{See help page to the \code{density} function.}
\value{Either a contour map or a perspective plot. Invisibly returns
\item{x }{The numeric vector supplied in argument \code{x}.}
\item{y }{The numeric vector supplied in argument \code{y}.}
\item{xd }{The density object tied to vector \code{x}.}
\item{yd }{The density object tied to vector \code{y}.}
\item{mat }{The outer product of the \code{x} and \code{y} densities in matrix format.}
\item{col }{The color used for plotting.}
}
\author{Antoine Tremblay, Statistics Canada, \email{trea26@gmail.com}.} 
\seealso{
\code{\link{contour}};
\code{\link{persp}};
\code{\link{density}};
\code{\link{outer}}
}
\examples{
# see example in LMERConvenienceFunctions help page.
}
\keyword{hplot}
\keyword{dynamic}
\keyword{multivariate}
