\name{plotRaw3d.fnc}
\alias{plotRaw3d.fnc}
\title{Visualize raw surface averages (3d)}
\description{For a specified response variable and interacting continuous predictors, visualize in 3d the surface average.}
\usage{
plotRaw3d.fnc(data = NULL, response = NULL, pred = NULL, intr = NULL, 
xy = TRUE, color = "topo", zlim = NULL, xlab = NULL, ylab = NULL, 
zlab = NULL, main = NULL, shift = 0, scale = 1, plot.type = "contour", 
theta = 30, phi = 30, ticktype = "detailed", 
contourstepsize = 1, legend.args = NULL, ...)
}
\arguments{
  \item{data}{A data frame.}
  \item{response}{The quoted name of a continuous response variable.}
  \item{pred}{The quoted name of a continuous predictor.}
  \item{intr}{The quoted name of an interacting continuous predictor.}
  \item{xy}{Whether to the \code{x} and \code{y} values from the data or to set them to \code{seq(0, 1, len = nrow(z))}. Defaults to \code{TRUE}.}
  \item{color}{The colour scheme to use. One of \code{"topo"}, \code{"heat"}, \code{"cm"}, \code{"terrain"}, \code{"gray"} or \code{"bw"}.}
  \item{zlim}{A two element vector specifying the plotting limits for the \emph{z}-axis.}
  \item{xlab, ylab, zlab}{Titles for the axes. N.B. These must be character strings; expressions are not accepted. Numbers will be coerced to character strings.}
  \item{main}{The main title on top of the plot.}
  \item{shift}{Constant to add to the smooth (on the scale of the linear predictor) before plotting. Defaults to 0.}
  \item{scale}{Constant by which to multiply the smooth before plotting. Defaults to 1.}
  \item{plot.type}{The type of plot to make. Can be any of \code{"contour"} or \code{"persp"}.}
  \item{theta}{Angle defining the viewing direction. \code{theta} gives the azimuthal direction.}
  \item{phi}{Angle defining the viewing direction. \code{phi} gives the colatitude.}
  \item{ticktype}{Character: \code{"simple"} draws just an arrow parallel to the axis to indicate direction of increase; \code{"detailed"} draws normal ticks as per 2D plots.}
  \item{contourstepsize}{The size of the steps from contour line to contour line. Defaults to 1. Used only if \code{plot.type = "contour"}.}
  \item{legend.args}{When \code{plot.type = "image.plot"}, arguments for a complete specification of the legend label. This is in the form of list and is just passed to the mtext function. (See example in \code{image.plot} help page). Defaults to \code{NULL}.}
  \item{...}{Further arguments passed to functions \code{image}, \code{image.plot}, \code{contour} or \code{persp}.}

}
\details{\code{NA}s will be set to \code{0}.}
\value{Either a dynamic 3d perspective plot, a perspective plot, or a contour plot. Also invisibly returns the plotting matrix and the color vector.}
\author{Antoine Tremblay, Statistics Canada \email{trea26@gmail.com}}
\examples{
# see example in LMERConvenienceFunctions help page.
}
\keyword{ hplot }
