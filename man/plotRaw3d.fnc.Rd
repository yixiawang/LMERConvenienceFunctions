\name{plotRaw3d.fnc}
\alias{plotRaw3d.fnc}
\title{Visualize raw surface averages (3d)}
\description{For a specified response variable and interacting continuous predictors, visualize in 3d the surface average.}
\usage{
plotRaw3d.fnc(data = NULL, response = NULL, pred = NULL, intr = NULL, 
color = "topo", zlim = NULL, xlab = NULL, ylab = NULL, zlab = NULL, 
main = NULL, shift = NULL, plot.type = "persp", theta = 30, phi = 30, 
ticktype = "detailed", contourstepsize = 1, ret = FALSE)
}
\arguments{
  \item{data}{A data frame.}
  \item{response}{The quoted name of a continuous response variable.}
  \item{pred}{The quoted name of a continuous predictor.}
  \item{intr}{The quoted name of an interacting continuous predictor.}
  \item{color}{The colour scheme to use. One of \code{"topo"}, \code{"heat"}, \code{"cm"}, \code{"terrain"}, \code{"gray"} or \code{"bw"}.}
  \item{zlim}{A two element vector specifying the plotting limits for the \emph{z}-axis.}
  \item{xlab, ylab, zlab}{Titles for the axes. N.B. These must be character strings; expressions are not accepted. Numbers will be coerced to character strings.}
  \item{main}{The main title on top of the plot.}
  \item{shift}{Constant to add to each smooth (on the scale of the linear predictor) before plotting.}
  \item{plot.type}{The type of plot to make. Can be any of \code{"contour"}, \code{"persp"}, the default, or, if package \code{rgl} is available, \code{"persp3d"}.}
  \item{theta}{Angle defining the viewing direction. \code{theta} gives the azimuthal direction.}
  \item{phi}{Angle defining the viewing direction. \code{phi} gives the colatitude.}
  \item{ticktype}{Character: \code{"simple"} draws just an arrow parallel to the axis to indicate direction of increase; \code{"detailed"} draws normal ticks as per 2D plots.}
  \item{contourstepsize}{The size of the steps from contour line to contour line. Defaults to 1. Used only if \code{plot.type = "contour"}.}
  \item{ret}{Whether the \emph{z} matrix and corresponding color matrix ought to be returned. Defaults to \code{FALSE}.}
}
\details{\code{NA}s will be set to \code{0}.}
\value{Either a dynamic 3d perspective plot, a perspective plot, or a contour plot.}
\author{Antoine Tremblay, IWK Health Center \email{trea26@gmail.com}}
\examples{
data(erpFz)
# contour plot example
devAskNewPage(ask=FALSE)
plotRaw3d.fnc(data = erpFz, response = "Amplitude", pred = "LengthBc", intr = "WMCc",
	plot.type = "contour", zlim = c(-10,2))

# perspective plot example
dev.new()
devAskNewPage(ask=FALSE)
plotRaw3d.fnc(data = erpFz, response = "Amplitude", pred = "LengthBc", intr = "WMCc",
	plot.type= "persp", theta = 50, phi = 20, zlim = c(-10,2))

# if package rgl is available, a dynamic perspective plot example
if("rgl" \%in\% .packages(all.available = TRUE)){
library(rgl)
open3d()
plotRaw3d.fnc(data = erpFz, response = "Amplitude", pred = "LengthBc", intr = "WMCc",
	plot.type = "persp3d", zlim = c(-10,2))
}
}
\keyword{ hplot }
