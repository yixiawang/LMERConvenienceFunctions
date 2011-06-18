\name{plotLMER3d.fnc}
\alias{plotLMER3d.fnc}
\title{Dynamic 3d plot for \code{mer} object.}
\description{Plot dynamic 3d partial effects of a (generalized) linear mixed-effects model fit with \code{LMER}.}
\usage{plotLMER3d.fnc(model, pred, intr, plot.type = "contour",
xlim = range(x, na.rm = TRUE), ylim = range(y, na.rm = TRUE), 
zlim = range(z, na.rm = TRUE), xlab = "", ylab = "", zlab = "", 
cex = 1, fun = NA, n = 30, color = "topo", lit = TRUE, theta = 0, 
phi = 0, contourstepsize = 0.2, play3d = FALSE, ref.surf = FALSE,
xaxt = "s", yaxt = "s")}
\arguments{
  \item{model}{A \code{mer} object.}
  \item{pred}{The quoted name of a model predictor.}
  \item{intr}{The quoted name of a continuous model predicor.}
  \item{plot.type}{The type of plot to make. Can be any of \code{"contour"} (default), \code{"persp"}, or, if package \code{rgl} is available, \code{"persp3d"}.}
  \item{xlim, ylim, zlim}{\emph{x-}, \emph{y-} and \emph{z-}limits. The plot is produced so that the rectangular volume defined by these limits is visible.}
  \item{xlab, ylab, zlab}{Titles for the axes. N.B. These must be character strings; expressions are not accepted. Numbers will be coerced to character strings. If \code{plot.type = "contour"}, the value attributed to \code{zlab} is passed to \code{main}.}
  \item{cex}{The size of label and main text.}
  \item{fun}{A function to be applied for transforming the dependent variable, if \code{NA}, no transformation is applied; for models with \code{family = "binomial"}, fun is set to \code{plogis} by default; this can be disabled by setting \code{fun=function(x)return(x)}.}
  \item{n}{Integer denoting number of points for the plot, chosen at equally spaced intervals across the empirical range of the predictor variable.}
  \item{color}{The colour scheme to use for plots. One of ``\code{topo}'', ``\code{heat}'', ``\code{cm}'', ``\code{terrain}'', ``\code{gray}'' or ``\code{bw}''.  Schemes ``\code{gray}'' and ``\code{bw}'' also modify the colors used.}
  \item{lit}{Logical, specifying if lighting calculation should take place on geometry.}
  \item{theta}{Angle defining the viewing direction. \code{theta} gives the azimuthal direction. Used only if \code{plot.type = "persp"}.}
  \item{phi}{Angle defining the viewing direction. \code{phi} gives the colatitude. Used only if \code{plot.type = "persp"}.}
  \item{contourstepsize}{The size of the steps from contour line to contour line.}
  \item{play3d}{If \code{plot.type = "persp3d"} and \code{play3d} is set to \code{TRUE}, the 3d plot will spin around \code{axis} \code{c(0, 0, 1)} at \code{rpm} \code{4} for \code{duration} \code{20} seconds. The \code{axis}, \code{rpm}, and \code{duration} can be changed by supplying a three-argument list where the first argument is a three-element vector for the rotation axis, the second argument is an interger for the rotations per minute (rpm), and the the third argument is a rotation duration time.}
  \item{ref.surf}{If \code{plot.type = "persp3d"}, whether a reference surface at 0 units ought to be plotted. Defaults to \code{FALSE}.}
  \item{xaxt, yaxt}{Used only when \code{plot.type = "contour"}.A character which specifies the \emph{x}- and \emph{y}-axes type.  Specifying \code{"n"} suppresses plotting of the axis.  The default is \code{"s"}.}
}
\details{See help page to Harald Baayen's \code{plotLMER.fnc} function as well as to Duncan Murdoch's \code{persp3d} function.}
\author{Antoine Tremblay, IWK Health Center, \email{trea26@gmail.com}.} 
\value{If \code{plot.type = "contour"} or \code{plot.type = "persp"}, a contour or perspective plot, respectively. If \code{plot.type = "persp3d"}, a 2d plot as created by \code{plotLMER.fnc} as well as a dynamic 3d plot as created by \code{persp3d}.}
\seealso{
\code{\link{plotLMER.fnc}};
\code{\link{persp}};
\code{\link{contour}};
}
\examples{
library(lme4)
library(languageR)
data(lexdec)
lexdec3 = lexdec[lexdec$Correct == "correct", ]

m1 = lmer(RT ~ Frequency * FamilySize^2 * Length * 
    NativeLanguage + (1|Subject), data = lexdec3)

plotLMER3d.fnc(m1, pred = "Frequency", intr = "Length",
	fun = exp, plot.type = "contour", contourstepsize = 7,
	xlab = "Frequency", ylab = "Length", zlab = "Reaction Time")

dev.new()

plotLMER3d.fnc(m1, pred = "Frequency", intr = "Length",
	fun = exp, plot.type = "persp", theta = 30, phi = 30,
	xlab = "Frequency", ylab = "Length", zlab = "Reaction Time")

if("rgl" \%in\% .packages(all.available = TRUE)){
library(rgl)
dev.new()
plotLMER3d.fnc(m1, pred = "Frequency", intr = "Length",
	fun = exp, plot.type = "persp3d",xlab = "Frequency", 
	ylab = "Length", zlab = "Reaction Time")
dev.off()
}
}
\keyword{hplot}
