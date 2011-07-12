\name{plotLMER3d.fnc}
\alias{plotLMER3d.fnc}
\title{Dynamic 3d plot for \code{mer} object.}
\description{Plot dynamic 3d partial effects of a (generalized) linear mixed-effects model fit with \code{LMER}.}
\usage{plotLMER3d.fnc(model, pred, intr, plot.type = "contour",
xlim = range(x, na.rm = TRUE), ylim = range(y, na.rm = TRUE), 
zlim = range(z, na.rm = TRUE), xlab = NULL, ylab = NULL, zlab = NULL, 
main = NULL, cex = 1, fun = NA, n = 30, color = "topo", alpha = 1, 
alpha.rs = 0.65, alpha.u = 1, lit = TRUE, theta = 0, phi = 0, 
contourstepsize = 0.2, play3d = FALSE, ref.surf = FALSE,
underneath = FALSE, rug = FALSE, rug.u = FALSE, ret = FALSE)}
\arguments{
  \item{model}{A \code{mer} object.}
  \item{pred}{The quoted name of a model predictor.}
  \item{intr}{The quoted name of a continuous model predicor.}
  \item{plot.type}{The type of plot to make. Can be any of \code{"contour"} (default), \code{"persp"}, or, if package \code{rgl} is available, \code{"persp3d"}.}
  \item{xlim, ylim, zlim}{\emph{x-}, \emph{y-} and \emph{z-}limits. The plot is produced so that the rectangular volume defined by these limits is visible.}
  \item{xlab, ylab, zlab}{Titles for the axes. N.B. These must be character strings; expressions are not accepted. Numbers will be coerced to character strings.}
  \item{main}{The main title on top of the plot.}
  \item{cex}{The size of label and main text.}
  \item{fun}{A function to be applied for transforming the dependent variable, if \code{NA}, no transformation is applied; for models with \code{family = "binomial"}, fun is set to \code{plogis} by default; this can be disabled by setting \code{fun=function(x)return(x)}.}
  \item{n}{Integer denoting number of points for the plot, chosen at equally spaced intervals across the empirical range of the predictor variable.}
  \item{color}{The colour scheme to use for plots. One of ``\code{topo}'', ``\code{heat}'', ``\code{cm}'', ``\code{terrain}'', ``\code{gray}'' or ``\code{bw}''.  Schemes ``\code{gray}'' and ``\code{bw}'' also modify the colors used.}
  \item{alpha, alpha.rs, alpha.u}{For \code{plot.type = "persp3d"}, alpha values between 0.0 (fully transparent) to 1.0 (opaque) for the main 3d surface, the reference surface, and the "underneath" surface, respectively.}
  \item{lit}{Logical, specifying if lighting calculation should take place on geometry.}
  \item{theta}{Angle defining the viewing direction. \code{theta} gives the azimuthal direction. Used only if \code{plot.type = "persp"}.}
  \item{phi}{Angle defining the viewing direction. \code{phi} gives the colatitude. Used only if \code{plot.type = "persp"}.}
  \item{contourstepsize}{The size of the steps from contour line to contour line.}
  \item{play3d}{If \code{plot.type = "persp3d"} and \code{play3d} is set to \code{TRUE}, the 3d plot will spin around \code{axis} \code{c(0, 0, 1)} at \code{rpm} \code{4} for \code{duration} \code{20} seconds. The \code{axis}, \code{rpm}, and \code{duration} can be changed by supplying a three-argument list where the first argument is a three-element vector for the rotation axis, the second argument is an interger for the rotations per minute (rpm), and the the third argument is a rotation duration time.}
  \item{ref.surf}{If \code{plot.type = "persp3d"}, whether a reference surface at the mean ought to be plotted. Defaults to \code{FALSE}.}
  \item{underneath}{If \code{plot.type = "persp3d"}, whether a flat miror image of the 3d surface ought to be plotted underneath it. Defaults to \code{FALSE}.}
  \item{rug}{Whether a rug ought to be plotted on the 3d surface. Defaults to \code{FALSE}.}
  \item{rug.u}{For \code{plot.type = "persp3d"}, whether a rug ought to be plotted on the flat miror image of the 3d surface. Defaults to \code{FALSE}.}
  \item{ret}{Whether the \emph{z} matrix and corresponding color matrix ought to be returned. Defaults to \code{FALSE}.}
}
\details{See help page to Harald Baayen's \code{plotLMER.fnc} function as well as to Duncan Murdoch's \code{persp3d} function. To save screenshots of "persp3d" plots (after plotting), use function \code{rgl.snapshot} (produces \code{png} files) or function \code{rgl.postscript} (produces \code{eps} files).}
\author{Antoine Tremblay, IWK Health Center, \email{trea26@gmail.com}.} 
\value{If \code{plot.type = "contour"} or \code{plot.type = "persp"}, a contour or perspective plot, respectively. If \code{plot.type = "persp3d"}, a 2d plot as created by \code{plotLMER.fnc} as well as a dynamic 3d plot as created by \code{persp3d}. If \code{ret = TRUE}, a two-element list is returned containing the \emph{matrix} and the matrix of corresponding colors is returned.}
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

devAskNewPage(ask=FALSE)
plotLMER3d.fnc(m1, pred = "Frequency", intr = "Length",
	fun = exp, plot.type = "contour", contourstepsize = 7,
	xlab = "Frequency", ylab = "Length", zlab = "Reaction Time")

dev.new()
devAskNewPage(ask=FALSE)
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
