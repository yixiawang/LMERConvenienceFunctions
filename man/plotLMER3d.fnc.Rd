\name{plotLMER3d.fnc}
\alias{plotLMER3d.fnc}
\title{Dynamic 3d plot for \code{mer} object.}
\description{Plot dynamic 3d partial effects of a (generalized) linear mixed-effects model fit with \code{LMER}.}
\usage{plotLMER3d.fnc(model = NULL, pred, intr, plot.type = "contour",
xlim = range(x, na.rm = TRUE), ylim = range(y, na.rm = TRUE), 
zlim = range(z, na.rm = TRUE), xlab = NULL, 
ylab = NULL, zlab = NULL, main = NULL, shift = 0, scale = 1, cex = 1, 
fun = NA, n = 30, color = "topo", alpha = 1, alpha.rs = 0.65, alpha.u = 1, 
lit = TRUE, theta = 0, phi = 0, contourstepsize = 0.2, legend.args = NULL, 
play3d = FALSE, ref.surf = FALSE, underneath = FALSE, add.raw = FALSE, 
color.raw = "grey", alpha.raw = 0.5, rug = FALSE, rug.u = FALSE, 
plot.dat = "default", path = "default", ...)}
\arguments{
  \item{model}{A \code{mer} object or \code{NULL} (the default) to plot from an existing data-plotting object returned by this function and saved as an \code{.rda} file.}
  \item{pred}{The quoted name of a model predictor.}
  \item{intr}{The quoted name of a continuous model predicor.}
  \item{plot.type}{The type of plot to make. Can be any of \code{"contour"} (default), \code{"image.plot"} if package \code{fields} is available, \code{"persp"}, or, if package \code{rgl} is available, \code{"persp3d"}.}
  \item{xlim, ylim, zlim}{\emph{x-}, \emph{y-} and \emph{z-}limits. The plot is produced so that the rectangular volume defined by these limits is visible.}
  \item{xlab, ylab, zlab}{Titles for the axes. N.B. These must be character strings; expressions are not accepted. Numbers will be coerced to character strings.}
  \item{main}{The main title on top of the plot.}
  \item{shift}{Constant to add to the smooth (on the scale of the linear predictor) before plotting. Defaults to 0. Passed to \code{plotRaw3d.fnc}.}
  \item{scale}{Constant by which to multiply the smooth before plotting. Defaults to 1. Passed to \code{plotRaw3d.fnc}.}
  \item{cex}{The size of label and main text.}
  \item{fun}{A function to be applied for transforming the dependent variable, if \code{NA}, no transformation is applied; for models with \code{family = "binomial"}, fun is set to \code{plogis} by default; this can be disabled by setting \code{fun=function(x)return(x)}.}
  \item{n}{Integer denoting number of points for the plot, chosen at equally spaced intervals across the empirical range of the predictor variable.}
  \item{color}{The colour scheme to use for plots. One of \code{topo}, \code{heat}, \code{cm}, \code{terrain}, \code{gray} or \code{bw}.  Schemes \code{gray}, \code{grey}, and \code{bw} also modify the colors used.}
  \item{alpha, alpha.rs, alpha.raw, alpha.u}{For \code{plot.type = "persp3d"}, alpha values between 0.0 (fully transparent) to 1.0 (opaque) for the main 3d surface, the reference surface, the added raw surface, and the "underneath" surface, respectively.}
  \item{lit}{Logical, specifying if lighting calculation should take place on geometry.}
  \item{theta}{Angle defining the viewing direction. \code{theta} gives the azimuthal direction. Used only if \code{plot.type = "persp"}.}
  \item{phi}{Angle defining the viewing direction. \code{phi} gives the colatitude. Used only if \code{plot.type = "persp"}.}
  \item{contourstepsize}{The size of the steps from contour line to contour line.}
  \item{legend.args}{When \code{plot.type = "image.plot"}, arguments for a complete specification of the legend label. This is in the form of list and is just passed to the mtext function. (See example in \code{image.plot} help page). Defaults to \code{NULL}.}
  \item{play3d}{If \code{plot.type = "persp3d"} and \code{play3d} is set to \code{TRUE}, the 3d plot will spin around \code{axis} \code{c(0, 0, 1)} at \code{rpm} \code{4} for \code{duration} \code{20} seconds. The \code{axis}, \code{rpm}, and \code{duration} can be changed by supplying a three-argument list where the first argument is a three-element vector for the rotation axis, the second argument is an interger for the rotations per minute (rpm), and the the third argument is a rotation duration time.}
  \item{ref.surf}{If \code{plot.type = "persp3d"}, whether a reference surface at the mean ought to be plotted. Defaults to \code{FALSE}.}
  \item{underneath}{If \code{plot.type = "persp3d"}, whether a flat miror image of the 3d surface ought to be plotted underneath it. Defaults to \code{FALSE}.}
  \item{add.raw}{If \code{plot.type = "persp3d"}, whether to add a surface representing the raw data. Defaults to \code{FALSE}.}
  \item{color.raw}{The colour scheme to use for the raw data surface. One of \code{topo}, \code{heat}, \code{cm}, \code{terrain}, \code{gray} or \code{bw}.  Schemes \code{gray}, \code{grey}, and \code{bw} also modify the colors used.}
  \item{rug}{Whether a rug ought to be plotted on the 3d surface. Defaults to \code{FALSE}.}
  \item{rug.u}{For \code{plot.type = "persp3d"}, whether a rug ought to be plotted on the flat miror image of the 3d surface. Defaults to \code{FALSE}.}
  \item{plot.dat, path}{Whether to cache the plotting data generated by a previous call to \code{plotLMER3d.fnc}. Generating the 3d plots can be time consuming. If the \code{plot.dat} argument is non-FALSE, the plotting information generated in the first call to the function will be saved so that in a second call to the function with exactly the same argument values, the plotting information will be retrieved and plotting will be significantly quicker. If \code{plotting.data = "default"} and \code{path = "default"}, the plotting information will be saved in a a temporary directory and the name of the file containing the information will equal to \code{paste("lmer___",model@call,pred,intr,".rda",sep="")}. The name of the file and the path where it will be saved can be set by the user in the \code{plot.dat} and \code{path} arguments. For example, \code{plot.dat = "my_plotting.data", path = "Documents"}. Note that \code{"lmer___"} will be appended to the begining of whatever is specified in \code{plot.dat} and \code{".rda"} to the end. Also note that if the user wants to save the plotting information returned by this function, the name of this object has to be \code{z}.}
  \item{...}{Further arguments to be passed to \code{image}, \code{contour}, \code{image.plot}, \code{persp}, or \code{persp3d}.}
}
\details{See help page to Harald Baayen's \code{plotLMER.fnc} function as well as to Duncan Murdoch's \code{persp3d} function and the help page to function \code{image.plot} from package \code{fields}. To save screenshots of "persp3d" plots (after plotting), use function \code{rgl.snapshot} (produces \code{png} files) or function \code{rgl.postscript} (produces \code{eps} files).}
\author{Antoine Tremblay, Dalhousie University, \email{trea26@gmail.com}.} 
\value{Invisibly returns plotting information (\code{x} and \code{y} vectors, \code{z} matrix, and colors, \code{col}). If \code{plot.type = "contour"}, \code{plot.type = "image.plot"}, or \code{plot.type = "persp"}, a contour or perspective plot, respectively. If \code{plot.type = "persp3d"}, a 2d plot as created by \code{plotLMER.fnc} as well as a dynamic 3d plot as created by \code{persp3d}. If \code{ret = TRUE}, a two-element list is returned containing the \emph{matrix} and the matrix of corresponding colors is returned. If argument \code{intel} in non-null, a file containing plotting information will be saved.}
\seealso{
\code{\link{persp}};
\code{\link{contour}};
\code{\link{plotLMER.fnc}}.
}
\examples{
if(try(require(LCFdata,quietly=TRUE))){
	data(z)
	temp.dir <- tempdir()
	save(z,file=file.path(temp.dir,"lmer___z.rda"))

  plotLMER3d.fnc(pred = "LengthBc", intr = "WMCc", 
    plot.dat = "z", path = temp.dir)
  plotLMER3d.fnc(pred = "LengthBc", intr = "WMCc",
    plot.type = "persp", phi = 25, plot.dat = "z", 
    path = temp.dir)
  if(try(require(rgl,quietly=TRUE))){
    require(rgl)
    open3d()
    plotLMER3d.fnc(pred = "LengthBc", intr = "WMCc",
      plot.type = "persp3d", plot.dat = "z", path = temp.dir)
  }
}
}
\keyword{hplot}
