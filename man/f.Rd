\name{f}
\alias{f}
\title{List files and directories in current directory.}
\description{List files and directories in current directory in matrix format. Each row is preceded by a row number.}
\usage{f(path = ".", pattern = NULL, all.files = FALSE, 
full.names = FALSE, recursive = FALSE, ignore.case = FALSE)}
\arguments{
  \item{path}{A character vector of full path names; the default corresponds to the working directory \code{getwd()}. Missing values will be ignored.}
  \item{pattern}{An optional regular expression. Only file names which match the regular expression will be returned.}
  \item{all.files}{Logical. If \code{FALSE}, only the names of visible files are returned.  If \code{TRUE}, all file names will be returned.}
  \item{full.names}{Logical. If \code{TRUE}, the directory path is prepended to the file names.  If \code{FALSE}, only the file names are returned.}
  \item{recursive}{Logical. Should the listing recurse into directories?}
  \item{ignore.case}{Logical. Should pattern-matching be case-insensitive?}
}
\value{A matrix containing the names of the files and directories, preceded by a row number, in the specified directories. If a path does not exist or is not a directory or is unreadable it is skipped, with a warning.  

The files are sorted in alphabetical order, on the full path if \code{full.names = TRUE}. Directories are included only if \code{recursive = FALSE}.}
\author{Antoine Tremblay, Statistics Canada, \email{trea26@gmail.com}}
\note{File naming conventions are platform dependent.  \code{recursive = TRUE} is not supported on all platforms and may be ignored (with a warning).}
\seealso{
\code{\link{list.files}}
}
\examples{
f()
}
