\name{f}
\alias{f}
\title{
List files and directories in current directory.
}
\description{
List files and directories in current directory in matrix format. Each row is preceded by a row number.
}
\usage{
f(path = ".", pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
}
\arguments{
  \item{path}{
A character vector of full path names; the default corresponds to the working directory `getwd()'. Missing values will be ignored.
}
  \item{pattern}{
An optional regular expression. Only file names which match the regular expression will be returned.
}
  \item{all.files}{
A logical value. If `FALSE', only the names of visible files are returned.  If `TRUE', all file names will be returned.
}
  \item{full.names}{
A logical value. If `TRUE', the directory path is prepended to the file names.  If `FALSE', only the file names are returned.
}
  \item{recursive}{
Logical. Should the listing recurse into directories?
}
  \item{ignore.case}{
Logical. Should pattern-matching be case-insensitive?
}
}
\value{
A matrix containing the names of the files and directories, preceded by a row number, in the specified directories, or `""' if there were no files. If a path does not exist or is not a directory or is unreadable it is skipped, with a warning.

The files are sorted in alphabetical order, on the full path if `full.names = TRUE'. Directories are included only if `recursive = FALSE'.
}
\author{
Antoine Tremblay, IWK Health Center, <trea26@gmail.com>
}
\note{
File naming conventions are platform dependent.

`recursive = TRUE' is not supported on all platforms and may be ignored (with a warning).
}

\seealso{
\code{\link{list.files}}
}
\examples{
f()
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
