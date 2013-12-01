\name{mrmrPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
%\alias{summary.cmim}
\title{
MRMR PCA
}
\description{
PCA using information theoretic measures.
}

\usage{
mrmrPCA(x, y = c(), method="CI")
}


\arguments{
  \item{x} {
    A matrix of variables.
  }

  \item{y} {
    A vector describing the class for each row.
  }

  \item(method) {
    Method for measuring variance and covariance between columns.
  }

}

\value{
  \code{mrmrPCA} returns an list with class "princomp" as returned by
  the \code{princomp} function with the following components:

  \item{sdev} { e$values + abs(min(e$values)) }
  \item{loadings} { the matrix of variable loadings (i.e., a matrix
  the columns of which contain the eigenvectors). TODO: make it of class
  \code{loadings}}
  \item{n.obs} { the number of observations. }
  \item{scale} { currently all 1s because I am not sure what this is. }
  \item{scores} { the data rotated along the principal components. }
  \item{call} { the matched call. }
}

\details{
%  ~~ If necessary, more details than the description above ~~
TODO: list all supported methods.
}

\references{
% ~put references to the literature/web site here ~
TODO: list references
}

\author{
Diman Todorov
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
