\name{as.mitml.list}
\alias{as.mitml.list}
\title{Convert a list of data sets to \code{mitml.list}}
\description{This function adds a \code{mitml.list} class attribute to a list of data frames. The resulting object can be used with the other functions of this package.}
\usage{

as.mitml.list(x)

}
\arguments{

  \item{x}{A list of data frames.}

}
\value{

The original list with an additional class attribute \code{mitml.list}.

}
\author{Simon Grund}
\seealso{\code{\link{is.mitml.list}}}

\examples{
# data frame with 'imputation' indicator
dat <- data.frame(imputation=rep(1:10,each=20), x=rnorm(200))

# split into a list and convert to 'mitml.list'
l <- split(dat, dat$imputation)
l <- as.mitml.list(l)

is.mitml.list(l)
# TRUE
}
