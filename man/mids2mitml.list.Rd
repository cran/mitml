\name{mids2mitml.list}
\alias{mids2mitml.list}
\title{Convert objects of class \code{mids} to \code{mitml.list}}
\description{This function converts a \code{mids} class object (from the \code{mice} package) to \code{mitml.list}, such that it may be managed and analyzed using the functions of this package.}
\usage{

mids2mitml.list(x)

}
\arguments{

  \item{x}{An object of class \code{mids} as produced by \code{mice} (see the \code{mice} package).}

}
\value{

A list of imputed data sets with an additional class attribute \code{mitml.list}.

}
\author{Simon Grund}
\seealso{\code{\link{mitmlComplete}}}

\examples{
data(studentratings)

# imputation using mice
require(mice)
imp <- mice(studentratings)

implist <- mids2mitml.list(imp)
}
