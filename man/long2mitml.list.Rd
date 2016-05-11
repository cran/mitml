\name{long2mitml.list}
\alias{long2mitml.list}
\alias{jomo2mitml.list}
\title{Convert imputations from long format to \code{mitml.list}}
\description{These functions convert data sets containing multiple imputations in long format to objects of class \code{mitml.list} for use in further analyses.}
\usage{

long2mitml.list(x, split, exclude=NULL)

jomo2mitml.list(x)

}
\arguments{

  \item{x}{A data frame in long format containing multiple imputations (see details).}
  \item{split}{A character string naming the column in \code{x} that denotes different imputatins (see details).}
  \item{exclude}{A vector denoting values of \code{split} which should be excluded from the list.}

}
\details{

The general function \code{long2mitml.list} is intedended for converting data frames from the long format to \code{mitml.list} (i.e., a list of imputed data sets).
In this format, imputations are enclosed in a single data frame where some \code{split} variables is used to denote different imputations.

The function \code{jomo2mitml.list} is a special case of \code{long2mitml.list} and is intended for converting imputations generated with the native functions of the \code{jomo} package.

}
\value{

A list of imputed data sets with an additional class attribute \code{mitml.list}.

}
\author{Simon Grund}
\seealso{\code{\link{mitmlComplete}}}

\examples{
data(studentratings)
require(jomo)

# impute data using jomo (native functions)
clus <- studentratings[,"ID"]
imp <- jomo(Y=studentratings[,c("ReadAchiev","ReadDis")], clus=clus)

# split imputations
impList <- long2mitml.list(imp, split="Imputation", exclude=0)
impList <- jomo2mitml.list(imp)
}
