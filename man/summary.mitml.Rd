\name{summary.mitml}
\alias{summary.mitml}
\title{Summary measures for imputation models}
\description{
Provides summary statistics and additional information regarding the imputation process.
}

\usage{

\method{summary}{mitml}(object, n.Rhat=3, goodness.of.appr=FALSE, ...)

}

\arguments{

  \item{object}{An object of class \code{mitml} as produced by \code{panImpute}.}
  \item{n.Rhat}{(optional) An integer denoting the number of sequences used for calculating the potential scale reduction factor. Default is to \code{3}.}
  \item{goodness.of.appr}{(optional) A logical flag indicating if the goodness of approximation should be printed. Default is to \code{FALSE} (see Note).}
  \item{\dots}{Not being used.}

}

\details{

The \code{summary} method calculates summary statistics for objects of class \code{mitml}, gives information on the imputation process, and the amount of missing data per variable.

The output includes the potential scale reduction factor (PSRF, or \eqn{\hat{R}}), which is calculated for each parameter of the imputation model as a measure of convergence (Gelman and Rubin, 1992).
Calculation of the PSRFs can be suppressed by setting \code{n.Rhat=NULL}.

}

\value{

Returns an object of class \code{summary.mitml}. A print method is used for better readable console output.

}

\note{

The PSRFs are not computed from different chains, but by dividing each chain from the imputation phase into a number of sequences as denoted by \code{n.Rhat}.
This is slightly different from the original method proposed by Gelman and Rubin.

The goodness of approximation indicates what proportion of the posterior standard deviation is due to simulation error.
For multiple imputation, goodness of approximation is not essential; it should be considered only if posterior summaries, such as the EAP, are of interest.

}

\references{
Gelman, A., and Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. \emph{Statistical Science, 7}, 457-472.

Hoff, P. D. (2009). \emph{A first course in Bayesian statistical methods}. New York, NY: Springer.
}

\author{Simon Grund}
\seealso{\code{\link{panImpute}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

# print summary
summary(imp)
}
\keyword{methods}
