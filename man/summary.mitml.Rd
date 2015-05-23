\name{summary.mitml}
\alias{summary.mitml}
\title{Summary measures for imputation models}
\description{
Provides summary statistics and additional information regarding the imputation process.
}

\usage{

\method{summary}{mitml}(object, n.Rhat=3, ...)

}

\arguments{

  \item{object}{An object of class \code{mitml} as produced by \code{panImpute}.}
  \item{n.Rhat}{An integer denoting the number of sequences used for calculating the potential scale reduction factor. The imputation phase is divided into a number of sequences for which the potential scale reduction is is calculated. Default is to \code{3}.}
  \item{\dots}{Not being used.}

}

\details{

The \code{summary} method calculates summary statistics for objects of class \code{mitml}, gives information on the imputation process, and the amount of missing data per variable.

The console output features the potential scale reduction factor (PSRF, also called \eqn{\hat{R}}), which is calculated for each parameter of the imputation model to assess convergence of the Gibbs sampler (Gelman and Rubin, 1992).
Due to the potentially large number of parameters in the imputation model, the console output features a summary of all PSRFs rather than the individual values.

Note that the PSRFs are not computed from different chains, but by dividing the imputation phase into separate sequences.
This is different from the original method proposed by Gelman and Rubin.

}

\value{

Returns an object of class \code{summary.mitml}. A print method is used for better readable console output.

}

\references{
Gelman, A., and Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. \emph{Statistical Science, 7}, 457-472.
}

\author{Simon Grund}
\seealso{\code{\link{panImpute}}}
\examples{
data(studentratings)

fml <- LangDis + SES ~ LangAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

# print summary
summary(imp)

# internal structure of the summary-object
str(summary(imp))
}
\keyword{methods}
