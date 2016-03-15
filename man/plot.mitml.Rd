\name{plot.mitml}
\alias{plot.mitml}
\title{Print diagnostic plots}
\description{
Generates diagnostic plots for assessing the convergence and autocorrelation behavior of \code{pan}'s and \code{jomo}'s MCMC algorithms.
}

\usage{

\method{plot}{mitml}(x, print=c("beta","psi","sigma"), group="all",
  trace=c("imputation","burnin","all"), smooth=3,
  export=c("none","png","pdf"), dev.args=list(), ...)

}

\arguments{

  \item{x}{An object of class \code{mitml} as produced by \code{panImpute} and \code{jomoImpute}.}
  \item{print}{A character vector containing one or several of \code{"beta"}, \code{"psi"} or \code{"sigma"} denoting which parameters should be plotted. Default is to plot all parameters.}
  \item{group}{Either \code{"all"} or an integer denoting for which group plots should be generated. Used only when \code{group} has been specified in \code{panImpute} for group-wise imputation.}
  \item{trace}{One of \code{"imputation"}, \code{"burnin"} or \code{"all"} denoting which part of the parameter chain should be used for the trace plot. Default plots only the iterations after burn-in.}
  \item{smooth}{A numeric value denoting the smoothing factor for the trend line in trace plots. Higher values correspond to less smoothing. Default is \code{3}. If set to \code{0} or \code{NULL}, the trend line is suppressed.}
  \item{export}{(optional) A character string specifying if plots should be exported to file. If \code{"png"} or \code{"pdf"}, then plots are printed into a folder named "panPlots" in the current directory using either the \code{png} or \code{pdf} device. Default is to \code{"none"}, which does not export files.}
  \item{dev.args}{(optional) A named list containing additional arguments that are passed to the graphics device.}
  \item{\dots}{Parameters passed to the plotting functions.}

}

\details{

The \code{plot} method generates a series of plots for the parameters of the imputation model which can be used for diagnostic purposes.

Setting \code{print} to \code{"beta"}, \code{"psi"} and \code{"sigma"} will plot the fixed effects, the variances and covariances of random effects, and the variances and covariances of residuals, respectively.
Each plotting window contains a trace plot (upper left), an autocorrelation plot (lower left), a kernel density approximation of the posterior distribution (upper right), and a posterior summary (lower right).
The \code{trace} and \code{smooth} arguments can be used to influence how the trace plot is drawn, and what part of the chain should be used for it.

The plots are presented one at a time.
To proceed with the next plot, the user may left-click in the plotting window or press "enter" while in the R console.
No plots are displayed when exporting to file.

}

\value{

None (invisible \code{NULL}).

}

\author{Simon Grund}
\seealso{\code{\link{panImpute}}, \code{\link{jomoImpute}}}
\examples{
\dontrun{
data(studentratings)

# * Example 1: simple imputation

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

plot all parameters (default)
plot(imp)
plot(imp, print=c("beta","psi","sigma"))

plot fixed effects only
plot(imp, print="beta")

# export plots to file (using pdf device)
plot(imp, export="pdf", dev.args=list(width=9, height=4, pointsize=12))

# * Example 2: groupwise imputation

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, group=FedState, n.burn=1000,
       n.iter=100, m=5)

plot fixed effects for all groups (default for 'group')
plot(imp, print="beta", group="all")

plot fixed effects for first group only
plot(imp, print="beta", group=1)
}
}
\keyword{methods}
