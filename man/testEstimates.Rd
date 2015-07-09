\name{testEstimates}
\alias{testEstimates}
\title{Compute final estimates and inferences}
\description{
Computes final parameter estimates and inferences from multiply imputed data sets.
}

\usage{

testEstimates(model, var.comp=FALSE, df.com=NULL)

}
\arguments{

  \item{model}{A list of fitted statistical models (see examples).}
  \item{var.comp}{A logical flag indicating if estimates for variance components should be calculated. Default is to \code{FALSE}.}
  \item{df.com}{(optional) A numeric vector denoting the complete-data degrees of freedom for the hypothesis test.}

}
\details{

This function calculates final parameter estimates and inferences as suggested by Rubin (1987, "Rubin's rules") for each parameter of the fitted model.
In other words, \code{testEstimates} aggregates estimates and standard errors across multiply imputed data sets.

Rubin's original method assumes that the complete-data degrees of freedom are infinite, which is reasonable in larger samples.
Alternatively, the degrees of freedom can be adjusted for smaller samples by specifying \code{df.com} (Barnard & Rubin, 1999).
The \code{df.com} argument can either be a single number if the degrees of freedom are equal for all tests, or a numeric vector with one element per test.

Using the \code{var.comp} argument, final estimates for variance components and related parameters can be requested.
These will be shown as a separate table within the console output.
Accessing variance components is highly dependent on the model being estimated and not implemented for all models. 
Users may prefer calculating these estimates manually using \code{with.mitml.list} (see Example 3).
No inferences are calculated for variance components.

Currently, the procedure supports statistical models that define \code{coef} and \code{vcov} methods (e.g., \code{lm}) and multilevel models estimated with \code{lme4} or \code{nlme}.
Support for further models may be added in future releases.

}
\value{

Returns a list containing the final parameter and inferences as well as the relative increase in variance due to nonresponse and the fraction of missing information (Rubin, 1987).
A \code{print} method is used for better readable console output. 

}

\references{
Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with multiple imputation. \emph{Biometrika, 86}, 948-955.

Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}. Hoboken, NJ: Wiley.
}


\author{Simon Grund}

\seealso{\code{\link{with.mitml.list}}}
\examples{
data(studentratings)

fml <- ReadDis + SES ~ ReadAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

implist <- mitmlComplete(imp, print=1:5)

# fit multilevel model using lme4
require(lme4)
fit.lmer <- with(implist, lmer(SES ~ (1|ID)))

# final estimates and inferences sperately for each parameter (Rubin's rules)
testEstimates(fit.lmer)

# ... adjusted df for finite samples
testEstimates(fit.lmer, df.com=49)

# ... with additional table for variance components and ICCs
testEstimates(fit.lmer, var.comp=TRUE)
}
