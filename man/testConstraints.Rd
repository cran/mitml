\name{testConstraints}
\alias{testConstraints}
\title{Test functions and constraints of model parameters}
\description{
Performs hypothesis tests for arbitrary functions of a parameter vector using the Delta method.
}

\usage{

testConstraints(model, constraints, method=c("D1","D2"), df.com=NULL)

}
\arguments{

  \item{model}{A list of fitted statistical models (see examples).}
  \item{constraints}{A character vector specifying constraints or functions of the vector of model parameters to be tested.}
  \item{method}{A character string denoting the method by which the test is performed. Can be either \code{"D1"} or \code{"D2"} (see details). Default is to \code{"D1"}.}
  \item{df.com}{(optional) A single number or a numeric vector denoting the complete-data degrees of freedom for the hypothesis test. Only used if \code{method="D1"}.}

}
\details{

This function is similar in functionality to \code{testModels} but extended to arbitrary functions (or constraints) of the model parameters.
The function is based on the Delta method (e.g., Casella & Berger, 2002) according to which any function of the model parameters can be tested using Wald-like methods.
It is assumed that the parameters can be extracted using \code{coef} and \code{vcov} methods (or similar; e.g., regression coefficients, fixed effects in multilevel models)

Constraints and functions of the model parameters can be specified in the \code{constraints} argument.
The constraints must be supplied as a character vector, where each string denotes a function or a constraint to be tested (see examples).

As in \code{testModels}, the parameter vector is assumed to follow a multivariate normal distribution.
The Wald-like tests can be aggregated across data sets either by method \eqn{D_1} (Li, Raghunathan & Rubin, 1991) or \eqn{D_2} (Li, Meng, Raghunathan & Rubin, 1991), where \eqn{D_1} operates on the constrained estimates and standard errors, and \eqn{D_2} operates on the Wald-statistics.

For \eqn{D_1}, the complete-data degrees of freedom can be adjusted for smaller samples by specifying \code{df.com}.

Currently, the procedure supports statistical models that define \code{coef} and \code{vcov} methods (e.g., \code{lm}) and multilevel models estimated with \code{lme4} or \code{nlme}.
Support for further models may be added in future releases.

}
\value{

Returns a list containing the results of the model comparison, the constrained estimates and standard errors, and the relative increase in variance due to nonresponse (Rubin, 1987).
A \code{print} method is used for better readable console output. 

}

\references{
Casella, G., & Berger, R. L. (2002). \emph{Statistical inference (2nd. Ed.)}. Pacific Grove, CA: Duxbury.

Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991). Significance levels from repeated p-values with multiply-imputed data. \emph{Statistica Sinica, 1}, 65-92.

Li, K. H., Raghunathan, T. E., & Rubin, D. B. (1991). Large-sample significance levels from multiply imputed data using moment-based statistics and an F reference distribution. \emph{Journal of the American Statistical Association, 86}, 1065-1073.
}

\author{Simon Grund}

\seealso{\code{\link{testModels}}, \code{\link{with.mitml.list}}}
\examples{
data(studentratings)

fml <- MathDis + LangDis + SchClimate ~ (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

implist <- mitmlComplete(imp, print=1:5)

# fit simple regression model
fit.lm <- with(implist, lm(SchClimate ~ LangDis + MathDis))

# apply Rubin's rules
testEstimates(fit.lm)

# * Example 1: test 'identity' function of two parameters
# multi-parameter hypothesis test, thus equivalent to model comparison

constraints <- c("LangDis","MathDis")
testConstraints(fit.lm, constraints)

# ... adjusting for finite samples
testConstraints(fit.lm, constraints, df.com=749)

# ... using D2
testConstraints(fit.lm, constraints, method="D2")

# * Example 2: test for equality of two effects
# tests the hypothesis that the effects of 'LangDis' and 'MathDis'
# are equal in magnitude (LangDis=MathDis)

constraints <- c("LangDis-MathDis")
testConstraints(fit.lm, constraints)

# * Example 3: test against a fixed value
# tests if the effect of "LangDis" is different from one (LangDis=1)

constraints <- c("LangDis-1")
testConstraints(fit.lm, constraints)

# * Example 4: test function of effects
# tests if the cube of 'LangDis' is different from zero (LangDis^3=0)

constraints <- c("LangDis^3")
testConstraints(fit.lm, constraints)
}
