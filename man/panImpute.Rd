\name{panImpute}
\alias{panImpute}
\title{Impute multilevel missing data using \code{pan}}
\description{This function provides a user-friendly interface to the \code{pan} package.
Imputations can be generated using \code{type} or \code{formula}, which offer different options for model specification.}

\usage{

panImpute(data, type, formula, n.burn=5000, n.iter=100, m=10, group=NULL, 
  prior=NULL, seed=NULL, save.pred=FALSE)

}

\arguments{

  \item{data}{A data frame containing incomplete and auxiliary variables, the cluster indicator variable, and any other variables that should be present in the imputed datasets.}
  \item{type}{An integer vector specifying the role of each variable in the imputation model (see details).}
  \item{formula}{A formula specifying the role of each variable in the imputation model. The basic model is constructed by \code{model.matrix}, thus allowing to include derived variables in the imputation model using \code{I()} (see details and examples).}
  \item{n.burn}{The number of burn-in iterations before any imputations are drawn. Default is to 5,000.}
  \item{n.iter}{The number of iterations between imputations. Default is to 100.}
  \item{m}{The number of imputed data sets to generate.}
  \item{group}{(optional) A character string denoting the name of an additional grouping variable to be used with the \code{formula} argument. When specified, the imputation model is run separately within each of these groups.}
  \item{prior}{(optional) A list with components \code{a}, \code{Binv}, \code{c}, and \code{Dinv} for specifying prior distributions for the covariance matrix of random effects and the covariance matrix of residuals (see details). Default is to using least-informative priors.}
  \item{seed}{(optional) An integer value initializing \code{pan}'s random number generator for reproducible results. Default is to using random seeds.}
  \item{save.pred}{(optional) Logical flag indicating if variables derived using \code{formula} should be included in the imputed data sets. Default is to \code{FALSE}.}

}

\details{

This function serves as the main interface to the \code{pan} algorithm.
The imputation model can be specified using either the \code{type} or the \code{formula} argument.

The \code{type} interface is designed to provide quick-and-easy imputations using \code{pan}. The \code{type} argument must be an integer vector denoting the role of each variable in the imputation model:
\itemize{
  \item{\code{1}: target variables containing missing data}
  \item{\code{2}: predictors with fixed effect on all targets (completely observed)}
  \item{\code{3}: predictors with random effect on all targets (completely observed)}
  \item{\code{-1}: grouping variable within which the imputation is run separately}
  \item{\code{-2}: cluster indicator variable}
  \item{\code{0}: variables not featured in the model}
}
At least one target variable and the cluster indicator must be specified.
The intercept is automatically included both as a fixed and random effect.
If a variable of type \code{-1} is found, then imputations are performed separately within each level of that variable.
This is useful if the cluster variable (e.g., schools) is contained in an even larger grouping variable for which imputations models are not deemed comparable (e.g., federal states, educational systems).

The \code{formula} argument is intended as more flexible and feature-rich interface to \code{pan}. Specifying the \code{formula} argument is similar to specifying other formulae in R.
Given below is a list of operators that \code{panImpute} currently understands:
\itemize{
  \item{\code{~}: separates the target (left-hand) and predictor (right-hand) side of the model}
  \item{\code{+}: adds target or predictor variables to the model}
  \item{\code{*}: adds an interaction term of two or more predictors}
  \item{\code{|}: denotes cluster-specific random effects and specifies the cluster indicator (i.e., \code{1|ID})}
  \item{\code{I()}: defines functions to be interpreted by \code{model.matrix}}
}
Predictors are allowed to have fixed effects, random effects or both on all target variables.
The intercept is automatically included both as a fixed and random effect, but it can be constrained if necessary (see examples).
Any predictors defined by \code{I()} will be used for imputation but not included in the data set unless \code{save.pred=TRUE}.

In order to run separate imputation models for an additional grouping variable, the \code{group} argument may be used.
The variable name must be specified without quotation marks and must be present in the data set.

As a default prior, \code{panImpute} uses least informative inverse-Wishart priors for the covariance matrices of random effects and of residuals, that is, with minimum degrees of freedom (largest dispersion) and identity matrices for scale.
For better control, the \code{prior} argument may be used for specifying alternative prior distributions.
These must be supplied as a list containing the following components:
\itemize{
  \item{\code{a}: degrees of freedom for the residual covariance matrix}
  \item{\code{Binv}: scale matrix for the residual covariance matrix}
  \item{\code{c}: degrees of freedom for the covariance matrix of random effects}
  \item{\code{Dinv}: scale matrix for the covariance matrix of random effects}
}
A sensible choice for a diffuse non-default prior is to set the degrees of freedom to the lowest value possible, and the scale matrices according to a prior guess of the corresponding covariance matrices (see Schafer & Yucel, 2002).
}

\value{

Returns an object of class \code{mitml}.
A \code{mitml} class object is a list, each containing the following components:

  \item{data}{The original (incomplete) data set that has been sorted according to the cluster variable and (if given) the grouping variable. An attribute \code{"sort"} contains the original row order. An attribute \code{"group"} contains the optional grouping variable.}
  \item{replacement.mat}{A matrix containing the multiple replacements (i.e., imputations) for each missing value. The replacement matrix contains one row for each missing value and one one column for each imputed data set.}
  \item{index.mat}{A matrix containing the row and column index for each missing value. The index matrix is used to \emph{link} the missing values in the data set with their corresponding rows in the replacement matrix.}
  \item{call}{The matched function call.}
  \item{model}{A list containing the names of the cluster variable, the target variables, and the predictor variables with fixed and random effects, respectively.}
  \item{prior}{The prior parameters used in the imputation model.}
  \item{iter}{A list containing the number of burn-in iterations, the number of iterations between imputations, and the number of imputed data sets.}
  \item{par.burnin}{A multi-dimensional array containing the parameters of the imputation model from the burn-in phase.}
  \item{par.imputation}{A multi-dimensional array containing the parameters of the imputation model from the imputation phase.}

}

\note{

For objects of class \code{mitml}, methods for the generic functions \code{print}, \code{summary} and \code{plot} have been defined.
\code{mitmlComplete} is used to extract the imputed data sets.

}

\references{
Schafer, J. L., and Yucel, R. M. (2002). Computational strategies for multivariate linear mixed-effects models with missing values. \emph{Journal of Computational and Graphical Statistics, 11}, 437-457
}

\author{Simon Grund, Alexander Robitzsch, Oliver Luedtke}
\seealso{\code{\link{mitmlComplete}}, \code{\link{summary.mitml}}, \code{\link{plot.mitml}}}

\examples{
data(studentratings)

# *** ................................
# the 'type' interface
# 

# * Example 1.1: imputation of 'LangDis' and 'SES', predicted by 'LangAchiev'

type <- c(-2,0,0,0,0,0,2,1,0,0)
names(type) <- colnames(studentratings)
type

fml <- LangDis + SES ~ LangAchiev + (1|ID)
imp <- panImpute(studentratings, type=type, n.burn=1000, n.iter=100, m=5)

# * Example 1.2: 'LangDis' and 'SES', predicted by 'LangAchiev' and 
# 'CognAbility', with random slope for 'LangAchiev'

type <- c(-2,0,0,0,0,0,3,1,2,0)
names(type) <- colnames(studentratings)
type

fml <- LangDis + SES ~ LangAchiev + (1|ID)
imp <- panImpute(studentratings, type=type, n.burn=1000, n.iter=100, m=5)

# * Example 1.3: 'LangDis' and 'SES' groupwise for 'FedState',
# and predicted by 'LangAchiev'

type <- c(-2,-1,0,0,0,0,2,1,0,0)
names(type) <- colnames(studentratings)
type

fml <- LangDis + SES ~ LangAchiev + (1|ID)
imp <- panImpute(studentratings, type=type, n.burn=1000, n.iter=100, m=5)


# *** ................................
# the 'formula' interface
# 

# * Example 2.1: imputation of 'LangDis', predicted by 'LangAchiev'
# (random intercept)

fml <- LangDis ~ LangAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

# ... the intercept can be suppressed using '0' or '-1' (here for fixed intercept)
fml <- LangDis ~ 0 + LangAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

# * Example 2.2: imputation of 'LangDis', predicted by 'LangAchiev'
# (random slope)

fml <- LangDis ~ LangAchiev + (1+LangAchiev|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

# * Example 2.3: imputation of 'LangDis', predicted by 'LangAchiev',
# groupwise for 'FedState'

fml <- LangDis ~ LangAchiev + (1|ID)
imp <- panImpute(studentratings, formula=fml, group="FedState", n.burn=1000,
n.iter=100, m=5)

# * Example 2.4: imputation of 'LangDis', predicted by 'LangAchiev'
# including the cluster mean of 'LangAchiev' as an additional predictor

fml <- LangDis ~ LangAchiev + I(clusterMeans(LangAchiev,ID)) + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5)

# ... using 'save.pred' to save the calculated cluster means in the data set
fml <- LangDis ~ LangAchiev + I(clusterMeans(LangAchiev,ID)) + (1|ID)
imp <- panImpute(studentratings, formula=fml, n.burn=1000, n.iter=100, m=5,
save.pred=TRUE)

head(mitmlComplete(imp,1))
}
\keyword{models}
