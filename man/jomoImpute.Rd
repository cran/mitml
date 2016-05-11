\name{jomoImpute}
\alias{jomoImpute}
\title{Impute multilevel missing data using \code{jomo}}
\description{This function provides a interface to the \code{jomo} package using the MCMC algorithms presented in Carpenter & Kenward (2013).
Using this wrapper function, \code{jomo} supports imputation of (mixed) categorical and continuous variables (Goldstein et al., 2009) as well as imputation using random (residual) covariance matrices (Yucel, 2011).
Imputations can be generated using \code{type} or \code{formula}, which offer different options for model specification.}

\usage{

jomoImpute(data, type, formula, random.L1=c("none","mean","full"), n.burn=5000,
  n.iter=100, m=10, group=NULL, prior=NULL, seed=NULL, save.pred=FALSE,
  silent=FALSE)

}

\arguments{

  \item{data}{A data frame containing incomplete and auxiliary variables, the cluster indicator variable, and any other variables that should be present in the imputed datasets.}
  \item{type}{An integer vector specifying the role of each variable in the imputation model (see details).}
  \item{formula}{A formula specifying the role of each variable in the imputation model. The basic model is constructed by \code{model.matrix}, thus allowing to include derived variables in the imputation model using \code{I()} (see details and examples).}
  \item{random.L1}{A character string denoting if the covariance matrix of residuals should vary across groups and how the values of these matrices are stored (see details). Default is to \code{"none"}, assuming a common covariance matrix across clusters.}
  \item{n.burn}{The number of burn-in iterations before any imputations are drawn. Default is to 5,000.}
  \item{n.iter}{The number of iterations between imputations. Default is to 100.}
  \item{m}{The number of imputed data sets to generate. Default is to 10.}
  \item{group}{(optional) A character string denoting the name of an additional grouping variable to be used with the \code{formula} argument. When specified, the imputation model is run separately within each of these groups.}
  \item{prior}{(optional) A list with components \code{Binv}, \code{Dinv}, and \code{a} for specifying prior distributions for the covariance matrix of random effects and the covariance matrix of residuals (see details). Default is to using least-informative priors.}
  \item{seed}{(optional) An integer value initializing R's random number generator for reproducible results. Default is to using the global seed.}
  \item{save.pred}{(optional) Logical flag indicating if variables derived using \code{formula} should be included in the imputed data sets. Default is to \code{FALSE}.}
  \item{silent}{(optional) Logical flag indicating if console output should be suppressed. Default is to \code{FALSE}.}

}

\details{

This function serves as an interface to the \code{jomo} package.
The function supports imputation of multilevel continuous and categorical data.
In order for categorical variables to be detected correctly, these must be formatted as a \code{factor} variables (see examples).
The imputation model can be specified using either the \code{type} or the \code{formula} argument.

The \code{type} interface is designed to provide quick-and-easy imputations using \code{jomo}. The \code{type} argument must be an integer vector denoting the role of each variable in the imputation model:
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

The \code{formula} argument is intended as more flexible and feature-rich interface to \code{jomo}. Specifying the \code{formula} argument is similar to specifying other formulae in R.
Given below is a list of operators that \code{jomoImpute} currently understands:
\itemize{
  \item{\code{~}: separates the target (left-hand) and predictor (right-hand) side of the model}
  \item{\code{+}: adds target or predictor variables to the model}
  \item{\code{*}: adds an interaction term of two or more predictors}
  \item{\code{|}: denotes cluster-specific random effects and specifies the cluster indicator (i.e., \code{1|ID})}
  \item{\code{I()}: defines functions to be interpreted by \code{model.matrix}}
}
Predictors are allowed to have fixed effects, random effects or both on all target variables.
The intercept is automatically included both as a fixed and random effect, but it can be constrained if necessary (see examples).
Note that, when specifying random effects other than the intercept, these will \emph{not} be automatically added as fixed effects and must be included explicitly.
Any predictors defined by \code{I()} will be used for imputation but not included in the data set unless \code{save.pred=TRUE}.

Using \code{jomo} it is possible to model the covariance matrix of residials as random across clusters using the \code{random.L1} argument (see also Yucel, 2011; Carpenter & Kenward, 2013).
In addition, the \code{random.L1} argument determines how the values of these matrices are stored.
If set to "none", a common covariance matrix is assumed across groups (the default, similar to \code{panImpute}).
If set to "mean", the covariance matrices are modeled as random, but only the average covariance matrix across clusters is stored.
If set to "full", the covariance matrices are modeled as random, and all variances and covariances from all clusters are stored.

In order to run separate imputation models for an additional grouping variable, the \code{group} argument may be used.
The variable name must be specified without quotation marks and must be present in the data set.

As a default prior, \code{jomoImpute} uses least informative inverse-Wishart priors for the covariance matrices of random effects and of residuals, that is, with minimum degrees of freedom (largest dispersion) and identity matrices for scale.
For better control, the \code{prior} argument may be used for specifying alternative prior distributions.
These must be supplied as a list containing the following components:
\itemize{
  \item{\code{Binv}: scale matrix for the residual covariance matrix}
  \item{\code{Dinv}: scale matrix for the covariance matrix of random effects}
  \item{\code{a}: starting value for the degrees of freedom of random covariance matrices of residuals (only used with \code{random.L1="mean"} or \code{random.L1="full"})}
}
Note that \code{jomo} does not allow for the degrees of freedom for the inverse-Wishart priors to be specified by the user.
These are always set to the lowest value possible (largest dispersion) or determined iteratively for the random covariance matrices of residuals (see above).

}

\value{

Returns an object of class \code{mitml}.
A \code{mitml} class object is a list, each containing the following components:

  \item{data}{The original (incomplete) data set that has been sorted according to the cluster variable and (if given) the grouping variable. An attribute \code{"sort"} contains the original row order. An attribute \code{"group"} contains the optional grouping variable.}
  \item{replacement.mat}{A matrix containing the multiple replacements (i.e., imputations) for each missing value. The replacement matrix contains one row for each missing value and one one column for each imputed data set.}
  \item{index.mat}{A matrix containing the row and column index for each missing value. The index matrix is used to \emph{link} the missing values in the data set with their corresponding rows in the replacement matrix.}
  \item{call}{The matched function call.}
  \item{model}{A list containing the names of the cluster variable, the target variables, and the predictor variables with fixed and random effects, respectively.}
  \item{random.L1}{A character string denoting the handling of random residual covariance matrices (see details).}
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
Carpenter, J. R., & Kenward, M. G. (2013). \emph{Multiple imputation and its application}. Hoboken, NJ: Wiley.

Goldstein, H., Carpenter, J., Kenward, M. G., & Levin, K. A. (2009). Multilevel models with multivariate mixed response types. \emph{Statistical Modelling}, 9, 173-197.

Yucel, R. M. (2011). Random covariances and mixed-effects models for imputing multivariate multilevel continuous data. \emph{Statistical Modelling}, 11, 351-370.
}

\author{Simon Grund, Alexander Robitzsch, Oliver Luedtke}
\seealso{\code{\link{panImpute}}, \code{\link{mitmlComplete}}, \code{\link{summary.mitml}}, \code{\link{plot.mitml}}}

\examples{
# NOTE: The number of iterations in these examples is much lower than it
# should be! This is done in order to comply with CRAN policies, and more
# iterations are recommended for applications in practice!

data(studentratings)

# *** ................................
# the 'type' interface
# 

# * Example 1.1: 'ReadDis' and 'SES', predicted by 'ReadAchiev' and 
# 'CognAbility', with random slope for 'ReadAchiev'

type <- c(-2,0,0,0,0,0,3,1,2,0)
names(type) <- colnames(studentratings)
type

imp <- jomoImpute(studentratings, type=type, n.burn=100, n.iter=10, m=5)

# * Example 1.2: 'ReadDis' and 'SES' groupwise for 'FedState',
# and predicted by 'ReadAchiev'

type <- c(-2,-1,0,0,0,0,2,1,0,0)
names(type) <- colnames(studentratings)
type

imp <- jomoImpute(studentratings, type=type, n.burn=100, n.iter=10, m=5)

# *** ................................
# the 'formula' interface
# 

# * Example 2.1: imputation of 'ReadDis', predicted by 'ReadAchiev'
# (random intercept)

fml <- ReadDis ~ ReadAchiev + (1|ID)
imp <- jomoImpute(studentratings, formula=fml, n.burn=100, n.iter=10, m=5)

# ... the intercept can be suppressed using '0' or '-1' (here for fixed intercept)
fml <- ReadDis ~ 0 + ReadAchiev + (1|ID)
imp <- jomoImpute(studentratings, formula=fml, n.burn=100, n.iter=10, m=5)

# * Example 2.2: imputation of 'ReadDis', predicted by 'ReadAchiev'
# (random slope)

fml <- ReadDis ~ ReadAchiev + (1+ReadAchiev|ID)
imp <- jomoImpute(studentratings, formula=fml, n.burn=100, n.iter=10, m=5)

# * Example 2.3: imputation of 'ReadDis', predicted by 'ReadAchiev',
# groupwise for 'FedState'

fml <- ReadDis ~ ReadAchiev + (1|ID)
imp <- jomoImpute(studentratings, formula=fml, group="FedState", n.burn=100,
n.iter=10, m=5)

# * Example 2.4: imputation of 'ReadDis', predicted by 'ReadAchiev'
# including the cluster mean of 'ReadAchiev' as an additional predictor

fml <- ReadDis ~ ReadAchiev + I(clusterMeans(ReadAchiev,ID)) + (1|ID)
imp <- jomoImpute(studentratings, formula=fml, n.burn=100, n.iter=10, m=5)

# ... using 'save.pred' to save the calculated cluster means in the data set
fml <- ReadDis ~ ReadAchiev + I(clusterMeans(ReadAchiev,ID)) + (1|ID)
imp <- jomoImpute(studentratings, formula=fml, n.burn=100, n.iter=10, m=5,
save.pred=TRUE)

head(mitmlComplete(imp,1))

# * Example 2.5: imputation of 'ReadAchiev' and 'MathAchiev' using random
# covariances matrices at level 1 (residuals)

fml <- ReadAchiev + MathAchiev ~ (1|ID)
imp <- jomoImpute(studentratings, formula=fml, random.L1="full", n.burn=100,
n.iter=10, m=5)

# * Example 2.6: imputation of 'Sex' (categorical) and 'MathAchiev' (continuous),
# predicted by 'ReadAchiev' (random slopes)

# induce some artificial missing data for 'Sex'
studentratings <- within(studentratings,{ Sex[!duplicated(ID)] <- NA
Sex <- as.factor(Sex) })

fml <- Sex + MathAchiev ~ ReadAchiev + (1+ReadAchiev|ID)
imp <- jomoImpute(studentratings, formula=fml, n.burn=100, n.iter=10, m=5)

}
\keyword{models}
