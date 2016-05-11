\name{mitml-package}
\alias{mitml-package}
\docType{package}
\title{mitml: Tools for multiple imputation in multilevel modeling}
\description{Provides tools for multiple imputation of missing data in multilevel modeling.
This package includes a user-friendly interface to the algorithms implemented in the \code{pan} and the \code{jomo} package, as well as several functions for visualizing, managing and the analyzing multiply imputed data sets.

The main interface to \code{pan} is the function \code{\link{panImpute}}, which allows specification of a wide array of imputation models.
In addition, the function \code{\link{jomoImpute}} provides an interface to \code{jomo}, further extending the functionality of \code{pan}.
Imputations are efficiently stored in objects of class \code{mitml}.
In order to obtain the completed (i.e., imputed) data sets from those objects, \code{\link{mitmlComplete}} is used, which produces a list of imputed data sets (class \code{mitml.list}) that can be used in further analyses.

Several functions have been added to allow for convenient analysis of multiply imputed data sets.
The functions \code{\link[=with.mitml.list]{with}} and \code{\link[=with.mitml.list]{within}} may be used for manipulating the data sets and for fitting statistical models.
Parameter estimates can be extracted using \code{\link{testEstimates}}.
Single- and multi-parameter hypotheses tests can be performed using \code{\link{testConstraints}} and \code{\link{testModels}}.
In addition, the \code{\link{anova}} method provides a simple interface to model comparisons.

Data sets can be imported and exported from or to different statistical software packages. Currently, \code{\link{write.mitmlMplus}}, \code{\link{write.mitmlSAV}} and \code{\link{write.mitmlSPSS}} export data sets to M\emph{plus} and SPSS, while other functions import previously imputed data from within R.

Finally, the package supplies tools for summarizing and visualizing imputation models, which is useful for the assessment of convergence and the reporting of results.
}
\author{
Authors: Simon Grund, Alexander Robitzsch, Oliver Luedtke

Maintainer: Simon Grund <grund@ipn.uni-kiel.de>
}
\keyword{package}
