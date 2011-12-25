\name{UPShclus}
\alias{UPShclus}
\title{Hierarchical Clustering of Patients on X-covariates for Unsupervised Propensiy Scoring}
\description{Derive a full, hierarchical clustering tree (dendrogram) for all patients (regardless
 of treatment received) using Mahalonobis between-patient distances computed from specified
 baseline X-covariate characteristics.}
\usage{
  UPShclus(dframe, xvars, method="diana")
}
\arguments{
  \item{dframe}{Name of data.frame containing baseline X covariates.}
  \item{xvars}{List of names of X variable(s).}
  \item{method}{Hierarchical Clustering Method: "diana", "agnes" or "hclus".}
}
\details{The first step in an Unsupervised Propensity Scoring alalysis is always
  to hierarchically cluster patients in baseline X-covariate space.  UPShclus uses
  a Mahalabobis metric and clustering methods from the R "cluster" library for this
  key initial step.
}
\value{ An output list object of class UPShclus:
  \item{dframe}{Name of data.frame containing baseline X covariates.}
  \item{xvars}{List of names of X variable(s).}
  \item{method}{Hierarchical Clustering Method: "diana", "agnes" or "hclus".}
  \item{upshcl}{Hierarchical clustering object created by choice between three possible methods.}
}
\references{
  Kaufman L, Rousseeuw PJ.  (1990) \bold{Finding Groups in Data.  An Introduction to
  Cluster Analysis}.  New York: John Wiley and Sons.

  Kereiakes DJ, Obenchain RL, Barber BL, et al. (2000) Abciximab provides
  cost effective survival advantage in high volume interventional practice.
  \emph{Am Heart J} \bold{140}: 603-610.

  Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
  \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.

  Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.

  Rubin DB. (1980) Bias reduction using Mahalanobis metric matching.
  \emph{Biometrics} \bold{36}: 293-298.
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{UPSaccum}}, \code{\link{UPSnnltd}} and \code{\link{UPSgraph}}.}
\examples{
  data(lindner)
  UPSxvars <- c("stent", "height", "female", "diabetic", "acutemi", "ejecfrac", "ves1proc")
  UPSharch <- UPShclus(lindner, UPSxvars)
  plot(UPSharch)
}
\keyword{cluster}
\keyword{design}