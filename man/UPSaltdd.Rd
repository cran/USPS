\name{UPSaltdd}
\alias{UPSaltdd}
\title{Artificial Distribution of LTDs from Random Clusters}
\description{For a given number of clusters, UPSaltdd() characterizes the potentially biased
  distribution of "Local Treatment Differences" (LTDs) in a continuous outcome y-variable between
  two treatment groups due to Random Clusterings.  When the NNobj argument is not NA and specifies
  an existing UPSnnltd() object, UPSaltdd() also computes a smoothed CDF for the NN/LTD 
  distribution for direct comparison with the Artificial LTD distribution. 
}
\usage{UPSaltdd(dframe, trtm, yvar, faclev=3, scedas="homo", NNobj=NA, clus=50, reps=10, seed=12345)}
\arguments{
  \item{dframe}{Name of data.frame containing a treatment-factor and the outcome y-variable.}
  \item{trtm}{Name of treatment factor variable with two levels.}
  \item{yvar}{Name of continuous outcome variable.}
  \item{faclev}{Maximum number of different numerical values an outcome variable can assume
   without automatically being converted into a "factor" variable; faclev=1 causes a binary
   indicator to be treated as a continuous variable determining an average or proportion.}
  \item{scedas}{Scedasticity assumption: "homo" or "hete"}
  \item{NNobj}{Name of an existing UPSnnltd object or NA.}
  \item{clus}{Number of Random Clusters requested per Replication; ignored when NNobj is not NA.}
  \item{reps}{Number of overall Replications, each with the same number of requested clusters.}
  \item{seed}{Seed for Monte Carlo random number generator.}
}
\details{Multiple calls to UPSaltdd() for different UPSnnltd objects or different numbers of
  clusters are typically made after first invoking UPSgraph().
}
\value{ An output list object of class UPSaltdd:
  \item{dframe}{Name of data.frame containing X, t & Y variables.}
  \item{trtm}{Name of treatment factor variable.}
  \item{yvar}{Name of outcome Y variable.}
  \item{faclev}{Maximum number of different numerical values an outcome variable can assume without
   automatically being converted into a "factor" variable; faclev=1 causes a binary indicator to be
   treated as a continuous variable determining an average or proportion.}
  \item{scedas}{Scedasticity assumption: "homo" or "hete"}
  \item{NNobj}{Name of an existing UPSnnltd object or NA.}
  \item{clus}{Number of Random Clusters requested per Replication.}
  \item{reps}{Number of overall Replications, each with the same number of requested clusters.}
  \item{pats}{Number of patients with no NAs in their yvar outcome and trtm factor.}
  \item{seed}{Seed for Monte Carlo random number generator.}
  \item{altdd}{Matrix of LTDs and relative weights from artificial clusters.}
  \item{alxmin}{Minimum artificial LTD value.}
  \item{alxmax}{Maximum artificial LTD value.}
  \item{alymax}{Maximum weight among artificial LTDs.}  
  \item{altdcdf}{Vector of artificial LTD x-coordinates for smoothed CDF.}
  \item{qq}{Vector of equally spaced CDF values from 0.0 to 1.0.}
  \item{nnltdd}{Optional matrix of relevant NN/LTDs and relative weights.}
  \item{nnlxmin}{Optional minimum NN/LTD value.}
  \item{nnlxmax}{Optional maximum NN/LTD value.}
  \item{nnlymax}{Optional maximum weight among NN/LTDs.}  
  \item{nnltdcdf}{Optional vector of NN/LTD x-coordinates for smoothed CDF.}
  \item{nq}{Optional vector of equally spaced CDF values from 0.0 to 1.0.}
}
\references{
  Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
  \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.

  Obenchain RL. (2006) \bold{USPSinR.pdf}  ../R\_HOME/library/USPS  41 pages.

  Rosenbaum PR, Rubin RB. (1983) The Central Role of the Propensity Score
  in Observational Studies for Causal Effects. \emph{Biometrika} \bold{70}:
  41--55.

  Rubin DB. (1980) Bias reduction using Mahalanobis metric matching.
  \emph{Biometrics} \bold{36}: 293--298.
}
\author{Bob Obenchain <softrx@iquest.net>}
\seealso{\code{\link{UPSnnltd}}, \code{\link{UPSaccum}} and \code{\link{UPSgraph}}.}
\examples{
    data(lindner)
    abcdf <- UPSaltdd(lindner, abcix, lifepres, faclev=1)
    abcdf
    plot(abcdf)
}
\keyword{nonparametric}