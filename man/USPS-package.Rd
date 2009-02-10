\name{USPS-package}
\alias{USPS-package}
\docType{package}
\title{Unsupervised and Supervised Propensity Scoring adjustments for Bias and Confounding}
\description{Identify and display the distribution of Local Treatment Differences (LTDs)
  and Local Average Treatment Effects (LATEs) across Clusters of patients chosen so as
  to be relatively well matched within Clusters defined by patient baseline X-covariates.
}
\details{
  \tabular{ll}{
  Package: \tab USPS\cr
  Type: \tab Package\cr
  Version: \tab 1.2-0\cr
  Date: \tab 2009-02-10\cr
  License: \tab GNU GENERAL PUBLIC LICENSE, Version 2, June 1991\cr
  }

  SUPERVISED OUTCOME BINNING AND SMOOTHING WITH ESTIMATED PROPENSITY SCORES:
  
  Once one has fitted a somewhat smooth curve through scatters of observed outcomes, Y,
  versus fitted propensity scores, X, for the patients in each of the two treatment groups,
  one can consider the question: "Over the range where both smooth curves are defined (i.e. their
  common support), what is the (weighted) average signed difference between these two curves?"

  UNSUPERVISED NEAREST NEIGHBORS / LOCAL TREATMENT DIFFERENCES:
  
  Multiple calls to UPSnnltd(n) for varying numbers of clusters, n, are typically made
  after first invoking UPShclus() to hierarchically cluster patients in X-space and then
  invoking UPSaccum() to specify a Y outcome variable and a two-level treatment factor t.
  UPSnnltd(n) then determines the LTD Distribution corresponding to n clusters and,
  optionally, displays this distribution in a "Snowball" plot.
  
  UNSUPERVISED INSTRUMENTAL VARIABLES / LOCAL AVERAGE TREATMENT EFFECTS:
    
  Multiple calls to UPSivadj(n) for varying numbers of clusters, n, yield alternative
  linear smoothes of LATE estimates plotted versus within cluster propensity score
  (observed treatment fraction) percentages.
}
\author{Bob Obenchain <wizbob@att.net>}
\references{
  Green PJ, Silverman BW. (1994) \bold{Nonparametric Regression and Generalized Linear Models:
  A Roughness Penalty Approach}. \emph{Chapman and Hall}.
  
  McClellan M, McNeil BJ, Newhouse JP. (1994) Does More Intensive Treatment of
  Myocardial Infarction in the Elderly Reduce Mortality?: Analysis Using Instrumental
  Variables. \emph{JAMA} \bold{272}: 859-866.  
  
  Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
  \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.
 
  Obenchain RL.  \bold{USPSinR.pdf}  ../R\_HOME/library/USPS  2009; 40 pages.
  
  Rosenbaum PR, Rubin RB. (1983) The Central Role of the Propensity Score
  in Observational Studies for Causal Effects. \emph{Biometrika} \bold{70}:
  41-55.

  Rubin DB. (1980) Bias reduction using Mahalanobis metric matching.
  \emph{Biometrics} \bold{36}: 293-298.
}
\keyword{package}
\examples{
  demo(abcix)
}