\name{SPSnbins}
\alias{SPSnbins}
\title{Change the Number of Bins in Supervised Propensiy Scoring}
\description{Change the Number of Bins in Supervised Propensiy Scoring}
\usage{dframe2 <- SPSnbins(dframe1, prnk, qbin, bins=8)}
\arguments{
  \item{dframe1}{Name of data.frame of the form output by SPSlogit().}
  \item{prnk}{Name of PS tied-rank variable from previous call to SPSlogit().}
  \item{qbin}{Name of variable to contain the re-assigned bin number for each patient.}
  \item{bins}{Number of PS bins desired.}}
\details{Part or all of the first phase of Supervised Propensity Scoring will need to be redone
  if SPSbalan() detects dependence of within-bin X-covariate distributions upon treatment choice.
  Use SPSnbins() to change (increase) the number of adjacent PS bins.  If this does not achieve
  balance, invoke SPSlogit() again to modify the form of your PS logistic model, typically by
  adding interaction and/or curvature terms in continuous X-covariates.}
\value{ An output data.frame with new variables inserted:
  \item{dframe2}{Modified version of the data.frame specified as the first argument to SPSnbins().}
}
\references{
  Cochran WG. (1968) The effectiveness of adjustment by subclassification
  in removing bias in observational studies. \emph{Biometrics} \bold{24}:
  205--213.

  Obenchain RL. (2006) \bold{USPSinR.pdf}  ../R\_HOME/library/USPS  40 pages.

  Rosenbaum PR, Rubin DB. (1984) Reducing Bias in Observational Studies
  Using Subclassification on a Propensity Score. \emph{J Amer Stat Assoc}
  \bold{79}: 516--524.
}
\author{Bob Obenchain <softrx@iquest.net>}
\seealso{\code{\link{SPSlogit}}, \code{\link{SPSbalan}} and \code{\link{SPSoutco}}.}
\examples{
  data(lindner)
  PStreat <- abcix~stent+height+female+diabetic+acutemi+ejecfrac+ves1proc
  logtSPS <- SPSlogit(lindner, PStreat, PSfit, PSrnk, PSbin, appn="lindSPS")
  logtSPS

  # If imbalance within the SPSlogit() default of bins=5 was detected, then ...

  lindSPS <- SPSnbins(lindSPS, PSrnk, PSbin6, bins=6)
  table(lindSPS$PSbin,lindSPS$PSbin6)
}
\keyword{design}