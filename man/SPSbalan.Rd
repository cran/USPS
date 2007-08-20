\name{SPSbalan}
\alias{SPSbalan}
\title{Test for Within-Bin X-covariate Balance in Supervised Propensiy Scoring}
\description{Test for Conditional Independence of X-covariate Distributions from Treatment
  Selection within Given, Adjacent PS Bins.}
\usage{
  xcvobj <- SPSbalan(dframe, trtm, qbin, xvar, faclev=3)
}
\arguments{
  \item{dframe}{Name of augmented data.frame written to the appn="" argument of SPSlogit().}
  \item{trtm}{Name of the two-level treatment factor variable.}
  \item{qbin}{Name of variable containing bin numbers.}
  \item{xvar}{Name of one baseline covariate X variable used in the SPSlogit() PS model.}
  \item{faclev}{Maximum number of different numerical values an X-covariate can assume without
   automatically being converted into a "factor" variable; faclev=1 causes a binary indicator
   to be treated as a continuous variable determining a proportion.}
}
\details{The second step in Supervised Propensity Scoring analyses is to verify that baseline
  X-covariates have the same distribution, regardless of treatment, within each fitted PS bin.
}
\value{ An output list object of class SPSbalan.

 "contin"uous xvar => only the following four outputs...
  \item{aovdiff}{ANOVA output for marginal test.}
  \item{form2}{Formula for differences in X due to bins and to treatment nested within bins.}
  \item{bindiff}{ANOVA output for the nested within bin model.}
  \item{df3}{Output data.frame containing 3 variables: X-covariate, treatment and bin.}

 "factor" xvar => only the following four outputs... 
  \item{factab}{Marginal table of counts by X-factor level and treatment.}
  \item{tab}{Three-way table of counts by X-factor level, treatment and bin.}
  \item{cumchi}{Cumulative Chi-Square statistic for interaction in the three-way, nested table.}
  \item{cumdf}{Degrees of-Freedom for the Cumulative Chi-Squared.}
}
\references{
  Cochran WG. (1968) The effectiveness of adjustment by subclassification
  in removing bias in observational studies. \emph{Biometrics} \bold{24}:
  205--213.

  Obenchain RL. (2006) \bold{USPSinR.pdf}  ../R\_HOME/library/USPS  40 pages.

  Rosenbaum PR, Rubin RB. (1983) The Central Role of the Propensity Score
  in Observational Studies for Causal Effects. \emph{Biometrika} \bold{70}:
  41--55.

  Rosenbaum PR, Rubin DB. (1984) Reducing Bias in Observational Studies
  Using Subclassification on a Propensity Score. \emph{J Amer Stat Assoc}
  \bold{79}: 516--524.
}
\author{Bob Obenchain <softrx@iquest.net>}
\seealso{\code{\link{SPSlogit}}, \code{\link{SPSnbins}} and \code{\link{SPSoutco}}.}
\examples{
  data(lindner)
  PStreat <- abcix~stent+height+female+diabetic+acutemi+ejecfrac+ves1proc
  logtSPS <- SPSlogit(lindner, PStreat, PSfit, PSrnk, PSbin, appn="lindSPS")

  SPSbalvs <- SPSbalan(lindSPS, abcix, PSbin, ves1proc)
  SPSbalvs
  plot(SPSbalvs)
}
\keyword{htest}