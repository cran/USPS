\name{SPSoutco}
\alias{SPSoutco}
\title{Examine Treatment Differences on an Outcome Measure in Supervised Propensiy Scoring}
\description{Examine Within-Bin Treatment Differences on an Outcome Measure and Average these
  Differences across Bins.}
\usage{
  outobj <- SPSoutco(dframe, trtm, qbin, yvar, faclev=3)
}
\arguments{
  \item{dframe}{Name of augmented data.frame written to the appn="" argument of SPSlogit().}
  \item{trtm}{Name of treatment factor variable.}
  \item{qbin}{Name of variable containing the PS bin number for each patient.}
  \item{yvar}{Name of an outcome Y variable.}
  \item{faclev}{Maximum number of different numerical values an X-covariate can assume without
   automatically being converted into a "factor" variable; faclev=1 causes a binary indicator to
   be treated as a continuous variable determining an average or proportion.}
}
\details{Once the second phase of Supervised Propensity Scoring confirms, using SPSbalan(), that
  X-covariate Distributions have been Balanced Within-Bins, the third phase can start: Examining
  Within-Bin Outcome Difference due to Treatment and Averaging these Differences across Bins.
  Graphical displays of SPSoutco() results feature R barplot() invocations.
}
\value{ An output list object of class SPSoutco:
  \item{dframe}{Name of augmented data.frame written to the appn="" argument of SPSlogit().}
  \item{trtm}{Name of the two-level treatment factor variable.}
  \item{yvar}{Name of an outcome Y variable.}
  \item{bins}{Number of variable containing bin numbers.}
  \item{PStdif}{Character string describing the treatment difference.}
  \item{rawmean}{Unadjusted outcome mean by treatment group.}
  \item{rawvars}{Unadjusted outcome variance by treatment group.}
  \item{rawfreq}{Number of patients by treatment group.}
  \item{ratdif}{Unadjusted mean outcome difference between treatments.}
  \item{ratsde}{Standard error of unadjusted mean treatment difference.}
  \item{binmean}{Unadjusted mean outcome by cluster and treatment.}
  \item{binvars}{Unadjusted variance by cluster and treatment.}
  \item{binfreq}{Number of patients by bin and treatment.}
  \item{awbdif}{Across cluster average difference with cluster size weights.}
  \item{awbsde}{Standard error of awbdif.}
  \item{wwbdif}{Across cluster average difference, inverse variance weights.}
  \item{wwbsde}{Standard error of wwbdif.}
  \item{form}{Formula for overall, marginal treatment difference on X-covariate.}
  \item{faclev}{Maximum number of different numerical values an X-covariate can assume without
   automatically being converted into a "factor" variable; faclev=1 causes a binary indicator to be
   treated as a continuous variable determining an average or proportion.}
  \item{youtype}{"contin"uous => only next six outputs; "factor" => only last four outputs.}
  \item{aovdiff}{ANOVA output for marginal test.}
  \item{form2}{Formula for differences in X due to bins and to treatment nested within bins.}
  \item{bindiff}{ANOVA summary for treatment nested within bin.}
  \item{pbindif}{Unadjusted treatment difference by cluster.}
  \item{pbinsde}{Standard error of the unadjusted difference by cluster.}
  \item{pbinsiz}{Cluster radii measure: square root of total number of patients.}
  \item{factab}{Marginal table of counts by Y-factor level and treatment.}
  \item{tab}{Three-way table of counts by Y-factor level, treatment and bin.}
  \item{cumchi}{Cumulative Chi-Square statistic for interaction in the three-way, nested table.}
  \item{cumdf}{Degrees of-Freedom for the Cumulative Chi-Squared.}
}
\references{
  Cochran WG. (1968) The effectiveness of adjustment by subclassification
  in removing bias in observational studies. \emph{Biometrics} \bold{24}:
  205--213.

  Obenchain RL. (2005) \bold{USPSinR.pdf}  ../R\_HOME/library/USPS  40 pages.

  Rosenbaum PR, Rubin RB. (1983) The Central Role of the Propensity Score
  in Observational Studies for Causal Effects. \emph{Biometrika} \bold{70}:
  41--55.

  Rosenbaum PR, Rubin DB. (1984) Reducing Bias in Observational Studies
  Using Subclassification on a Propensity Score. \emph{J Amer Stat Assoc}
  \bold{79}: 516--524.
}
\author{Bob Obenchain <sunsetstats@earthlink.net>}
\seealso{\code{\link{SPSlogit}}, \code{\link{SPSbalan}} and \code{\link{SPSnbins}}.}
\examples{
  data(lindner)
  PStreat <- abcix~stent+height+female+diabetic+acutemi+ejecfrac+ves1proc
  logtSPS <- SPSlogit(lindner, PStreat, PSfit, PSrnk, PSbin, appn="lindSPS")

  SPSlifeo <- SPSoutco(lindSPS, abcix, PSbin, lifepres, faclev=1)
  SPSlifeo
  plot(SPSlifeo)    
}
\keyword{nonparametric}
\keyword{hplot}
