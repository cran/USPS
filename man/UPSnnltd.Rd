\name{UPSnnltd}
\alias{UPSnnltd}
\title{Nearest Neighbor Distribution of LTDs in Unsupervised Propensiy Scoring}
\description{For a given number of patient clusters in baseline X-covariate space, UPSnnltd()
  characterizes the distribution of Nearest Neighbor "Local Treatemnt Differences" (LTDs) on
  a specified Y-outcome variable.}
\usage{
  nnobj <- UPSnnltd(numclust)
}
\arguments{
  \item{numclust}{Number of clusters in baseline X-covariate space.}
}
\details{Multiple calls to UPSnnltd(n) for varying numbers of clusters, n, are typically made
  after first invoking UPShclus() to hierarchically cluster patients in X-space and then
  invoking UPSaccum() to specify a Y outcome variable and a two-level treatment factor t.
  UPSnnltd(n) then determines the LTD Distribution corresponding to n clusters and,
  optionally, displays this distribution in a "Snowball" plot.
}
\value{ An output list object of class UPSnnltd:
  \item{hiclus}{Name of clustering object created by UPShclus().}
  \item{dframe}{Name of data.frame containing X, t & Y variables.}
  \item{trtm}{Name of treatment factor variable.}
  \item{yvar}{Name of outcome Y variable.}
  \item{numclust}{Number of clusters requested.}
  \item{actclust}{Number of clusters actually produced.}
  \item{scedas}{Scedasticity assumption: "homo" or "hete"}
  \item{PStdif}{Character string describing the treatment difference.}
  \item{nnhbindf}{Vector containing cluster number for each patient. }
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
  \item{faclev}{Maximum number of different numerical values an outcome variable can assume without
   automatically being converted into a "factor" variable; faclev=1 causes a binary indicator to be
   treated as a continuous variable determining an average or proportion.}
  \item{youtype}{"contin"uous => only next eight outputs; "factor" => only last three outputs.}
  \item{aovdiff}{ANOVA summary for treatment main effect only.}
  \item{form2}{Formula for outcome differences due to bins and to treatment nested within bins.}
  \item{bindiff}{ANOVA summary for treatment nested within cluster.}
  \item{sig2}{Estimate of error mean square in nested model.}
  \item{pbindif}{Unadjusted treatment difference by cluster.}
  \item{pbinsde}{Standard error of the unadjusted difference by cluster.}
  \item{pbinsiz}{Cluster radii measure: square root of total number of patients.}
  \item{symsiz}{Symbol size of largest possible Snowball in a UPSnnltd() plot with 1 cluster.}
  \item{factab}{Marginal table of counts by Y-factor level and treatment.}
  \item{cumchi}{Cumulative Chi-Square statistic for interaction in the three-way, nested table.}
  \item{cumdf}{Degrees of-Freedom for the Cumulative Chi-Squared.}
}
\references{
  Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
  \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.

  Obenchain RL. (2006) \bold{USPSinR.pdf}  ../R\_HOME/library/USPS  40 pages.

  Rosenbaum PR, Rubin RB. (1983) The Central Role of the Propensity Score
  in Observational Studies for Causal Effects. \emph{Biometrika} \bold{70}:
  41--55.

  Rubin DB. (1980) Bias reduction using Mahalanobis metric matching.
  \emph{Biometrics} \bold{36}: 293--298.
}
\author{Bob Obenchain <softrx@iquest.net>}
\seealso{\code{\link{UPSivadj}}, \code{\link{UPSaccum}} and \code{\link{UPSgraph}}.}
\examples{
    data(lindner)
    UPSxvars <- c("stent", "height", "female", "diabetic", "acutemi", "ejecfrac", "ves1proc")
    UPSharch <- UPShclus(lindner, UPSxvars)
    UPSaccum(UPSharch, lindner, abcix, lifepres, faclev=1, scedas="homo", accobj="ABClife")
    
    lif070nn <- UPSnnltd(70)
    lif070nn
    plot(lif070nn)
}
\keyword{nonparametric}