\name{UPSgraph}
\alias{UPSgraph}
\title{Display Sensitivity Analysis Graphic in Unsupervised Propensiy Scoring}
\description{Plot summary of results from multiple calls to UPSnnltd() and/or UPSivadj() after an
  initial setup call to UPSaccum().  The UPSgraph() plot displays any sensitivity of the LTD and
  LATE Distributions to choice of Number of Clusters in X-space.}
\usage{
  UPSgraph(nncol = "red", nwcol = "green3", ivcol = "blue", \dots)
}
\arguments{
  \item{nncol}{optional; string specifying color for display of the Mean of the LTD
   distribution when weighted by cluster size from any calls to UPSnnltd().}
  \item{nwcol}{optional; string specifying color for display of the Mean of the LTD
   distribution when weighted inversely proportional to variance from any calls to UPSnnltd().}
  \item{ivcol}{optional; string specifying color for display of the Difference in LATE
   prediction, at PS = 100\% minus that at PS = 0\%, from any calls to UPSivadj().}
  \item{\dots}{Optional parameter(s) passed on to plot().}
}
\details{The third phase of Unsupervised Propensity Scoring is a graphical Sensitivity
  Analysis that depicts how the Overall Means of the LTD and LATE distributions change with
  the number of clusters.
}
\value{ NULL }
\references{
  Kaufman L, Rousseeuw PJ. (1990) \bold{Finding Groups in Data.  An Introduction to
  Cluster Analysis}.  \emph{New York: John Wiley and Sons}.

  Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
  \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.

  Obenchain RL. (2005) \bold{USPSinR.pdf}  ../R\_HOME/library/USPS  40 pages.

  Rubin DB. (1980) Bias reduction using Mahalanobis metric matching.
  \emph{Biometrics} \bold{36}: 293--298.
}
\author{Bob Obenchain <sunsetstats@earthlink.net>}
\seealso{\code{\link{UPSnnltd}}, \code{\link{UPSivadj}} and \code{\link{UPSaccum}}.}
\examples{
  data(lindner)
  UPSxvars <- c("stent", "height", "female", "diabetic", "acutemi", "ejecfrac", "ves1proc")
  UPSharch <- UPShclus(lindner, UPSxvars)
  UPSaccum(UPSharch, lindner, abcix, lifepres, faclev=1, scedas="homo", accobj="ABClife")
  lif001nn <- UPSnnltd(1)
  lif020nn <- UPSnnltd(20)
  lif070nn <- UPSnnltd(70)
  lif120nn <- UPSnnltd(120)

  UPSgraph()

  ABClife
}
\keyword{hplot}