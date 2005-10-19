\name{plot.UPSnnltd}
\alias{plot.UPSnnltd}
\title{Display a Snowball plot of the Nearest Neighbor Distribution of LTDs in Unsupervised Propensiy Scoring}
\description{Make a "Snowball" plot of the Local Treatment Difference (LTD) distribution for a
  given number of patient clusters in baseline X-covariate space.
}
\usage{
  \method{plot}{UPSnnltd}(x, pballs = TRUE, \dots)
}
\arguments{
  \item{x}{Object of class UPSnnltd.}
  \item{pballs}{Logical; Display cluster-size weighted (red) and inverse-variance weighted (green)
     "Power Balls" to represent across cluster LTD averages within the "Snowball" plot.}
  \item{\dots}{Optional argument(s) to pass on to plot().}
}
\value{ NULL }
\references{
  Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
  \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.

  Obenchain RL. (2005) \bold{USPSinR.pdf}  ../R\_HOME/library/USPS  40 pages.
}
\author{Bob Obenchain <sunsetstats@earthlink.net>}
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
\keyword{methods}
\keyword{hplot}
