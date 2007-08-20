\name{plot.UPSnnltd}
\alias{plot.UPSnnltd}
\title{Display plots of the NN/LTD Distribution in Unsupervised Propensiy Scoring}
\description{Make a Snow Ball ("snob"), Probability Density ("dens") and/or Cumulative Distribution
  Function ("cdf") plot of the Local Treatment Difference (LTD) Distribution across all Informative
  Clusters within a requested total number clusters in patient baseline X-covariate space.
}
\usage{
  \method{plot}{UPSnnltd}(x, pballs = TRUE, nnplot = "snob", nnalpha = 1.4, \dots)
}
\arguments{
  \item{x}{Object of class UPSnnltd.}
  \item{pballs}{Logical; Display cluster-size weighted (red) and inverse-variance weighted (green)
     "Power Balls" to represent across cluster LTD averages within the "Snowball" plot.}
  \item{nnplot}{String; Display "all" plots together, "seq" = one after the other, "snob" = only Snow Balls,
     "dens" = only fitted gss probability density, or "cdf" = only the gss cumulative distribution.}
  \item{nnalpha}{Numeric; alpha argument to ssden() function.}
  \item{\dots}{Optional argument(s) to pass on to plot().}
}
\value{ NULL }
\details{Plots of UPSnnltd objects describe the LTD Distribution across all informative clusters.  The
  generalized smoothing spline (gss) displays corresponding to the nnplot = "dens" and "cdf" options
  may fail to exist when the number of informative clusters is less than 60.
}
\references{
  Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
  \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.

  Obenchain RL. (2006) \bold{USPSinR.pdf}  ../R\_HOME/library/USPS  40 pages.
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
\keyword{methods}
\keyword{hplot}
