\name{plot.SPSsmoot}
\alias{plot.SPSsmoot}
\title{Display Spline Smooth of Outcome by Treatment in Supervised Propensiy Scoring}
\description{Express Expected Outcome by Treatment as Spline Functions of Fitted
 Propensity Scores.
}
\usage{
  \method{plot}{SPSsmoot}(x, tcol = "blue", ucol = "red", dcol = "green3", \dots)
}
\arguments{
  \item{x}{output list object of class SPSsmoot.}
  \item{tcol}{optional; quoted name of color for treated patient smooth.}
  \item{ucol}{optional; quoted name of color for untreated patient smooth.}
  \item{dcol}{optional; quoted name of color for combined patient density.}
  \item{\dots}{optional; argument(s) passed on to plot().}
}
\details{Plots of SPSsmoot objects display each patient's propensity score versus his/her observed
 (continuous) outcome.  Patients receiving the "standard" treatment (trtm=0) are represented by
 cyan circles, while patients receiving the "new" treatment (trtm=1) are represented by magenta
 triangles.  The smooth fits of outcome to propensity score within treatment cohorts are show as
 cyan (trtm=0) and magenta (trtm=1) curves, respectively, superimposed upon the scatter.

 Because smooth fits can be difficult to see when the scatters contain many points, a second
 plot rescaled to show only the two smooth (lowess or spline) fits, again using cyan (trtm=0)
 and magenta (trtm=1) curves.  For details, see the returned ssfit data frame.

 Finally, a third plot shows total patient frequencies (black circles) within a 100-cell
 histogram along the propensity score axis as well as the corresponding density()
 smooth in red.  For details, see the returned ssgrid data frame.

 SPSloess() fits can tend to look rather "rough" compared to SPSsmoot() fits.  Cubic spline smoothing
 appears to give answers that are interpretable as smoothed mean values for highly skewed distributions.  Loess
 smoothing, at least when fam="symmetric," tends to give answers more easily interpretable as modes or medians
 of highly skewed distributions.  This median versus mean analogy may help explain why the weighted average
 signed treatment differences from SPSloess() tend to seem more precise than those from SPSsmoot() for highly
 skewed distributions.
}
\value{ NULL }
\references{
  Chambers JM, Hastie T. (1992) \bold{Statistical Models in S} \emph{Wadsworth & Brooks/Cole}.

  Green PJ, Silverman BW. (1994) \bold{Nonparametric Regression and Generalized Linear Models:
  A Roughness Penalty Approach}. \emph{Chapman and Hall}.

  Hastie TJ, Tibshirani RJ. (1990) \bold{Generalized Additive Models}. \emph{Chapman and Hall}.

  Obenchain RL. (2005) \bold{USPSinR.pdf}  ../R\_HOME/library/USPS  40 pages.

  Sheather SJ, Jones MC. (1991) A reliable data-based bandwidth selection method for kernel
  density estimation. \emph{J Roy Statist Soc B} \bold{53}: 683--690.

  R implementation 0f smooth.spline() by Ripley BD and Maechler M. ('spar/lambda', etc).
}
\author{Bob Obenchain <sunsetstats@earthlink.net>}
\seealso{\code{\link{SPSloess}}, \code{\link{SPSbalan}} and \code{\link{SPSoutco}}.}
\examples{
  data(lindner)
  PStreat <- abcix~stent+height+female+diabetic+acutemi+ejecfrac+ves1proc
  logtSPS <- SPSlogit(lindner, PStreat, PSfit, PSrnk, PSbin, appn="lindSPS")

  SPScbss7 <- SPSsmoot(lindSPS, abcix, PSfit, cardbill, df=7)
  SPScbss7
  plot(SPScbss7)
}
\keyword{methods}
\keyword{hplot}

