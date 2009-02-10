\name{plot.SPSloess}
\alias{plot.SPSloess}
\title{Display LOESS Smooth of Outcome by Treatment in Supervised Propensiy Scoring}
\description{Express Expected Outcome by Treatment as LOESS Smooths of Fitted Propensity Scores.
}
\usage{
  \method{plot}{SPSloess}(x, tcol="blue", ucol="red", dcol="green3", \dots)
}
\arguments{
  \item{x}{output list object of class SPSloess.}
  \item{tcol}{optional; quoted name of color for treated patient smooth.}
  \item{ucol}{optional; quoted name of color for untreated patient smooth.}
  \item{dcol}{optional: quoted name of color for combined patient density.}
  \item{\dots}{optional; argument(s) passed on to plot().}
}
\details{Plots of SPSloess objects display each patient's propensity score versus his/her observed
 (continuous) outcome.  Patients receiving the "standard" treatment (trtm=0) are represented by
 cyan circles, while patients receiving the "new" treatment (trtm=1) are represented by magenta
 triangles.  The smooth fits of outcome to propensity score within treatment cohorts are show as
 cyan (trtm=0) and magenta (trtm=1) curves, respectively, superimposed upon the scatter.

 Because smooth fits can be difficult to see when the scatters contain many points, a second
 plot rescaled to show only the two smooth (lowess or spline) fits, again using cyan (trtm=0)
 and magenta (trtm=1) curves.  For details, see the returned lofit data frame.

 Finally, a third plot shows total patient frequencies (black circles) within a 100-cell
 histogram along the propensity score axis as well as the corresponding density()
 smooth in red.  For details, see the returned logrid data frame.

 Winsorizing Cost data: PSframe$TRIMBILL <- pmin( PSframe$cardbill, 50000)

 The fam="symmetric" default option of SPSloess tends to be fairly robust to outlying outcomes, at least when
 the loess span is wide enough.  Thus reducing (Winsorizing) outlying cardbill values to \$50K (as illustrated
 above) should have little effect on a fitted loess smooth with an appropriate span.  Looking for the effects
 of Winsorizing on SPSloess() or SPSsmoot() constitutes "sensitivity analysis."

 The original lowess() function of Cleveland and Devlin (1988) could be used here because only one X variable
 (namely, fitted propensity score) is involved, but I choose loess() instead to give users flexibility to choose
 between fam="gaussian" and fam="symmetric" option, which provides some resistance to outlying outcome values.

 SPSloess() fits can tend to look rather "rough" compared to SPSsmoot() fits.  Cubic spline smoothing
 appears to give answers that are interpretable as smoothed mean values for highly skewed distributions.  Loess
 smoothing, at least when fam="symmetric," tends to give answers more easily interpretable as modes or medians
 of highly skewed distributions.  This median versus mean analogy may help explain why the weighted average
 signed treatment differences from SPSloess() tend to seem more precise than those from SPSsmoot() for highly
 skewed distributions.
}
\value{ NULL }
\references{
  Cleveland WS, Devlin SJ. (1988) Locally-weighted regression: an
  approach to regression analysis by local fitting. \emph{J Amer Stat Assoc}
  \bold{83}: 596-610.

  Cleveland WS, Grosse E, Shyu WM. (1992) Local regression models. Chapter 8 of
  \bold{Statistical Models in S} eds Chambers JM and Hastie TJ. \emph{Wadsworth & Brooks/Cole}.

  Obenchain RL. (2009) \bold{USPSinR.pdf}  ../R\_HOME/library/USPS  40 pages.

  Ripley BD, loess() based on the 'cloess' package of Cleveland, Grosse and Shyu.
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{SPSlogit}}, \code{\link{SPSsmoot}} and \code{\link{SPSoutco}}.}
\examples{
  data(lindner)
  PStreat <- abcix~stent+height+female+diabetic+acutemi+ejecfrac+ves1proc
  logtSPS <- SPSlogit(lindner, PStreat, PSfit, PSrnk, PSbin, appn="lindSPS")

  SPScbls5 <- SPSloess(lindSPS, abcix, PSfit, cardbill, span=.5)
  SPScbls5 
  plot(SPScbls5)   
}
\keyword{methods}
\keyword{hplot} 